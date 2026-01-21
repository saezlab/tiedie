import os
import re
import sys
import math
import atexit
import random
import operator
import warnings

import networkx as nx


# Track linker cutoff errors to avoid flooding warnings
_linker_cutoff_error_count = 0
_linker_cutoff_error_shown = False


def _report_linker_cutoff_errors():
    """Report total linker cutoff errors at exit."""
    if _linker_cutoff_error_count > 1:
        warnings.warn(
            f'find_linker_cutoff failed {_linker_cutoff_error_count} '
            'times total during execution',
            stacklevel=1,
        )


atexit.register(_report_linker_cutoff_errors)


def parse_heats(file, network_nodes=None):
    """Parse input heats file in form:
        <gene> <heat> <perturbation/activity sign (+/-)>

    Returns:
        - Two hashes: one indexing by gene and storing the input heats, and one storing the input signs
    """

    heats = {}
    signs = {}
    fh = None
    try:
        fh = open(file)
    except OSError as e:
        raise Exception("Error: can't open file: " + file) from e

    lineno = 1
    for line in fh:
        parts = line.rstrip().split('\t')
        if len(parts) > 2:
            prot, heat, sign = line.rstrip().split('\t')

            # provide a warning if node not in the network
            if network_nodes and prot not in network_nodes:
                sys.stderr.write(
                    'Warning: input heat node '
                    + prot
                    + ' not in the network and will be ignored...\n'
                )
                continue

            # input validation for heat values
            try:
                heats[prot] = float(heat)
            except ValueError as e:
                raise Exception(
                    'Error: non float heat value on line '
                    + str(lineno)
                    + ' gene '
                    + prot
                ) from e

            # input validation for input signs
            if sign != '+' and sign != '-':
                raise Exception(
                    'Error: invalid value for heat sign on line '
                    + str(lineno)
                    + sign
                )

            signs[prot] = sign
        else:
            heats[parts[0]] = float(parts[1])

        lineno += 1

    fh.close()
    return (heats, signs)


def edgelist_to_nodes(list):
    """Input:
        A list of edges in (source, interaction, target) string form.

    Returns:
        A set object of nodes in the input network

    >>> edgelist_to_nodes([("A","i>","B"),("B","-a>","C")])
    set(['A', 'C', 'B'])
    """

    nodes = set()
    for source, i, target in list:
        nodes.add(source)
        nodes.add(target)

    return nodes


def classify_interaction(i):
    """Returns the edge activation type (-1,0,1), and the textual description

    >>> classify_interaction("component>")
    (0, 'component')
    >>> classify_interaction("-a>")
    (1, 'a')
    >>> classify_interaction("-t>")
    (1, 't')
    >>> classify_interaction("-t|")
    (-1, 't')
    >>> classify_interaction("-a|")
    (-1, 'a')
    >>> classify_interaction("HPRD>")
    (1, 'INTERACTS')
    >>> classify_interaction("REWIRED>")
    (1, 'REWIRED')
    """
    componentRE = re.compile('^-?component>$')
    activatingRE = re.compile(r'^-?(\S)>$')
    inactivatingRE = re.compile(r'^-?(\S)\|$')
    rewiredAC = re.compile('^-?REWIRED>$')
    rewiredIN = re.compile(r'^-?REWIRED\|$')

    if componentRE.match(i):
        return (0, 'component')
    elif activatingRE.match(i):
        type = activatingRE.match(i)
        return (1, type.group(1))
    elif inactivatingRE.match(i):
        type = inactivatingRE.match(i)
        return (-1, type.group(1))
    elif rewiredAC.match(i):
        type = 'REWIRED'
        return (1, type)
    elif rewiredIN.match(i):
        type = 'REWIRED'
        return (-1, type)
    else:
        # default to activating links for HPRD or other protein
        # component links. These are bi-directional
        return (1, 'INTERACTS')


def get_out_degrees(network):
    """Get the out-degree of each node in the network

    Input:
        network:
            { [source]: (interaction, target) }

    Returns:
        a hash of node out-degrees

    >>> network = {}
    >>> network['S1'] = set()
    >>> network['S2'] = set()
    >>> network['T1'] = set()
    >>> network['S1'].add(('a>','T1'))
    >>> network['S2'].add(('a>','T2'))
    >>> network['S2'].add(('a|','T3'))
    >>> network['T1'].add(('t|','T2'))
    >>> get_out_degrees(network)
    {'S2': 2, 'S1': 1, 'T2': 0, 'T3': 0, 'T1': 1}

    """
    outDegrees = {}
    for s in network:
        outDegrees[s] = len(network[s])
        for i, t in network[s]:
            if t not in outDegrees:
                outDegrees[t] = 0

    return outDegrees


def edges_to_degrees(edges):
    """Takes simple edges in (source, target) format, and returns a hash of the
    total degree of each node.

    >>> edges_to_degrees([("A","B"),("B","C")])
    {'A': 1, 'C': 1, 'B': 2}
    """

    nodes = {}
    for s, t in edges:
        if s not in nodes:
            nodes[s] = {}
        if t not in nodes:
            nodes[t] = {}

        nodes[s][t] = 1
        nodes[t][s] = 1

    sizes = {}
    for n in nodes:
        sizes[n] = len(nodes[n])

    return sizes


def is_rewired(i):
    """Check if interaction type indicates a rewired edge."""
    rewiredRE = re.compile('.*REWIRED.*')
    rewiredComponentRE = re.compile(r'.*\-component.*')

    if rewiredRE.match(i):
        return True
    elif rewiredComponentRE.match(i):
        return True

    return False


# do a depth-first search by following directional links
# until we hit another source
# find edges
def search_dfs(
    source,
    action,
    discovered,
    linker_nodes,
    target_set,
    net,
    gene_states,
    transcriptional_signs,
    depth,
    truePaths,
    falsePaths,
    falsePathStatus,
):
    """Perform a depth-first search by following directional links
    until we hit another source. Validate link interactions along the way.
    Recursive calls.

    Input:
        source: source node by name
        action: +1/-1 binary action
        discovered: store validated 'discovered' paths
        linker_nodes: build the list of linker nodes as we recurse through the function. Add them to validation
        list if they lead to a known target
        net: network in hash-format {'source':(interaction, target), ...}
        gene_states: hash of 'action' states for each gene in the network
        transcriptional_signs: equivalent to 'gene_states' for transcriptionally active nodes
        depth: level of recursion (stop if it hits zero)
        ...additional: counts for real/false paths if using the REWIRED link test

    Returns:
        None
    """

    if depth == 0:
        return

    if source not in net:
        return

    for interaction, target in net[source]:
        # if we arrived here through a bad link, any continued path is counted as false
        pathStatus_ThisTarget = falsePathStatus
        (i_type, post_t_type) = classify_interaction(interaction)
        if is_rewired(interaction):
            pathStatus_ThisTarget = True

        # don't follow component links
        if i_type == 0:
            continue

        # activating nodes keep the signal, inactivating nodes reverse the signal
        action_this_target = None
        if i_type == 1:
            action_this_target = action
        elif i_type == -1:
            action_this_target = -action

        # for transcriptional states: the expression activity is what we want to measure
        this_state = None
        if target in gene_states:
            this_state = gene_states[target]
        # if post_t_type == "t":
        # this depends on weather we monitor the activities of downstream genes, or just the transcription
        # leave commented out for the former
        #   if target not in transcriptional_signs:
        #       continue
        #   this_state = transcriptional_signs[target]

        # we hit a target that has a matching action/signal from the original source
        if (
            (target in gene_states)
            and (target in target_set)
            and (action_this_target == this_state)
        ):
            for s, i, t in linker_nodes:
                discovered.add((s, i, t))
            discovered.add((source, interaction, target))
            linker_nodes = set()
            new_linkers = set()
            # and keep going

            # add this to our TP or FP score, depending on the path
            if pathStatus_ThisTarget:
                falsePaths.append(target)
            else:
                truePaths.append(target)

        # search the target, but with any previous linkers
        else:
            new_linkers = set()
            new_linkers.add((source, interaction, target))
            new_linkers = new_linkers.union(linker_nodes)

        # if we come from a transcriptionally activating link, this cuts the cycle. Gene must
        # be upregulated
        if post_t_type == 't':
            continue

        # add this link and keep searching from the target
        search_dfs(
            target,
            action_this_target,
            discovered,
            new_linkers,
            target_set,
            net,
            gene_states,
            transcriptional_signs,
            depth - 1,
            truePaths,
            falsePaths,
            pathStatus_ThisTarget,
        )


def classify_state(up_signs, down_signs):
    """Build a hash of putative effects of perturbations,
    and inferred transcription activity.

    >>> classify_state({'A':"+",'B':"+"}, {'B':"-",'C':"-"})
    ({'A': 1, 'C': -1, 'B': 1}, {'C': -1, 'B': -1})
    """

    c = {}
    t_states = {}
    # The order matters here:
    for gene in down_signs:
        if down_signs[gene] == '+':
            c[gene] = 1
            t_states[gene] = 1
        else:
            c[gene] = -1
            t_states[gene] = -1

    # The order matters here:
    for gene in up_signs:
        if up_signs[gene] == '+':
            c[gene] = 1
        else:
            c[gene] = -1

    return (c, t_states)


# build an index, source to targets fro the directed graph
def parse_net(network):
    """Build a directed network from a .sif file.

    Inputs:
        A network in .sif format, tab-separated (<source> <interaction> <target>)

    Returns:
        A network in hash key format, i.e. convert two lines of a file:
            <source>    <interaction1>  <target1>
            <source>    <interaction2>  <target2>
        To:
            {'source': set( (interaction, target1), (interaction, target2) )
    """
    net = {}
    for line in open(network):
        parts = line.rstrip().split('\t')
        source = parts[0]
        interaction = parts[1]
        target = parts[2]

        if source not in net:
            net[source] = set()

        net[source].add((interaction, target))

    return net


def find_linker_cutoff(
    source_set, target_set, up_heat_diffused, down_heat_diffused, size
):
    """For a given set of source, target, and diffused heats for each, find a threshold value
    that yeilds a "linker" set of the given size (relative to the input set size).

    Returns:
        The cutoff/threshold to use, and the Relevance Score at that cutoff

    >>> find_linker_cutoff( set(["A", "B"]), set(["X", "Y"]), {"A":1.0, "B":1.1, "C":0.5, "D":0.4}, {"X":2.0, "Y":2.1, "C":0.7, "D":0.5}, 0.2)
    (0.4999, 0.16666666666666666)
    >>> find_linker_cutoff( set(["A", "B"]), set(["X", "Y"]), {"A":1.0, "B":1.1, "C":0.5, "D":0.4}, {"X":2.0, "Y":2.1, "C":0.7, "D":0.5}, 1.0)
    (0, 0)
    >>> find_linker_cutoff( set(["A", "B"]), set(["X", "Y"]), {"A":1.0, "B":1.1, "C":0.5, "D":0.4}, {"X":2.0, "Y":2.1, "C":0.7, "D":0.5}, 0.0)
    (1000000, 0)

    """
    if down_heat_diffused is None:
        # diffusing from a single source (i.e. not TieDIE but the HotNet algorithm, for comparison)
        cutoff, score = find_linker_cutoff_single(
            source_set, up_heat_diffused, size
        )
    else:
        try:
            cutoff, score = find_linker_cutoff_multi(
                source_set,
                target_set,
                up_heat_diffused,
                down_heat_diffused,
                size,
            )
        except Exception as e:  # noqa: BLE001 - intentional fallback with warning
            global _linker_cutoff_error_count, _linker_cutoff_error_shown
            _linker_cutoff_error_count += 1
            if not _linker_cutoff_error_shown:
                warnings.warn(
                    f'find_linker_cutoff failed: {e}',
                    stacklevel=2,
                )
                _linker_cutoff_error_shown = True
            return (0, 0)

    return (cutoff, score)


def find_linker_cutoff_single(source_set, up_heat_diffused, size):
    """If diffusing from a single source (i.e. not TieDIE but the HotNet algorithm, the implementation is trivial"""

    source_set = set(source_set)
    up_sorted = sorted(up_heat_diffused, key=up_heat_diffused.get, reverse=True)

    EPSILON = 0.0001

    # we want to find this many diffused genes in total
    target_size = size * len(source_set)

    i = 1
    cutoff = None
    for gene, heat in sorted(
        up_heat_diffused.items(), key=operator.itemgetter(1), reverse=True
    ):
        # add an epsilon to the cutoff so that the threshold falls just above this gene
        cutoff = heat + EPSILON
        # above this cutoff, the remaining set of linker genes is stored here
        diffused_set = set(up_sorted[0:i])
        # if the unique set of linker genes is of the desired size, stop at this cutoff...
        if len(diffused_set.difference(source_set)) > target_size:
            break
        i += 1

    return (cutoff, 0)


def find_linker_cutoff_multi(
    source_set, target_set, up_heat_diffused, down_heat_diffused, size
):
    """Find linker cutoff using bidirectional heat diffusion."""
    if size == 0:
        return (1000000, 0)

    target_set = set(target_set)
    source_set = set(source_set)
    # reverse-sort both sets by the diffused heat values
    up_sorted = sorted(up_heat_diffused, key=up_heat_diffused.get, reverse=True)
    down_sorted = sorted(
        down_heat_diffused, key=down_heat_diffused.get, reverse=True
    )

    # rank the min heats, and decrement the cutoff, adding an additional gene at each step
    EPSILON = 0.0001

    f, min_heats = filter_linkers(up_heat_diffused, down_heat_diffused, 1)
    # Iterate through the reverse-sorted list of heats. Stop when the exclusive set of nodes is below the desired size
    for cutoff in [
        h - EPSILON
        for (l, h) in sorted(
            min_heats.items(), key=operator.itemgetter(1), reverse=True
        )
    ]:
        score, size_frac = score_linkers(
            up_heat_diffused,
            up_sorted,
            down_heat_diffused,
            down_sorted,
            source_set,
            target_set,
            cutoff,
            size,
        )

        # reached the desired size: return the score & cutoff
        if size_frac > 1:
            return (cutoff, score)


def score_linkers(
    heats1, sorted1, heats2, sorted2, sourceSet, targetSet, cutoff, size
):
    """Get linkers greater than this cutoff according to reverse-sorted list.

    Inputs:
        source and target sets, diffused heats for each and the heat-sorted
        order for each.
        The linker cutoff chosen.
    """

    # find the genes in the first set that fall above this cutoff
    filtered_h1 = {}
    for l in sorted1:
        s = heats1[l]
        if s < cutoff:
            break

        filtered_h1[l] = s

    # genes in second set above this cutoff
    filtered_h2 = {}
    for l in sorted2:
        s = heats2[l]
        if s < cutoff:
            break

        filtered_h2[l] = s

    # make sets of both 'relevance neighborhoods' R_s and R_t
    f1 = set(filtered_h1)
    f2 = set(filtered_h2)

    # the union are all the genes in the relevance neighborhoods of each set R_s U R_t
    union = f1.union(f2)
    intersection = f1.intersection(f2)
    # connecting genes are linkers not in the source or target.
    # intutively, this is the heat that flows to the same
    connecting = intersection.difference(sourceSet).difference(targetSet)
    # the score is the number of connecting 'linker' genes over the size of the entire
    # relevance neighborhoods
    score = len(connecting) / float(len(union))
    # the relative size of the connecting genes, compared to the input set sizes
    size_frac = (
        len(connecting) / float(len(sourceSet.union(targetSet)))
    ) / float(size)

    return (score, size_frac)


def score_linkers_multi(input_heats, min_heats, cutoff, size):
    """Get linkers greater than this cutoff according to reverse-sorted list.
    This version takes an arbitrary number of inputs.

    Inputs:
        input_heats: a dictionary of an arbitrary number of input heat sets.
        min_heats: pre-processed 'linker' heat values according to any particular
        linker function.
    """

    # get the set of all input genes
    all_inputs = set()
    for name in input_heats:
        all_inputs = all_inputs.union(input_heats[name].keys())

    # generate the set of linker genes according to the supplied heat cutoff.
    all_linkers = set()
    for gene, heat in sorted(
        min_heats.items(), key=operator.itemgetter(1), reverse=True
    ):
        if heat < cutoff:
            break
        all_linkers.add(gene)

    # generate the exclusive 'connecting' set of linker/non-input genes
    # and score based on the fractional criterion
    connecting = all_linkers.difference(all_inputs)
    score = len(connecting) / float(len(all_linkers))
    # the relative size of the connecting genes, compared to the input set sizes
    size_frac = (len(connecting) / float(len(all_inputs))) / float(size)

    return (score, size_frac)


def get_min_heats(diffused):
    """Gets the minimum heats for all genes, from a number of diffused heat vectors.

    Input:
        diffused = { 'set':{'gene1':heat1, 'gene2':...}

    Returns:
        A minimum-heat vector over all genes

    """

    mins = {}
    for file in diffused:
        # a hash of hashes: file is the index
        for gene, heat in diffused[file].items():
            if gene in mins:
                if mins[gene] > heat:
                    mins[gene] = heat
            else:
                mins[gene] = heat

    return mins


def get_max_heats(diffused):
    """Gets the maximum heats for all genes, from a number of diffused heat vectors.
    Input:
        diffused = { 'set':{'gene1':heat1, 'gene2':...}

    Returns:
        A max-heat vector over all genes

    """
    max = {}
    for file in diffused:
        for gene, heat in diffused[file].items():
            if gene in max:
                if max[gene] < heat:
                    max[gene] = heat
            else:
                max[gene] = heat

    return max


def filter_linkers(up_heats_diffused, down_heats_diffused, cutoff):
    """Use the min(diffused1, diffused2) function to return a list of genes
    that fall above that cutoff.
    Input:
        diffused heats for each set, and the numeric cutoff value

    Returns:
        a list of genes above the cutoff, a hash of minimum heat values
    """
    linkers = {}
    filtered = []
    if down_heats_diffused is None:
        # trivially: if this is a single list of diffused values, just return it
        return (up_heats_diffused.keys(), up_heats_diffused)

    for node in up_heats_diffused:
        if node not in down_heats_diffused:
            # it doesn't make the cut if it's not in both sets
            continue
        min_heat = min(up_heats_diffused[node], down_heats_diffused[node])
        linkers[node] = min_heat
        if min_heat > cutoff:
            filtered.append(node)

    return (filtered, linkers)


def map_ugraph_to_network(edge_list, network):
    """Map undirected edges to the network to form a subnetwork
    in the hash-key directed network format

    Input:
        edge_list: edges in (s,t) format
        network: network in {source:set( (int, target), ... )

    Returns:
        Subnetwork in the data structure format of network input
    """

    subnetwork = {}

    for s, t in edge_list:
        # find this equivalent edge(s) in the directed network
        # edges:
        if s in network:
            for i, nt in network[s]:
                if nt == t:
                    if s not in subnetwork:
                        subnetwork[s] = set()
                    subnetwork[s].add((i, t))

        if t in network:
            for i, nt in network[t]:
                if nt == s:
                    if t not in subnetwork:
                        subnetwork[t] = set()
                    subnetwork[t].add((i, s))

    return subnetwork


def connected_subnets(network, subnet_nodes):
    """Input:
        A network in hash[source] = set( (interaction, target), ... ) Form
        A set of nodes to use for edge selection

    Returns:
        An edgelist set (source, target)
        where both nodes are in the subset of interest

    >>> network = {}
    >>> network['S1'] = set()
    >>> network['S2'] = set()
    >>> network['T2'] = set()
    >>> network['T1'] = set()
    >>> network['T3'] = set()
    >>> network['S1'].add(('a>','T1'))
    >>> network['S2'].add(('a>','T2'))
    >>> network['T1'].add(('t|','T2'))
    >>> network['T2'].add(('a>','T1'))
    >>> network['T3'].add(('t>','G5'))
    >>> connected_subnets(network, set(['S1','T1','T2','T3','G5']))
    set([('S1', 'T1'), ('T1', 'T2'), ('T2', 'T1')])
    """
    edgelist = set()
    ugraph = set()

    for s in network:
        for i, t in network[s]:
            # ignore self-links
            if s == t:
                continue
            if s in subnet_nodes and t in subnet_nodes:
                edgelist.add((s, t))
                if (t, s) not in edgelist:
                    ugraph.add((s, t))

    # use networkx to find the largest connected sub graph
    G = nx.Graph()
    G.add_edges_from(list(ugraph))
    # get the biggest connected component, add edges between all
    validated_edges = set()
    for component in nx.connected_components(G):
        validated_nodes = component
        for s, t in edgelist:
            # validate both nodes
            if s in validated_nodes and t in validated_nodes:
                validated_edges.add((s, t))

        break

    return validated_edges


def connected_nodes(network, hot_nodes):
    """Call connected_subnets to restrict to connected nodes, and return just the nodes
    filtered in this step
    """

    nodes = set()
    for s, t in connected_subnets(network, hot_nodes):
        nodes.add(s)
        nodes.add(t)
    return nodes


def run_pcst(up_heats, down_heats, linker_genes, network_file):
    """Convert input to format used for PCST program.
    Requires BioNet R package to be installed
    """

    # convert up/down heats to p-values
    # find the maximum heat for any value
    # the BioNet package requires p-values for an input, so we have to 'fake' these
    # here, converting them from heats.
    s_up = sorted([v for k, v in up_heats.items()], reverse=True)
    s_down = sorted([v for k, v in down_heats.items()], reverse=True)

    if len(up_heats) > 0:
        max_heat = s_up[0]
        min_heat = s_up[-1]

        if len(s_down) > 0:
            if s_down[0] > max_heat:
                max_heat = s_down[0]
                min_heat = s_up[-1]
            if s_down[-1] > min_heat:
                min_heat = s_down[-1]
    else:
        max_heat = 1
        min_heat = 1

    # take the sqrt of the fold difference over the min
    normalized_max = math.sqrt(max_heat / min_heat)
    scores = {}
    # the order is important here: there may be overlap between the source, target
    # and linker sets. The linkers are the highest priority, over the source/target.
    for node in down_heats:
        heat = down_heats[node]
        normalized_heat = math.sqrt(heat / min_heat)
        pval = math.exp(
            normalized_heat * math.log(float('1e-10')) / normalized_max
        )
        scores[node] = str(pval)
    for node in up_heats:
        heat = up_heats[node]
        normalized_heat = math.sqrt(heat / min_heat)
        pval = math.exp(
            normalized_heat * math.log(float('1e-10')) / normalized_max
        )
        scores[node] = str(pval)
    for node in linker_genes:
        scores[node] = '1e-10'

    pid = str(os.geteuid())

    tmp_act = open('/tmp/tmp_act_' + pid + '.tab', 'w')
    for node in scores:
        tmp_act.write(node + '\t' + scores[node] + '\n')
    tmp_act.close()

    # PCST is implemented in the BioNet package, and requires R to run. Python will call this script and collect the output
    os.system(
        sys.path[0]
        + '/span.R --activities /tmp/tmp_act_'
        + pid
        + '.tab --network '
        + network_file
        + ' > /tmp/pcst_'
        + pid
        + '.tab 2>/dev/null'
    )

    pcst_network = []
    pcst_line = re.compile(r'\[\d+\]\s+(\S+)\s+\-\-\s+(\S+)\s+')
    pcst = open('/tmp/pcst_' + pid + '.tab')
    for line in pcst:
        m = pcst_line.match(line)
        if m:
            pcst_network.append((m.group(1), m.group(2)))
    pcst.close()

    return pcst_network


def write_network(net, out_file):
    """Write network to SIF file format."""
    out = open(out_file, 'w')
    for source in net:
        for int, target in net[source]:
            out.write('\t'.join([source, int, target]) + '\n')

    out.close()


def random_subnet(network, num_sources):
    """Take a random sample of nodes, of the specified size
    from the supplied network
    """
    sub = {}
    for source in random.sample(network, num_sources):
        sub[source] = network[source]

    return sub


def write_el(el, so, down_set, out_file):
    """Write edge list to file."""
    out = open(out_file, 'w')
    for source, int, target in el:
        out.write('\t'.join([source, int, target]) + '\n')
    out.close()

    if so is None or down_set is None:
        return

    out = open(out_file + '.txt', 'w')
    set2 = set()
    for source, int, target in el:
        if target in down_set:
            set2.add(target)
    out.write(so + '\t' + '\t'.join(set2) + '\n')
    out.close()


def write_na_file(file_name, hash_values, attr_name):
    """Write out a node-attribute file. Include the header
    attr_name, and use the supplied hash values.
    """
    fh = None
    try:
        fh = open(file_name, 'w')
    except OSError as e:
        raise Exception(
            "Error: couldn't open output NA file for writing:" + file_name
        ) from e

    fh.write(attr_name + '\n')
    for key in hash_values:
        # check data type: hash values should be numbers for .NA file
        try:
            float(hash_values[key])
        except ValueError as e:
            raise Exception('Error: bad input value') from e

        fh.write(key + ' = ' + str(hash_values[key]) + '\n')

    fh.close()


def sample_heats(heats):
    """Randomly sample 80% of heats for bootstrapping."""
    ss = int(len(heats) * 0.8)
    keys = random.sample(heats, ss)
    subset = {}
    for k in keys:
        subset[k] = heats[k]

    return subset


def get_network_nodes(network):
    """Take a network in hash-key format and return a set containing the
    nodes in it.
    """
    nodes = set()
    for s in network:
        nodes.add(s)
        for i, t in network[s]:
            nodes.add(t)
    return nodes


def normalize_heats(data):
    """Normalize absolute value sum of data hash to 1000"""
    FACTOR = 1000
    normalized = {}
    signs = {}
    sum = 0.0
    for event, val in data.items():
        sum += abs(val)

    for event, val in data.items():
        sign = '+'
        if val < 0:
            sign = '-'
        normalized[event] = FACTOR * abs(val) / sum
        signs[event] = sign

    return normalized
