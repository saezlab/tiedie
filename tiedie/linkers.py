import sys
import operator

from .util import connected_subnets, map_ugraph_to_network


def min(vals):
    """Get minimum value from a list."""
    min = vals[0]
    for v in vals:
        if v < min:
            min = v

    return min


def get_product(diffused):
    """Compute product of heat values across diffused vectors."""
    gene_scores = {}
    for file in diffused:
        # a hash of hashes: file is the index
        for gene, heat in diffused[file].items():
            if gene not in gene_scores:
                gene_scores[gene] = []
            gene_scores[gene].append(heat)

    gene_products = {}
    for gene in gene_scores:
        product = 1
        for v in gene_scores[gene]:
            product *= v
        gene_products[gene] = product

    return gene_products


def get_min_heats(consider_top, diffused):
    """Gets the minimum heats for all genes, from a number of diffused heat vectors.

    Input:
        diffused = { 'set':{'gene1':heat1, 'gene2':...}

    Returns:
        A minimum-heat vector over all genes

    """

    gene_scores = {}
    for file in diffused:
        # a hash of hashes: file is the index
        for gene, heat in diffused[file].items():
            if gene not in gene_scores:
                gene_scores[gene] = []
            gene_scores[gene].append(heat)

    min_gene_values = {}
    for gene in gene_scores:
        values = gene_scores[gene]
        # get the top X
        min_gene_values[gene] = min(
            sorted(values, reverse=True)[0:consider_top]
        )

    return min_gene_values


def get_max_heats(consider_top, diffused):
    """Gets the maximum heats for all genes, from a number of diffused heat vectors.

    Input:
        diffused = { 'set':{'gene1':heat1, 'gene2':...}

    Returns:
        A minimum-heat vector over all genes

    """

    gene_scores = {}
    for file in diffused:
        # a hash of hashes: file is the index
        for gene, heat in diffused[file].items():
            if gene not in gene_scores:
                gene_scores[gene] = []
            gene_scores[gene].append(heat)

    max_gene_values = {}
    for gene in gene_scores:
        values = gene_scores[gene]
        # get the top X
        max_gene_values[gene] = max(
            sorted(values, reverse=True)[0:consider_top]
        )

    return max_gene_values


def extract_subnetwork(
    network, input_heats, diffused_heats, size_control, opts
):
    """Generate a spanning subnetwork from input heats and diffused heats.

    Args:
        network: Dict mapping source nodes to sets of (interaction, target) tuples.
        input_heats: Dict of input heat sets by name (e.g., 'source', 'target').
        diffused_heats: Dict of diffused heat values by input set name.
        size_control: Size control factor for linker cutoff.
        opts: Additional options dict.

    Returns:
        Tuple of (subnetwork, node_set, linker_heats, linker_cutoff).

    Example:
        ```python
        input_heats = {
            'source': {'s1': 0.01, 's2': 0.5},
            'target': {'e1': 0.49}
        }
        diffused_heats = {
            'source': {'s1': 0.05, 's2': 0.4},
            'target': {'e1': 0.4, 't2': 0.3, 't1': 0.1}
        }
        network = {
            's1': {('-t>', 't1')},
            's2': {('-t>', 't2')},
            't2': {('-t>', 'e1')}
        }
        subnet, nodes, heats, cutoff = extract_subnetwork(
            network, input_heats, diffused_heats, 0.25, {}
        )
        ```
    """

    linker_cutoff = None
    linker_scores = None

    # get linker heats as a function of input sets
    consider = len(input_heats)

    linker_heats = get_min_heats(consider, diffused_heats)
    # linker_heats = get_product(diffused_heats)

    EPSILON = 0.0001
    linker_cutoff = None
    score = None
    linkers = set()
    linker_scores = {}
    for l, h in sorted(
        linker_heats.items(), key=operator.itemgetter(1), reverse=True
    ):
        c = h - EPSILON
        score, size_frac = linker_score(
            input_heats, linker_heats, c, size_control
        )
        linker_cutoff = c
        linkers.add(l)
        linker_scores[l] = h
        if size_frac > 1:
            break

    input_genes = set()
    for input in input_heats:
        input_genes = input_genes.union(input_heats[input].keys())
    # set of input heats
    ugraph = None
    # USE LINKER GENES AND INPUT GENES
    active_nodes = set(linkers)
    active_nodes = active_nodes.union(input_genes)
    ugraph = connected_subnets(network, active_nodes)
    if len(ugraph) == 0:
        sys.stderr.write(
            "Couldn't find any linking graph at this size setting!\n"
        )
    subnet_soln = map_ugraph_to_network(ugraph, network)

    subnet_soln_nodes = set()
    for s in subnet_soln:
        subnet_soln_nodes.add(s)
        for i, t in subnet_soln[s]:
            subnet_soln_nodes.add(t)

    return (subnet_soln, subnet_soln_nodes, linker_heats, linker_cutoff)


def linker_score(input_heats, min_heats, cutoff, size):
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
