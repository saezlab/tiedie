"""TieDIE: Tied Diffusion for Network Discovery

Command-line interface for the TieDIE algorithm.

Authors:
    Evan Paull (epaull@soe.ucsc.edu)

Minimum Inputs:
    - separate source/target input heat files: tab-separated, 3 columns each
      with <gene> <input heat> <sign (+/-)>
    - a search pathway in .sif format (geneA <interaction> geneB)

Outputs:
    Creates a directory in the current working directory, and writes all
    output to that. Information and warnings are logged to standard error.
"""

import os
import sys
import argparse

from .ppr import PPrDiffuser
from .util import (
    run_pcst,
    write_el,
    parse_net,
    search_dfs,
    parse_heats,
    write_na_file,
    write_network,
    classify_state,
    filter_linkers,
    get_out_degrees,
    normalize_heats,
    connected_subnets,
    get_network_nodes,
    find_linker_cutoff,
    map_ugraph_to_network,
)
from .kernel import Kernel
from .permute import NetBalancedPermuter
from .master_reg import ActivityScores
from .kernel_scipy import SciPYKernel


def extract_subnetwork(
    up_heats,
    down_heats,
    up_heats_diffused,
    down_heats_diffused,
    size_control,
    set_alpha,
    network,
    use_pcst=False,
    network_file=None,
):
    """Generate a spanning subnetwork from the supplied inputs, diffused heats
    and size control cutoff.

    Args:
        up_heats: upstream heats
        down_heats: downstream heats
        up_heats_diffused: diffused upstream heats
        down_heats_diffused: diffused downstream heats
        size_control: size control factor
        set_alpha: optional linker cutoff override
        network: parsed network
        use_pcst: use Prize-Collecting Steiner Tree formulation
        network_file: path to network file (required if use_pcst=True)

    Returns:
        Tuple of (spanning network, nodes in network, alpha score, linker scores)
    """

    linker_cutoff = None
    alpha_score = None

    if set_alpha:
        alpha_score = None
        linker_cutoff = float(set_alpha)
    else:
        linker_cutoff, alpha_score = find_linker_cutoff(
            up_heats,
            down_heats,
            up_heats_diffused,
            down_heats_diffused,
            size_control,
        )

    linker_nodes, linker_scores = filter_linkers(
        up_heats_diffused,
        down_heats_diffused,
        linker_cutoff,
    )

    ugraph = None

    if use_pcst:
        ugraph = run_pcst(up_heats, down_heats, linker_nodes, network_file)
    else:
        nodes = set(up_heats).union(set(down_heats)).union(set(linker_nodes))
        ugraph = connected_subnets(network, nodes)

    if len(ugraph) == 0:
        sys.stderr.write(
            "Couldn't find any linking graph at this size setting!\n"
        )
        return (None, None, None, None)

    subnet_soln = map_ugraph_to_network(ugraph, network)

    subnet_soln_nodes = set()

    for s in subnet_soln:
        subnet_soln_nodes.add(s)

        for i, t in subnet_soln[s]:
            subnet_soln_nodes.add(t)

    return (subnet_soln, subnet_soln_nodes, alpha_score, linker_scores)


def find_consistent_paths(
    up_signs,
    down_signs,
    search_network,
    output_folder,
    search_depth,
    output=True,
):
    """Filter the heat-generated network by searching for all directed paths
    from each source to each target gene.

    Args:
        up_signs: hash with up/down signs for each upstream node
        down_signs: hash with up/down signs for each downstream node
        search_network: heat derived subnetwork for DFS
        output_folder: write output networks under this directory
        search_depth: maximum search depth
        output: flag to write output to the specified folder

    Returns:
        Tuple of (TP count, FP count, validated edge list)
    """

    gene_states, t_states = classify_state(up_signs, down_signs)
    validated = set()
    down_set = set(down_signs.keys())
    TP = 0
    FP = 0

    for source in up_signs:
        action = gene_states[source]
        falsePaths = []
        truePaths = []
        edges_this_source = set()

        search_dfs(
            source,
            action,
            edges_this_source,
            set(),
            down_set,
            search_network,
            gene_states,
            t_states,
            search_depth,
            truePaths,
            falsePaths,
            False,
        )

        for edge in edges_this_source:
            validated.add(edge)

        TP += len(truePaths)
        FP += len(falsePaths)

        if output:
            out_file = output_folder + '/' + source + '.cn.sif'
            sys.stderr.write(
                'Writing Single Causal Neighborhood to ' + out_file + '\n'
            )
            write_el(edges_this_source, source, down_set, out_file)

    if output:
        out_file = output_folder + '/tiedie.cn.sif'
        sys.stderr.write(
            'Writing Full Causal Neighborhood to ' + out_file + '\n'
        )
        write_el(validated, 'ALL', down_set, out_file)

    return (TP, FP, validated)


def score_subnet(subnet_soln_nodes, up_heats, down_heats, report_fh):
    """Score sets according to a Compactness Score that weighs the coverage of
    source and target sets while penalizing for the number of linker nodes
    needed to connect them.
    """

    S = set(up_heats.keys())
    T = set(down_heats.keys())
    C = subnet_soln_nodes.difference(S).difference(T)
    U = S.union(T)
    Sr = S.intersection(subnet_soln_nodes)
    Tr = T.intersection(subnet_soln_nodes)
    PENALTY_CONST = 0.1
    penalty = (float(len(C)) / len(U)) * PENALTY_CONST
    score = (
        float(len(Sr)) / (len(S) * 2) + float(len(Tr)) / (len(T) * 2) - penalty
    )

    report_fh.write(
        str(float(len(Sr)) / len(S))
        + '\t'
        + 'of source nodes'
        + str(len(Sr))
        + ' out of '
        + str(len(S))
        + '\n'
    )
    report_fh.write(
        str(float(len(Tr)) / len(T))
        + '\t'
        + 'of target nodes'
        + str(len(Tr))
        + ' out of '
        + str(len(T))
        + '\n'
    )
    report_fh.write('And ' + str(len(C)) + ' connecting nodes\n')

    return score


def create_parser():
    """Create argument parser for TieDIE CLI."""

    parser = argparse.ArgumentParser(
        description='TieDIE: Tied Diffusion for Network Discovery',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        '-n',
        '--network',
        required=True,
        help='.sif network file for the curated pathway to search',
    )
    parser.add_argument(
        '-u',
        '--up_heats',
        required=True,
        help='File with upstream heats: <gene> <heat> <sign (+/-)>',
    )
    parser.add_argument(
        '-d',
        '--down_heats',
        help='File with downstream heats: <gene> <heat> <sign (+/-)>',
    )
    parser.add_argument(
        '--d_expr',
        help='Differential expression file (alternative to --down_heats)',
    )
    parser.add_argument(
        '-k',
        '--kernel',
        help='Pre-computed heat diffusion kernel file',
    )
    parser.add_argument(
        '-s',
        '--size',
        type=float,
        default=1.0,
        help='Network size control factor (default: 1.0)',
    )
    parser.add_argument(
        '-a',
        '--alpha',
        help='Linker cutoff (overrides size factor)',
    )
    parser.add_argument(
        '-c',
        '--depth',
        type=int,
        default=3,
        help='Search depth for causal paths (default: 3)',
    )
    parser.add_argument(
        '-p',
        '--permute',
        type=int,
        default=1000,
        help='Number of random permutations (default: 1000)',
    )
    parser.add_argument(
        '-m',
        '--min_hub',
        type=int,
        help='Minimum genes in regulon for TF (required with --d_expr)',
    )
    parser.add_argument(
        '--pagerank',
        action='store_true',
        help='Use Personalized PageRank to diffuse',
    )
    parser.add_argument(
        '--pcst',
        action='store_true',
        help='Use Prize-Collecting Steiner Tree formulation',
    )
    parser.add_argument(
        '--output_folder',
        default='TieDIE',
        help='Output folder (default: TieDIE)',
    )

    return parser


def main(args=None):
    """Main entry point for TieDIE CLI."""

    parser = create_parser()
    opts = parser.parse_args(args)

    # Input validation
    if opts.d_expr is None and opts.down_heats is None:
        parser.error('Must supply either --down_heats or --d_expr')

    if opts.down_heats is not None and opts.d_expr is not None:
        parser.error('Cannot supply both --down_heats and --d_expr')

    if opts.d_expr and not opts.min_hub:
        parser.error('--min_hub required when using --d_expr')

    if opts.kernel is None:
        sys.stderr.write(
            'Warning: No kernel file supplied, will use SCIPY to compute '
            'the matrix exponential, t=0.1...\n'
        )

    # Parse network
    sys.stderr.write('Parsing Network File..\n')
    network = parse_net(opts.network)
    network_nodes = get_network_nodes(network)

    # Parse upstream heats
    up_heats, up_signs = parse_heats(opts.up_heats, network_nodes)
    up_heats = normalize_heats(up_heats)

    # Parse or compute downstream heats
    down_heats = None
    down_signs = None

    if opts.d_expr:
        down_heats = ActivityScores.findRegulators(
            network,
            opts.d_expr,
            min_hub=opts.min_hub,
        )
        down_signs = {}

        for g, h in down_heats.items():
            if h < 0:
                down_signs[g] = '-'
            else:
                down_signs[g] = '+'

            down_heats[g] = abs(h)

    else:
        down_heats, down_signs = parse_heats(opts.down_heats, network_nodes)

    down_heats = normalize_heats(down_heats)

    # Create output folder
    output_folder = opts.output_folder

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Create diffuser
    if opts.pagerank:
        diffuser = PPrDiffuser(network)
    elif opts.kernel is not None:
        sys.stderr.write('Loading Heat Diffusion Kernel..\n')
        diffuser = Kernel(opts.kernel)
    else:
        sys.stderr.write(
            'Using SCIPY to compute the matrix exponential, t=0.1...\n'
        )
        diffuser = SciPYKernel(opts.network)

    # Validate kernel labels match network
    k_labels = diffuser.getLabels()

    if len(network_nodes) != len(k_labels) or len(
        network_nodes.intersection(k_labels)
    ) != len(k_labels):
        sys.stderr.write(
            'Error: the universe of gene/node labels in the network file '
            "doesn't match the supplied kernel file!\n"
        )
        sys.exit(1)

    # Diffuse heats
    sys.stderr.write('Diffusing Heats...\n')
    up_heats_diffused = diffuser.diffuse(up_heats, reverse=False)
    down_heats_diffused = diffuser.diffuse(down_heats, reverse=True)

    # Extract subnetwork
    subnet_soln, subnet_soln_nodes, alpha_score, linker_scores = (
        extract_subnetwork(
            up_heats,
            down_heats,
            up_heats_diffused,
            down_heats_diffused,
            opts.size,
            opts.alpha,
            network,
            use_pcst=opts.pcst,
            network_file=opts.network,
        )
    )

    # Generate linker stats
    out_degrees = get_out_degrees(subnet_soln)
    sys.stderr.write(
        'Writing network node stats to ' + output_folder + '/node.stats\n'
    )
    out_file = output_folder + '/node.stats'

    with open(out_file, 'w') as out:
        out.write('NODE\tCONNECTING\tMIN_HEAT\tOUT_DEGREE\n')
        node_types = {}

        for node in subnet_soln_nodes:
            out_deg = out_degrees[node]
            linker_heat = linker_scores[node]
            connecting = '0'

            if node in up_heats:
                node_types[node] = 1
            elif node in down_heats:
                node_types[node] = -1

            if node not in up_heats:
                if down_heats is not None and node not in down_heats:
                    connecting = '1'
                    node_types[node] = 0

            out.write(
                '\t'.join([node, connecting, str(linker_heat), str(out_deg)])
                + '\n'
            )

    # Write Cytoscape files
    write_na_file(output_folder + '/node_types.NA', node_types, 'NodeTypes')
    write_na_file(output_folder + '/heats.NA', linker_scores, 'LinkerHeats')

    # Write network
    sys.stderr.write('Writing ' + output_folder + '/tiedie.sif result\n')
    write_network(subnet_soln, output_folder + '/tiedie.sif')

    # Find consistent paths
    TP, FP, validated = find_consistent_paths(
        up_signs,
        down_signs,
        subnet_soln,
        output_folder,
        opts.depth,
        output=True,
    )

    # Write report
    report_file = output_folder + '/report.txt'

    with open(report_file, 'w') as report_fh:
        sys.stderr.write(
            'Writing Report to ' + report_file + ' :compactness analysis\n'
        )
        score = score_subnet(subnet_soln_nodes, up_heats, down_heats, report_fh)
        report_fh.write('Compactness Score:' + str(score) + '\n')

        # Permutation test
        sys.stderr.write(
            'Running permutation tests... (could take several minutes for '
            'inputs of hundreds of genes @1000 permutations)\n'
        )
        perObj = NetBalancedPermuter(network, up_heats)
        permutedHeats = perObj.permute(opts.permute)
        permuted_scores = []

        for heats in permutedHeats:
            diffused_heats = diffuser.diffuse(heats)
            cutoff, perm_score = find_linker_cutoff(
                heats,
                down_heats,
                diffused_heats,
                down_heats_diffused,
                opts.size,
            )
            permuted_scores.append(perm_score)

        # Write score files
        with open(output_folder + '/score.txt', 'w') as sig_fh:
            sig_fh.write(str(alpha_score) + '\n')

        with open(output_folder + '/permuted_scores.txt', 'w') as sig_fh:
            for val in sorted(permuted_scores, reverse=True):
                sig_fh.write(str(val) + '\n')

        # Calculate p-value
        no_gte = 0.0

        for val in sorted(permuted_scores, reverse=True):
            if val >= alpha_score:
                no_gte += 1
            else:
                break

        pval = (no_gte + 1) / (opts.permute + 1)
        sys.stderr.write(
            'Writing Report to ' + report_file + ' :empirical p-value...\n'
        )
        report_fh.write(
            'P-value: '
            + str(pval)
            + ' (with '
            + str(opts.permute)
            + ' random permutations)\n'
        )


if __name__ == '__main__':
    main()
