'''
This script is deprecated. Please refer to script09_test_refactor instead.
'''

import DTLReconGraph
import ReconciliationVisualization
import Diameter
import HistogramAlg
import HistogramAlgTools
from Histogram import Histogram

if __name__ == '__main__' :
    import pprint
    pp = pprint.PrettyPrinter()

    file_id = 936
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no%d.newick" % file_id, 2, 4, 2)
    assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))

    brute_force_histogram = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
    brute_force_diameter = HistogramAlgTools.BF_find_diameter(dtl_recon_graph, best_roots)

    # find normal diameter
    # The gene tree needs to be in node format, not edge format, so we find that now.
    # (This also puts the gene_tree into postorder, as an ordered dict)
    gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

    species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
    diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        False, False)
    diameter_alg_histogram = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        False, False)
    assert(brute_force_diameter == diameter_alg_diameter)
    assert(brute_force_histogram == diameter_alg_histogram)