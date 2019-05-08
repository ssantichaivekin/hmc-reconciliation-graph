import DTLReconGraph
import ReconciliationVisualization
import Diameter
import HistogramAlg
import HistogramAlgTools
from Histogram import Histogram

if __name__ == '__main__' :
    import pprint
    pp = pprint.PrettyPrinter()

    file_id = 40
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no%d.newick" % file_id, 2, 4, 2)

    gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
    species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
    
    diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        True, False)
    diameter_alg_hist = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        True, False)

    brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
    # if brute_force_hist != diameter_alg_hist :
    outname = './errorTrees/no5-id%d.png' % file_id
    ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
    expected_n_pairs = HistogramAlgTools.calculate_n_pairs(mpr_count)
    brute_force_n_pairs = HistogramAlgTools.count_mpr_pairs(brute_force_hist)
    diag_force_n_pairs = HistogramAlgTools.count_mpr_pairs(diameter_alg_hist)
    print "ID = ", file_id
    print "Expected pairs ", expected_n_pairs
    print "Brute Force: ", brute_force_hist, "pairs: ", brute_force_n_pairs
    print "Diameter Alg: ", diameter_alg_hist, "pairs: ", diag_force_n_pairs