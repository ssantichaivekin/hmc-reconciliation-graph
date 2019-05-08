'''
Test the correct of the current implementation of maximum distance reconciliation.
We do this by assuming that the DP algorithm generates a valid reconciliation graph,
then we use the diamater algorithm to find the max distance (diameter) among the 
reconciliation trees in that reconciliation graph. We then compare the max distance
which is calculated by the algorithm with the max distance that is computed by 
enumerating all the reconciliation trees in the graph, finding the pairwise distance
of each.
'''

import DTLReconGraph
import ReconciliationVisualization
import Diameter
import HistogramAlgDebugRefactor
import HistogramAlgTools
from Histogram import Histogram

if __name__ == '__main__' :
    D, T, L = 2, 4, 2
    tree_counts = {4:96, 5:1080}
    for tree_size in [4, 5]:
        print "=== Mismatch of tree size %d ===" % tree_size
        for file_id in range(0, tree_counts[tree_size]):
            filepath = "./newickSample/size%d/test-size%d-no%d.newick" % (tree_size, tree_size, file_id)
            edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
                = DTLReconGraph.reconcile(filepath, 2, 4, 2)
            assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))
            print "ID = ", file_id
            print "mpr", mpr_count

            brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
            gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

            species_tree, species_tree_root, species_node_count \
                = Diameter.reformat_tree(edge_species_tree, "hTop")

            diameter_alg_hist = HistogramAlgDebugRefactor.diameter_algorithm(
                species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                False, False)

            if brute_force_hist != diameter_alg_hist :
                outname = './errorTrees/no5-id%d-%d%d%d.png' % (file_id, 2, 4, 2)
                ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
                expected_n_pairs = HistogramAlgTools.calculate_n_pairs(mpr_count)
                brute_force_n_pairs = HistogramAlgTools.count_mpr_pairs(brute_force_hist)
                diag_force_n_pairs = HistogramAlgTools.count_mpr_pairs(diameter_alg_hist)
                print "ID = ", file_id
                print "Expected pairs ", expected_n_pairs
                print "Brute Force: ", brute_force_hist, "pairs: ", brute_force_n_pairs
                print "Diameter Alg: ", diameter_alg_hist, "pairs: ", diag_force_n_pairs