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
import itertools
import random
random.seed(1)

# If you just want to look at and/or debug the output of only one file,
# then you should just set the tree_size, D, T, L, file_id, and just run
# the histogram algorithm once. (see the commented out code below)
if __name__ == '__main__' :
    # TODO: Make the tester reads all trees in the corresponding folder instead
    # of requiring the user to specify the tree.
    tree_counts = {4:96, 5: 1080, 6: 25920, 7: 583, 8:250}
    # test from small to large tree sizes
    for tree_size in sorted(tree_counts.keys()):
        # test different D, T, L values in {1, 2, 3, 4}
        for D, T, L in itertools.product([1, 2, 3, 4], repeat=3):
            print "=== Mismatch of tree size %d : DTL %d %d %d ===" % (tree_size, D, T, L)
            file_ids = range(0, tree_counts[tree_size])
            for file_id in file_ids:
                # find the corresponding newick test file generated via script02
                filepath = "./newickSample/size%d/test-size%d-no%d.newick" % (tree_size, tree_size, file_id)
                # from the newick tree create the reconciliation graph
                edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
                    = DTLReconGraph.reconcile(filepath, D, T, L)
                # sanity check: the mpr_count returned is equal to the count generated via brute force
                assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))

                # Calculate the histogram via brute force
                brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)

                # Reformat the host and parasite tree to use it with the histogram algorithm
                gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
                species_tree, species_tree_root, species_node_count \
                    = Diameter.reformat_tree(edge_species_tree, "hTop")

                # Calculate the histogram via histogram algorithm
                diameter_alg_hist = HistogramAlgDebugRefactor.diameter_algorithm(
                    species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                    False, False)

                # If there is a mismatch, print the details and save the tree that causes
                # the error to a folder called errorTrees.
                if brute_force_hist != diameter_alg_hist :
                    outname = './errorTrees/no%d-id%d-%d%d%d.png' % (tree_size, file_id, D, T, L)
                    ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
                    expected_n_pairs = HistogramAlgTools.calculate_n_pairs(mpr_count)
                    brute_force_n_pairs = HistogramAlgTools.count_mpr_pairs(brute_force_hist)
                    diag_force_n_pairs = HistogramAlgTools.count_mpr_pairs(diameter_alg_hist)
                    print "ID = ", file_id, "DTL = ", D, T, L
                    print "Expected pairs ", expected_n_pairs
                    print "Brute Force: ", brute_force_hist, "pairs: ", brute_force_n_pairs
                    print "Diameter Alg: ", diameter_alg_hist, "pairs: ", diag_force_n_pairs
                    print ""

# if __name__ == '__main__' :
#     tree_size = 6
#     D, T, L = 1, 1, 1
#     print "=== Mismatch of tree size %d : DTL %d %d %d ===" % (tree_size, D, T, L)
#     file_id = 23258
#     filepath = "./newickSample/size%d/test-size%d-no%d.newick" % (tree_size, tree_size, file_id)
#     edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
#         = DTLReconGraph.reconcile(filepath, 2, 4, 2)
#     assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))
#     # print "ID = ", file_id, "DTL = ", D, T, L
#     # print "mpr", mpr_count

#     brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
#     gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

#     species_tree, species_tree_root, species_node_count \
#         = Diameter.reformat_tree(edge_species_tree, "hTop")

#     diameter_alg_hist = HistogramAlgDebugRefactor.diameter_algorithm(
#         species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
#         False, False)

#     if brute_force_hist != diameter_alg_hist :
#         outname = './errorTrees/no%d-id%d-%d%d%d.png' % (tree_size, file_id, D, T, L)
#         ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
#         expected_n_pairs = HistogramAlgTools.calculate_n_pairs(mpr_count)
#         brute_force_n_pairs = HistogramAlgTools.count_mpr_pairs(brute_force_hist)
#         diag_force_n_pairs = HistogramAlgTools.count_mpr_pairs(diameter_alg_hist)
#         print "ID = ", file_id, "DTL = ", D, T, L
#         print "Expected pairs ", expected_n_pairs
#         print "Brute Force: ", brute_force_hist, "pairs: ", brute_force_n_pairs
#         print "Diameter Alg: ", diameter_alg_hist, "pairs: ", diag_force_n_pairs
#         print ""