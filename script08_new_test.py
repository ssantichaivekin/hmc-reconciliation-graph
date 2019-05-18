'''
This script is deprecated. Please refer to script09_test_refactor instead.
'''

import DTLReconGraph
import ReconciliationVisualization
import Diameter
import HistogramAlgDebugRefactor
import HistogramAlgTools
from Histogram import Histogram

if __name__ == '__main__' :
    # compare size 4
    print("=== Mismatch of tree size 4 ===")
    for file_id in range(0, 95+1):
        edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size4/test-size4-no%d.newick" % file_id, 2, 4, 2)
        assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))

        brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
        # find normal diameter
        # The gene tree needs to be in node format, not edge format, so we find that now.
        # (This also puts the gene_tree into postorder, as an ordered dict)
        gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

        species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
        diameter_alg_hist = HistogramAlgDebugRefactor.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                           False, False)

        if brute_force_hist != diameter_alg_hist :
            outname = './errorTrees/no4-id%d-%d%d%d.png' % (file_id, 2, 4, 2)
            ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
            print("ID = ", file_id)
            print("Brute Force: ", brute_force_hist)
            print("Diameter Alg: ", diameter_alg_hist)
    # compare size 5
    print("=== Mismatch of tree size 5 ===")
    for file_id in range(0, 1080):
        print "ID = ", file_id
        edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no%d.newick" % file_id, 2, 4, 2)
        assert(mpr_count == sum(1 for _ in HistogramAlgTools.BF_enumerate_MPRs(dtl_recon_graph, best_roots)))
        print "mpr", mpr_count

        brute_force_hist = HistogramAlgTools.BF_find_histogram(dtl_recon_graph, best_roots)
        # find normal diameter
        # The gene tree needs to be in node format, not edge format, so we find that now.
        # (This also puts the gene_tree into postorder, as an ordered dict)
        gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

        species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
        diameter_alg_hist = HistogramAlgDebugRefactor.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
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
