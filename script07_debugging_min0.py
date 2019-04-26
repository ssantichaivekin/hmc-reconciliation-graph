from script05_test_correctness_new import *

if __name__ == '__main__' :
    import pprint
    pp = pprint.PrettyPrinter()

    file_id = 489
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no%d.newick" % file_id, 2, 4, 2)

    gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

    species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
    diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        True, False)
    diameter_alg_histogram = DiameterModified.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                        True, False)