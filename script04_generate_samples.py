'''
This script is deprecated. Please refer to script09_test_refactor instead.
'''

from script03_test_correctness_orig import *
import pprint

def generate_lots_of_trees_to_study():
    dup_cost = 1
    trans_cost = 4
    loss_cost = 1
    print("=== Mismatch of tree size 4 ===")
    for file_id in range(1, 90):
        edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile(
            "./newickSample/size4/test-size4-no%d.newick" % file_id, dup_cost, trans_cost, loss_cost)
        assert(mpr_count == sum(1 for _ in enumerate_recon_trees(dtl_recon_graph, best_roots)))

        brute_force_diameter = brute_force_find_diameter(dtl_recon_graph, best_roots)
        # find normal diameter
        # The gene tree needs to be in node format, not edge format, so we find that now.
        # (This also puts the gene_tree into postorder, as an ordered dict)
        gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

        species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
        diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                           False, False)

        if brute_force_diameter != 0 :
            outname = './errorTrees/no4-id%d-%d-%d-%d.png' % (file_id, dup_cost, trans_cost, loss_cost)
            ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
            print("ID = ", file_id)
            print("Brute Force: ", brute_force_diameter)
            print("Diameter Alg: ", diameter_alg_diameter)

def generate_recon_graph_info(tree_size, file_id, dup_cost, trans_cost, loss_cost):
    pp = pprint.PrettyPrinter()
    print("Tree Size = %d" % tree_size)
    print("File ID = %d" % file_id)
    print("DUP %d TRANS %d LOSS %s" % (dup_cost, trans_cost, loss_cost))

    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile(
            "./newickSample/size%d/test-size%d-no%d.newick" % (tree_size, tree_size, file_id), dup_cost, trans_cost, loss_cost)
    
    for recon_tree, root in enumerate_recon_trees(dtl_recon_graph, best_roots):
        print("Root:")
        pp.pprint(root)
        print("Tree")
        pp.pprint(recon_tree)

def generate_test_csv(tree_size):
    import csv
    with open('scratch_paper.csv', 'w') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(['p', 's1', 's2'])
        nodes = tree_size * 2 - 1
        for i in range(nodes):
            for j in range(nodes):
                for k in range(nodes):
                    csvwriter.writerow([i, j, k])

    

if __name__ == '__main__':
    # I have choosen size 4 id 55 with dup = 1, trans = 4, loss = 1 to test.
    # I think it tells good story about the algorithm and is very simple.
    # I will attempt it by hand and compare it with various algorithms.
    # generate_recon_graph_info(4, 55, 1, 4, 1)
    # Another example tree would be size 4 id 2
    # generate_recon_graph_info(4, 2, 3, 3, 1)

    generate_test_csv(4)
        