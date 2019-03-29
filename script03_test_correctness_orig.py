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

def eventNodeType(eventNode):
    '''
    Returns 'S', 'T', 'D', L', or 'C'
    'S': speciation
    'D': duplication
    'T': transfer
    'L': loss
    'C': end event
    '''
    return eventNode[0]

def mappingNodeToStr(mappingNode):
    return "{}-{}".format(mappingNode[0], mappingNode[1])

def firstChild(eventNode):
    # ['T', ('p3', 'h2'), ('p4', 'h4'), 0.5]
    # returns a mapping node
    assert(eventNode[1][0] and eventNode[1][1])
    return eventNode[1]

def secondChild(eventNode):
    # ['T', ('p3', 'h2'), ('p4', 'h4'), 0.5]
    # ['C', (None, None), (None, None), 1.0]
    # returns a mapping node
    # do not use on null entities
    assert(eventNode[2][0] and eventNode[2][1])
    return eventNode[2]

def eventNodeToStr(eventNode):
    # ('S', ('p7', 'h1'), ('p8', 'h2'), 0.5], 2)
    # ('S', ('n32', 'm143'), ('n31', 'm144'))
    if eventNodeType(eventNode) == 'S':
        return "spec | {}-{}, {}-{}".format(
            eventNode[1][0], eventNode[1][1],
            eventNode[2][0], eventNode[2][1]
        )
    elif eventNodeType(eventNode) == 'T':
        return "tran | {}-{}, {}-{}".format(
            eventNode[1][0], eventNode[1][1],
            eventNode[2][0], eventNode[2][1]
        )
    elif eventNodeType(eventNode) == 'D':
        return "dupl | {}-{}, {}-{}".format(
            eventNode[1][0], eventNode[1][1],
            eventNode[2][0], eventNode[2][1]
        )
    elif eventNodeType(eventNode) == 'L':
        return "loss | {}-{}".format(
            eventNode[1][0], eventNode[1][1],
        )
    if eventNodeType(eventNode) == 'C':
        return "END"

def enumerate_from_mapping_node(recongraph, mapping_node) :
    # for each mapping node, yield each of its event children
    for event_node in recongraph[mapping_node]:
        if eventNodeType(event_node) in ['S', 'T', 'D']:
            for left_mapping_dict in enumerate_from_mapping_node(recongraph, firstChild(event_node)):
                for right_mapping_dict in enumerate_from_mapping_node(recongraph, secondChild(event_node)):
                    recon_tree = {}
                    recon_tree[mapping_node] = [event_node]
                    recon_tree.update(left_mapping_dict)
                    recon_tree.update(right_mapping_dict)
                    yield recon_tree
        elif eventNodeType(event_node) == 'L':
            for child_mapping_dict in enumerate_from_mapping_node(recongraph, firstChild(event_node)):
                recon_tree = {}
                recon_tree[mapping_node] = [event_node]
                recon_tree.update(child_mapping_dict)
                yield recon_tree
        elif eventNodeType(event_node) == 'C':
            recon_tree = {}
            recon_tree[mapping_node] = [event_node]
            yield recon_tree
        

def enumerate_recon_trees(recongraph, roots) :
    '''
    Given a reconciliation graph, enumerate every reconciliation trees in the graph.
    Note that reconciliation trees are just like reconciliation graph, but each mapping
    node can only have one event node child.
    '''
    for root in roots:
        for recon_tree in enumerate_from_mapping_node(recongraph, root):
            yield recon_tree, root

def recon_trees_diff(recon_tree_A, recon_tree_B):
    '''
    Return the symmetric set difference between
    the two trees.
    '''
    diff_count = 0
    for mapping_node_key in recon_tree_A:
        if mapping_node_key not in recon_tree_B :
            diff_count += 1
        elif recon_tree_A[mapping_node_key] != recon_tree_B[mapping_node_key] :
            diff_count += 1
    
    for mapping_node_key in recon_tree_B:
        if mapping_node_key not in recon_tree_A :
            diff_count += 1
        elif recon_tree_A[mapping_node_key] != recon_tree_B[mapping_node_key] :
            diff_count += 1
    
    return diff_count



def brute_force_find_diameter(recongraph, roots) :
    '''
    Given a reconciliation graph, find the diameter of the graph via enumreating all
    its reconciliation trees.
    '''
    max_diff_count = 0
    recon_trees = [recon_tree for recon_tree, root in enumerate_recon_trees(recongraph, roots)]
    for recon_tree_i in range(1, len(recon_trees)):
        for recon_tree_j in range(recon_tree_i):
            recon_tree_A = recon_trees[recon_tree_i]
            recon_tree_B = recon_trees[recon_tree_j]
            diff_count = recon_trees_diff(recon_tree_A, recon_tree_B)
            max_diff_count = max(diff_count, max_diff_count)
    
    return max_diff_count

if __name__ == '__main__' :
    import pprint
    pp = pprint.PrettyPrinter()

    # compare size 4
    print("=== Mismatch of tree size 4 ===")
    for file_id in range(0, 95+1):
        edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size4/test-size4-no%d.newick" % file_id, 2, 4, 2)
        assert(mpr_count == sum(1 for _ in enumerate_recon_trees(dtl_recon_graph, best_roots)))

        brute_force_diameter = brute_force_find_diameter(dtl_recon_graph, best_roots)
        # find normal diameter
        # The gene tree needs to be in node format, not edge format, so we find that now.
        # (This also puts the gene_tree into postorder, as an ordered dict)
        gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

        species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
        diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                           False, False)

        if brute_force_diameter != diameter_alg_diameter :
            outname = './errorTrees/no4-id%d.png' % file_id
            ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
            print("ID = ", file_id)
            print("Brute Force: ", brute_force_diameter)
            print("Diameter Alg: ", diameter_alg_diameter)
    # compare size 5
    print("=== Mismatch of tree size 5 ===")
    for file_id in range(0, 1080):
        edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no%d.newick" % file_id, 2, 4, 2)
        assert(mpr_count == sum(1 for _ in enumerate_recon_trees(dtl_recon_graph, best_roots)))

        brute_force_diameter = brute_force_find_diameter(dtl_recon_graph, best_roots)
        # find normal diameter
        # The gene tree needs to be in node format, not edge format, so we find that now.
        # (This also puts the gene_tree into postorder, as an ordered dict)
        gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")

        species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(edge_species_tree, "hTop")
        diameter_alg_diameter = Diameter.diameter_algorithm(species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
                                           False, False)

        if brute_force_diameter != diameter_alg_diameter :
            outname = './errorTrees/no5-id%d.png' % file_id
            ReconciliationVisualization.visualizeAndSave(dtl_recon_graph, outname)
            print("ID = ", file_id)
            print("Brute Force: ", brute_force_diameter)
            print("Diameter Alg: ", diameter_alg_diameter)
