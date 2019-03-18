'''
Test the correct of the current implementation of maximum distance reconciliation.
We do this by assuming that the DP algorithm generates a valid reconciliation graph,
then we use the diamater algorithm to find the max distance (diameter) among the 
reconciliation trees in that reconciliation graph. We then compare the max distance
which is calculated by the algorithm with the max distance that is computed by 
enumerating all the reconciliation trees in the graph, finding the pairwise distance
of each.
'''


def enumerate_recon_trees(recongraph) :
    yield {}

def brute_force_find_diameter(recongraph) :
    return

def compare_alg_diameter_with_brute_force(recongraph) :
    return

if __name__ == '__main__' :
    # compare size 3
    pass
    # compare size 4
    pass
    # compare size 5
    pass
