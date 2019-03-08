Write test cases by generating all possible trees of fixed number of leaves and try different values of duplication, transfer, loss

1. We do this by generating trees using recusrive algorithm:
    A tree of size 6 can be constructed in one of the following ways
        1. The left subtree is of size 1, and the right is of size 5
        2. The left subtree is of size 2, and the right is of size 4
        3. The left subtree is of size 3, and the right is of size 3
    Note that we do not use (4,2) because we believe it would create 
    a similar tree to (2,4). We then do this recursively, creating different 
    trees, and then save them for future uses.

2. Then, we create two trees, and we tag them as appropriate. We then link the leaves in randomized manner, that is
we link it using all permutations. Note that this will create many same trees in the process.

3. We parse the trees to newick format

Now, we run it through the reconciliation graph creation. We generally want graphs with low diameter but with
many possible trees. That is, we want it to be complex and intertwined. We will probably sort all of them,
and then computationally (probably choosing via some matrix of diameter and possible trees -- possible_trees / diameter?)
pick 3 of size 4, and 1 of size 5, and 1 of size 6 per each dup-preferred, transfer-preferred, or loss-preferred category.
Maybe at the end, we will opt for mixed category instead.

I thikn we should end up with around 10 small tests and 3 bigger tests in the end.

My test reconciliations shall be saved in json format so that it is easy to load and rerun.