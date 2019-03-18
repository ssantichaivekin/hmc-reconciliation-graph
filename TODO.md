Write test cases by generating all possible trees of fixed number of leaves and try different values of duplication, transfer, loss

- [x] We do this by generating trees using recusrive algorithm:
    A tree of size 6 can be constructed in one of the following ways
        1. The left subtree is of size 1, and the right is of size 5
        2. The left subtree is of size 2, and the right is of size 4
        3. The left subtree is of size 3, and the right is of size 3
    Note that we do not use (4,2) because we believe it would create 
    a similar tree to (2,4). We then do this recursively, creating different 
    trees, and then save them for future uses.

- [x] Then, we create two trees, and we tag them as appropriate. We then link the leaves in randomized manner, that is
we link it using all permutations. Note that this will create many same trees in the process.

- [x] We parse the trees to newick format

- [ ] Now, we run it through the reconciliation graph creation. We generally want graphs with low diameter but with
many possible trees. We then test the corrctness of the original algorithm.