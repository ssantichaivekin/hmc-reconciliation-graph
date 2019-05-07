'''
This script generates all possible (1) newick host tree 
(2) newick parasite tree and (3) mapping from host to 
parasite.

This script will save the outputs as a .newick files
in the specified folder.

Author: Santi Santichaivekin
'''

import itertools
import os

def generateTreeTuples(numLeaves, prefix, startNum):
    '''
    Generate trees in Newick format
    '''
    root = "{}{}".format(prefix, startNum)
    if numLeaves == 1 :
        yield root

    leftTreeStartNum = startNum + 1
    mid = numLeaves // 2
    for leftTreeLeaves in range(1, mid+1): # from 1 to mid inclusive
        rightTreeLeaves = numLeaves - leftTreeLeaves
        rightTreeStartNum = leftTreeStartNum + (2**leftTreeLeaves - 1) # add the left tree size
        for leftTreeTuple in generateTreeTuples(leftTreeLeaves, prefix, leftTreeStartNum) :
            for rightTreeTuple in generateTreeTuples(rightTreeLeaves, prefix, rightTreeStartNum):
                yield (root, leftTreeTuple, rightTreeTuple)

def treeTupleStrings(treeTuple):
    '''
    Convert a python representation of tree to a string in
    newick format.
    '''
    if type(treeTuple) is str :
        return treeTuple
    else :
        return "({},{}){}".format(
            treeTupleStrings(treeTuple[1]), # left
            treeTupleStrings(treeTuple[2]), # right
            treeTuple[0] # root
            )

def treeTupleLeaves(treeTuple):
    '''
    Return the list of leaves of the tree.
    '''
    if type(treeTuple) is str :
        leaves = [treeTuple]
        return leaves
    else :
        leaves = []
        # ignore the root treeTuple[0] !
        leaves += treeTupleLeaves(treeTuple[1])
        leaves += treeTupleLeaves(treeTuple[2])
        return leaves

def generateMappings(parasiteLeaves, hostLeaves):
    '''
    Return an iterator to all possible mappings
    from the parasite leaves to host leaves.
    '''
    for hostLeavesPermuted in itertools.permutations(hostLeaves) :
        yield list(zip(parasiteLeaves, hostLeavesPermuted))

def mappingToString(parasiteHostMapping):
    outStr = ""
    for parasiteLeaf, hostLeaf in parasiteHostMapping :
        outStr += "{}:{}\n".format(parasiteLeaf, hostLeaf)
    return outStr

def generateNewickTests(numLeaves, destFolderName):
    '''
    Generate all possible trees in newick format and 
    save in the designated folder.
    '''
    filenum = 0
    # assume that the foldername exist.
    for parasiteTree in generateTreeTuples(numLeaves, 'n', 0):
        for hostTree in generateTreeTuples(numLeaves, 'm', 0):
            parasiteLeaves = treeTupleLeaves(parasiteTree)
            hostLeaves = treeTupleLeaves(hostTree)
            for mappingPhi in generateMappings(parasiteLeaves, hostLeaves):
                # print it to a file instead of to screen
                filename = "test-size{}-no{}.newick".format(numLeaves, filenum)
                filepath = os.path.join(destFolderName, filename)
                with open(filepath, 'w') as f :
                    f.write(treeTupleStrings(hostTree) + ';\n')
                    f.write(treeTupleStrings(parasiteTree) + ';\n')
                    f.write(mappingToString(mappingPhi))
                filenum += 1
                
if __name__ == '__main__':
    for tree_size in [2, 3, 4, 5, 6]:
        targetFolder = './newickSample/size%d' % tree_size
        os.makedirs(targetFolder)
        generateNewickTests(tree_size, targetFolder)