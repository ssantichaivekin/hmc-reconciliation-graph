'''
Generating Tests for Reconciliation Graphs
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
    if type(treeTuple) is str :
        return treeTuple
    else :
        return "({},{}){}".format(
            treeTupleStrings(treeTuple[1]), # left
            treeTupleStrings(treeTuple[2]), # right
            treeTuple[0] # root
            )

def treeTupleLeaves(treeTuple):
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
    os.makedirs('./newickSample/size2')
    generateNewickTests(2, './newickSample/size2')
    os.makedirs('./newickSample/size3')
    generateNewickTests(3, './newickSample/size3')
    os.makedirs('./newickSample/size4')
    generateNewickTests(4, './newickSample/size4')
    os.makedirs('./newickSample/size5')
    generateNewickTests(5, './newickSample/size5')