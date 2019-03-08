'''
Testing visualization
'''

import DTLReconGraph
import ReconciliationVisualization

result = DTLReconGraph.reconcile("./newickSample/size5/test-size5-no700.newick", 2, 4, 2)
import pprint
pp = pprint.PrettyPrinter()
pp.pprint(result)
ReconciliationVisualization.visualizeAndSave(result[2], './sampleVis700.png')

