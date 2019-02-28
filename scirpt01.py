import DTLReconGraph
import ReconciliationVisualization

result = DTLReconGraph.reconcile("./TreeLifeData/COG0001.newick", 4, 4, 1)
import pprint
pp = pprint.PrettyPrinter()
pp.pprint(result)
ReconciliationVisualization.visualizeAndSave(result[2])

