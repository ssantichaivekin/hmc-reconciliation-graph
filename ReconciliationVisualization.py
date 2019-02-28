

import networkx as nx

data = {('p1', 'h1'): [['C', (None, None), (None, None), 1.0], 0],
 ('p2', 'h3'): [['C', (None, None), (None, None), 1.0], 0],
 ('p3', 'h2'): [['C', (None, None), (None, None), 1.0], 0],
 ('p4', 'h4'): [['C', (None, None), (None, None), 1.0], 0],
 ('p6', 'h7'): [['S', ('p7', 'h1'), ('p8', 'h2'), 0.5], 2],
 ('p6', 'h8'): [['S', ('p7', 'h3'), ('p8', 'h4'), 0.5], 2],
 ('p7', 'h1'): [['T', ('p1', 'h1'), ('p2', 'h3'), 0.5], 1],
 ('p7', 'h3'): [['T', ('p2', 'h3'), ('p1', 'h1'), 0.5], 1],
 ('p8', 'h2'): [['T', ('p3', 'h2'), ('p4', 'h4'), 0.5], 1],
 ('p8', 'h4'): [['T', ('p4', 'h4'), ('p3', 'h2'), 0.5], 1]}


H = {('h6', 'h8'): ('h6', 'h8', ('h8', 'h3'), ('h8', 'h4')), 
     ('h8', 'h3'): ('h8', 'h3', None, None), 
     ('h6', 'h7'): ('h6', 'h7', ('h7', 'h1'), ('h7', 'h2')), 
     
     ('h7', 'h2'): ('h7', 'h2', None, None), 
     ('h8', 'h4'): ('h8', 'h4', None, None), 
     ('h7', 'h1'): ('h7', 'h1', None, None),
     'hTop': ('Top', 'h20', ('h20', 'h6'), ('h20', 'h16')), 
     ('h16', 'h18'): ('h16', 'h18', ('h18', 'h13'), ('h18', 'h14')), 
     ('h18', 'h13'): ('h18', 'h13', None, None), 
     ('h16', 'h17'): ('h16', 'h17', ('h17', 'h11'), ('h17', 'h12')), 
      
     ('h17', 'h12'): ('h17', 'h12', None, None), 
     ('h18', 'h14'): ('h18', 'h14', None, None), 
     ('h17', 'h11'): ('h17', 'h11', None, None),
     ('h20', 'h6'): ('h20', 'h6', ('h6', 'h8'), ('h6', 'h7')), 
     ('h20', 'h16'): ('h20', 'h16', ('h16', 'h18'), ('h16', 'h17')) 
     }

P = {('p6', 'p8'): ('p6', 'p8', ('p8', 'p3'), ('p8', 'p4')), 
     ('p7', 'p2'): ('p7', 'p2', None, None), 
     ('p6', 'p7'): ('p6', 'p7', ('p7', 'p1'), ('p7', 'p2')), 
     ('p8', 'p4'): ('p8', 'p4', None, None), 
     ('p8', 'p3'): ('p8', 'p3', None, None), 
     'pTop': ('Top', 'p20', ('p20', 'p6'), ('p20', 'p16')), 
     ('p7', 'p1'): ('p7', 'p1', None, None),
     
     ('p16', 'p18'): ('p16', 'p18', ('p18', 'p13'), ('p18', 'p14')), 
     ('p17', 'p12'): ('p17', 'p12', None, None), 
     ('p16', 'p17'): ('p16', 'p17', ('p17', 'p11'), ('p17', 'p12')), 
     ('p18', 'p14'): ('p18', 'p14', None, None), 
     ('p18', 'p13'): ('p18', 'p13', None, None),  
     ('p17', 'p11'): ('p17', 'p11', None, None),
     ('p20', 'p6'): ('p20', 'p6', ('p6', 'p8'), ('p6', 'p7')), 
     ('p20', 'p16'): ('p20', 'p16', ('p16', 'p18'), ('p16', 'p17')) }


phi =  {'p2': 'h3', 'p3': 'h2', 'p1': 'h1', 'p4': 'h4', 
        'p12': 'h13', 'p13': 'h12', 'p11': 'h11', 'p14': 'h14'} 

def findMaxLevel(DTL):
    maxlevel = float('-inf')
    for key, ele in DTL.items():
        if ele[-1] > maxlevel:
            maxlevel = ele[-1]           
    return maxlevel

def organize(DTL):
    level = findMaxLevel(DTL)
    output = []
    while level >= 0:
        for key, ele in DTL.items():
            if ele[-1] == level:
                eventToMap = (key, tuple(ele[0]))
                output.append(eventToMap)
                if level != 0:
                    mapToEvent1 = (tuple(ele[0]), ele[0][1])
                    mapToEvent2 = (tuple(ele[0]), ele[0][2])
                    output.append(mapToEvent1)
                    output.append(mapToEvent2)
        level -= 1
    return output
        


def graph(DTL):
    g=nx.DiGraph()
    g.add_edges_from(organize(DTL))
    p=nx.drawing.nx_pydot.to_pydot(g)
    p.write_png('DTL8.png')
                    

H = {('h6', 'h8'): ('h6', 'h8', ('h8', 'h3'), ('h8', 'h4')), ('h8', 'h3'): ('h8', 'h3', None, None), ('h6', 'h7'): ('h6', 'h7', ('h7', 'h1'), ('h7', 'h2')), 'hTop': ('Top', 'h6', ('h6', 'h7'), ('h6', 'h8')), ('h7', 'h2'): ('h7', 'h2', None, None), ('h8', 'h4'): ('h8', 'h4', None, None), ('h7', 'h1'): ('h7', 'h1', None, None)}
P = {('p6', 'p8'): ('p6', 'p8', ('p8', 'p3'), ('p8', 'p4')), ('p7', 'p2'): ('p7', 'p2', None, None), ('p6', 'p7'): ('p6', 'p7', ('p7', 'p1'), ('p7', 'p2')), ('p8', 'p4'): ('p8', 'p4', None, None), ('p8', 'p3'): ('p8', 'p3', None, None), 'pTop': ('Top', 'p6', ('p6', 'p7'), ('p6', 'p8')), ('p7', 'p1'): ('p7', 'p1', None, None)}
phi =  {'p3': 'h3', 'p2': 'h4', 'p4': 'h1', 'p1': 'h2'} 