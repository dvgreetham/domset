from __future__ import print_function
import time 
import scipy.optimize as sc
import os
import csv
import math
import networkx as nx

import random
#import cvxpy as cp
import numpy as np
from gurobipy import *




#Maximum degree for a graph
def max_delta(G):
    max_deg= 0
    for node in G.nodes():
        node_deg = G.degree(node)
        if node_deg > max_deg:
            max_deg = node_deg
        
    return max_deg

#Sum of weights of nodes nds in a graph  G
def sum_weights(G, nds):
   
    sum = 0
    for n in nds:
        sum+=float(G.node[n]['weight'])

    return sum


#returns viol, total no of neighbours that need to be dominated and 
#non-dominated nodes in a graph G, with the current dominating set doms 
def is_dominated(G, nodes,  doms, domd, alpha):
    for u in nodes:    
        intersect = [ii for ii in G.neighbors(u) if ii in doms ]                       
        rate=len(intersect)
        if rate < math.ceil(alpha*G.degree(u)):
            return  math.ceil(alpha*G.degree(u))-rate
        else: domd.add(u)
    return 0

def are_dominated(G, doms, alpha):
    viol=0
    nodes1=[]
    for u in G.nodes():
        intersect=[]
        difference=[]
        for ii in G.neighbors(u):
            if ii in doms:
                intersect.append(ii) 
            else:
                difference.append(ii)                       
        rate=len(intersect)
        if rate < math.ceil(alpha*G.degree(u)):
            r=math.ceil(alpha*G.degree(u))-rate
            viol+=r
            nodes1.extend(np.array(difference)[:r])
    return viol, nodes1


#solves LP min objfun*x, s.t adj_matrix *x>=b

def solve_ilp(A, objfun, constraint):
    n=len(objfun)
    indices=[]
    try:
        x=[]
        for i in range(n):
            x.append('x'+str(i))
        obj=0
       
    # Create a new model
        m = Model("mip1")
    
        # Create variables
        for i in range(n):
            x[i] = m.addVar(vtype=GRB.BINARY, name="x"+str(i))
            obj+=objfun[i]*x[i]
    
        # Set objective
        m.setObjective(obj, GRB.MINIMIZE)
    
        # Add constraint: x + 2 y + 3 z <= 4
        for i in range(n):
            l=[]
            for j in range(n):
                l+=A[i, j]*x[j]
            m.addConstr(l >= constraint[i], "c"+str(i))
    
        # Add constraint: x + y >= 1
    
    
        m.optimize()
    
        for v in m.getVars():
            print('%s %g' % (v.varName, v.x))
            if v.x==1: 
                ind=int(v.varname[1:])
                print(ind)
                indices.append(ind)
        print('Obj: %g' % m.objVal)
        return (m.objVal, indices)
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')

    

def ilp(alpha, file, path):
    
    G=nx.read_gexf(file)
   
    #=nx.read_gexf("bitcoinalpha.gexf")
    if (nx.number_connected_components(G)>1):
        print("G is disconnected")
        exit(1)
 
    n=G.number_of_nodes()
    m=G.number_of_edges()
    print(n,m)
    index={}
    
    i=0
    for node in G.nodes():
        index[node]=i
        i=i+1
    
   
    t=time.time()
    objfun=[]
    constraint=[]
    bounds1=[]

   # x= cp.Variable(n, boolean=True)
   
 #  A1=cp.Variable((n,n))
    #A=np.identity(n)
    A=np.zeros((n, n))
    for node, edata in G.nodes(data=True):
        objfun.append(edata['weight'])
        #print(edata)
        constraint.append(math.ceil(alpha*G.degree(node)))

        #print(math.ceil(alpha*G.degree(node)))
        bounds1.append((0,1))
        for ngh in sorted(G.neighbors(node)):
        #print index[node], index[ngh]
            A[index[node], index[ngh]]=1
    [objval, indices]=solve_ilp(A, objfun, constraint)
    dom_setG=set()
    for node in G.nodes():
        if index[node] in indices:
            dom_setG.add(node)
  
   
 
    print("Violations", is_dominated(G, dom_setG, alpha))
  
    
    print("Length of dominating set for alpha, ", alpha, " is ", len(dom_setG))
    print("Total sum of weights is", sum_weights(G, dom_setG), objval)      
       
    t1=time.time() - t
    print("Time taken for LP...", t1)
    w = csv.writer(open(path, 'a'))
    w.writerow([file, alpha, len(dom_setG), objval, t1])    
if __name__ == "__main__":

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, "ILP.csv")
    print(path)
    if ~os.path.isfile(path):
        w = csv.writer(open(path, 'w'))
        w.writerow(["File", "Threshold", "DomSet-Size", "DomSet-Weights","Time taken"])
#    for i in range(10):
#        graph="pp5-100-0.18-0.0001-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:
#            ilp(alpha, graph, path)   
#    for i in range(10):
#        graph="plc-500-10-0.8-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:
#            ilp(alpha, graph, path)   
#    for i in range(10):
#        graph="er-500-5000-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:        
#            ilp(alpha, graph, path)   
#    H=nx.Graph()
#    for i in range(5):
#        H.add_node(i, weight=4)
#    for i in range(1,5):
#        H.add_edge(i,0)
#    H.add_edge(3,4)
#    nx.write_gexf(H, "p1.gexf")
#    
    for  graph in ["dsjc250-5-rw.gexf", "r250-1-rw.gexf", "fpsol2-i-3-rw.gexf"]:# 
   # for  graph in ["fb1-rw.gexf", "fb2-rw.gexf", "bitcoinalpha.gexf"]:
        for alpha in [0.25, 0.5, 0.75]:
            ilp(alpha, graph, path)  
 
