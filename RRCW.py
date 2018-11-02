
import time 
import operator
import scipy.optimize as sc

import os
import csv
import math
import networkx as nx

import random

import numpy as np
import community





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
def solve_LP(adj_matrix,objfun,constraint,bounds1):
    #print(np.count_nonzero(A))
    res=sc.linprog(objfun, A_ub=adj_matrix, b_ub=constraint, bounds=bounds1, 
                   method='interior-point')# options={'lstsq':True, 'presolve':True} )
    #print(np.count_nonzero(res.x))
    if res.status ==2:
        print("Problem infeasible!")
    elif res.status==3:
        print("Problem unbounded!")
    elif res.status==1:
        print("Max iterations limit reached!")
    elif res.status==4:
        print("Serious numerical difficulties encountered!")
    return res.fun,res.x

# rrcw algoritm returns a low-weight dominating set total_domset for G, such that 
# each node in G has at least alpha% of neighbours in domset_total
def rrcw(alpha, file, path):
    
    G=nx.read_gexf(file)
    #G=nx.read_gexf("bitcoinalpha.gexf")
    if (nx.number_connected_components(G)>1):
        print("G is disconnected!!!")
        exit(1)
    t=time.time()
    index={}
    i=0
    for node in G.nodes():
        index[node]=i
        i=i+1
    
    partition = community.best_partition(G)
    size=float(len(set(partition.values())))
    print("Size_of_communities", size)
    
    domset_total=set()
    for com in set(partition.values()) :
            objfun=[]
            constraint=[]
            bounds1=[]
            index1={}
            list_nodes = [nodes for nodes in partition.keys()
                                        if partition[nodes] == com]
            G1=G.subgraph(list_nodes)
            n=G1.number_of_nodes()
            i=0
            for node in G1.nodes():
                index1[node]=i
                i=i+1
            A=np.zeros((n, n))
            for node, edata in G1.nodes(data=True):
                objfun.append(edata['weight'])    
                constraint.append(-math.ceil(alpha*G1.degree(node)))
                bounds1.append((0,1))
                for ngh in sorted(G1.neighbors(node)):
                    A[index1[node], index1[ngh]]=-1
            domsG1=set()
            dom_setG1=set()
            number_of_runs=math.ceil(math.log(max_delta(G1), 2))
            
            violation=1
            nruns=0
            
            while (violation and nruns<number_of_runs):
                
                r=random.uniform(0, 0.5)
                    #print r
                lp_result=solve_LP(A,objfun,constraint,bounds1)
                for i in range(0,len(lp_result[1])):
                    if r<lp_result[1][i]:
                        domsG1.add(i)
            
               
                
                for key, value in index1.items():
                    if value in domsG1:
                        dom_setG1.add(key)
            
                [viol, nds]=are_dominated(G1, dom_setG1, alpha)
                if viol==0: 
                    violation=0
                else:
                    nruns=nruns+1
             
            print("Length of dom set in ", com, " is ", len(dom_setG1)) 
            print("Total sum of weights is", sum_weights(G1, dom_setG1))      
            
            
   
            domset_total.update(dom_setG1)
           
            
    [viol, nodes]=are_dominated(G, domset_total, alpha)
    if viol>0:
        violation=1
        l=[]
        for u in G.nodes(data=True):
            if u[0] not in domset_total:
                l.append(u)
    
        weights=sorted(l, key=lambda x: x[1]['weight'])
     
  
        while (violation):
            [viol, nds]=are_dominated(G, domset_total, alpha)
          
   
            if viol==0: 
                violation=0
            else:
                avgrate=math.ceil(viol/len(nds))
                for ii in range(avgrate):
                    domset_total.add(weights.pop()[0])    
               
             
    print("Violations", are_dominated(G, domset_total, alpha))
    print("Length of dominating set for alpha, ", alpha, " is ", len(domset_total))
    print("Total sum of weights is", sum_weights(G, domset_total))      
       
    t1=time.time() - t
    print("Time taken for LP...", t1)
  
    w = csv.writer(open(path, 'a'))
    w.writerow([file, alpha, len(domset_total), sum_weights(G,domset_total), t1])      
    #print(domskeys)            
if __name__ == "__main__":

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, "RRCW-new1.csv")
    print(path)
    if ~os.path.isfile(path):
        w = csv.writer(open(path, 'w'))
        w.writerow(["File", "Threshold", "DomSet-Size", "DomSet-Weights","Time taken"])
    for i in range(10):
        graph="pp5-100-0.18-0.0001-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rrcw(alpha, graph, path)   
    for i in range(10):
        graph="plc-500-10-0.8-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rrcw(alpha, graph, path)   
    for i in range(10):
        graph="er-500-5000-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rrcw(alpha, graph, path)   
    for  graph in ["dsjc250-5-rw.gexf","r250-1-rw.gexf", "fpsol2-i-3-rw.gexf", "fb1-rw.gexf", "fb2-rw.gexf", "bitcoinalpha.gexf"]:
        for alpha in [0.25, 0.5, 0.75]:
            rrcw(alpha, graph, path)  
  