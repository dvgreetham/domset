
import time
import csv
import math
import networkx as nx
import random
import os
import numpy as np
import scipy.optimize as sc





#Maximum degree for a graph
def max_delta(G):
    max_deg= 0
    for node in G.nodes():
        node_deg = G.degree(node)
        if node_deg > max_deg:
            max_deg = node_deg
        
    return max_deg

#Sum of klouts for the whole graph
def sum_weights(*args):
    
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    sum = 0
  
    for n in nodes:
        
        sum+=float(args[0].node[str(n)]['weight'])

    return sum





def solve_LP(adj_matrix,objfun,constraint,bounds1):
   
    res=sc.linprog(objfun, A_ub=adj_matrix, b_ub=constraint, bounds=bounds1, 
                   method='interior-point', options={'lstsq':True, 'presolve':True} )
    print(np.count_nonzero(res.x))
    if res.status ==2:
        print("Problem infeasible!")
    elif res.status==3:
        print("Problem unbounded!")
    elif res.status==1:
        print("Max iterations limit reached!")
    elif res.status==4:
        print("Serious numerical difficulties encountered!")
    return res.fun,res.x
   
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

#alpha=float(sys.argv[2])
#f1=sys.argv[1]
#G=nx.dense_gnm_random_graph(n, m, seed=None)
def rr(alpha, file, path):
    #f1="/Users/nv902387/Documents/Python27/graphs/list_random_networks//"+str(i)+"_random_network.gml" 
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
    #A=np.identity(n)
    A=np.zeros((n, n))
    for node, edata in G.nodes(data=True):
        objfun.append(edata['weight'])
        #print(edata)
        constraint.append(-math.ceil(alpha*G.degree(node)))
        #print(math.ceil(alpha*G.degree(node)))
        bounds1.append((0,1))
        for ngh in sorted(G.neighbors(node)):
        #print index[node], index[ngh]
            A[index[node], index[ngh]]=-1
    print(A.size, np.count_nonzero(A), len(objfun),np.count_nonzero(np.array(objfun)), len(bounds1))
    
  
    domsG=set()
    dom_setG=set()
    number_of_runs=math.ceil(math.log(max_delta(G), 2))
            
            
    violation=1
    nruns=0
            
    while (violation and nruns<number_of_runs):
        r=random.uniform(0, 0.5)
            #print r
        lp_result=solve_LP(A,objfun,constraint,bounds1)
        for i in range(0,len(lp_result[1])):
            if r<lp_result[1][i]:
                domsG.add(i)
    
       
        
        for key, value in index.items():
            if value in domsG:
                dom_setG.add(key)
    
        [viol, nds]=are_dominated(G, dom_setG, alpha)
        if viol==0: 
            violation=0
        else:
            nruns=nruns+1
            #print('Nruns', nruns)
        
                #if(number_of_runs==flag) and (violation==1):
                    ##print("Violation!!")
            
    print("Length of dom set is ", len(dom_setG)) 
    print("Total sum of weights is", sum_weights(G, dom_setG))      
            
            
  
    

    [viol, nodes]=are_dominated(G, dom_setG, alpha)
    if viol>0:
        violation=1
        l=[u for u in G.nodes(data=True) if u[0] not in dom_setG]
    
        weights=sorted(l, key=lambda x: x[1]['weight'])
     
  
        while (violation):
            [viol, nds]=are_dominated(G, dom_setG, alpha)
          
   
            if viol==0: 
                violation=0
            else:
                avgrate=math.ceil(viol/len(nds))
                for ii in range(avgrate):
                    dom_setG.add(weights.pop()[0])    
               
             
         
    
 
        print("Violations", are_dominated(G, dom_setG, alpha))
  
    
    print("Length of dominating set for alpha, ", alpha, " is ", len(dom_setG))
    print("Total sum of weights is", sum_weights(G, dom_setG))      
       
    t1=time.time() - t
    print("Time taken for LP...", t1)
    w = csv.writer(open(path, 'a'))
    w.writerow([file, alpha, len(dom_setG), sum_weights(G,dom_setG), t1])    
if __name__ == "__main__":

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, "RR-new.csv")
    print(path)
    if ~os.path.isfile(path):
        w = csv.writer(open(path, 'w'))
        w.writerow(["File", "Threshold", "DomSet-Size", "DomSet-Weights","Time taken"])
    for i in range(10):
        graph="pp5-100-0.18-0.0001-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rr(alpha, graph, path)   
    for i in range(10):
        graph="plc-500-10-0.8-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rr(alpha, graph, path)   
    for i in range(10):
        graph="er-500-5000-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            rr(alpha, graph, path)   
    for  graph in ["dsjc250-5-rw.gexf","r250-1-rw.gexf", "fpsol2-i-3-rw.gexf", "fb1-rw.gexf", "fb2-rw.gexf", "bitcoinalpha.gexf"]:
        for alpha in [0.25, 0.5, 0.75]:
            rr(alpha, graph, path)  
  
