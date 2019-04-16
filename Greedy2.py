#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:24:27 2018

@author: dvg27
"""

import math
import networkx as nx
import csv
import os
import operator
import time
def sum_weights(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    sum = 0
    for n in nodes:
        sum+=float(args[0].node[n]['weight'])

    return sum


def cost_fun(G, nodes, doms, domd):
    cost_fun={}
    for node in nodes:
        ww=G.node[node]['weight']
        sum_of_neighbors_weights=0
        for u in G.neighbors(node):
            if u not in doms and  u not in domd:
                 sum_of_neighbors_weights=sum_of_neighbors_weights+G.node[u]['weight']
            
        cf=ww/max(0.1,sum_of_neighbors_weights)
                    #cf=cf+(1-s)*int(G.node[u]['weight'])/G.degree[u]
        cost_fun[node]=cf
    return cost_fun

def is_dominated(G, nodes,  doms, domd, alpha):
    for u in nodes:    
        intersect = [ii for ii in G.neighbors(u) if ii in doms ]                       
        rate=len(intersect)
        if rate < math.ceil(alpha*G.degree(u)):
            return  math.ceil(alpha*G.degree(u))-rate
        else: domd.add(u)
    return 0


def greedy3(alpha, file, path):

    t=time.time()
    #G=nx.read_gexf("lp1.gexf")
    G=nx.read_gexf(file)
    doms=set()
    domd=set()
    nodes=G.nodes()
    costfun=cost_fun(G, nodes, doms, domd)
    r=is_dominated(G, nodes, doms, domd, alpha)
    while r>0:
        sorted_rel_costfun = sorted(costfun.items(), key=operator.itemgetter(1))
        # add the first r nodes that maximise the cost function   
        
        for i in range(r):
            doms.add(sorted_rel_costfun[i][0])
        
        nodes=list(set(nodes)-doms)
        #nodes=list(set(nodes)-domd)
    
        costfun=cost_fun(G, nodes, doms, domd)
        r=is_dominated(G, G.nodes(), doms, domd, alpha)
    t1=time.time() - t
    sorted_rel_costfun = sorted(costfun.items(), key=operator.itemgetter(1))
    
    print("Time taken ...", t1)
    print("Length of dominating set  for alpha ", alpha, " is", len(doms))
    print("Total sum of weights is", sum_weights(G, doms))        
                
    
    w = csv.writer(open(path, 'a'))
    w.writerow([file, alpha, len(doms), sum_weights(G,doms), t1])    
    
if __name__ == "__main__":

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, "Greedy2-new11.csv")
    print(path)
    if os.path.isfile(path)==False:
        w = csv.writer(open(path, 'w'))
        w.writerow(["File", "Threshold", "DomSet-Size", "DomSet-Weights","Time taken"])

    for i in range(10):
        graph="pp5-100-0.18-0.0001-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            greedy3(alpha, graph, path)   
    for i in range(10):
        graph="plc-500-10-0.8-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            greedy3(alpha, graph, path)   
    for i in range(10):
        graph="er-500-5000-"+str(i)+".gexf"
        for alpha in [0.25, 0.5, 0.75]:
            greedy3(alpha, graph, path)   
   
    for  graph in [ "dsjc250-5-rw.gexf","r250-1-rw.gexf", "fpsol2-i-3-rw.gexf", "fb1-rw.gexf", "fb2-rw.gexf", "bitcoinalpha.gexf"]:
        for alpha in [0.25, 0.5, 0.75]:
            greedy3(alpha, graph, path)  
   