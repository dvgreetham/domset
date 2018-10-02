#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:24:27 2018

@author: dvg27
"""

import math
import networkx as nx
import operator

def sum_klouts(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    sum = 0
    for n in nodes:
        sum+=int(args[0].node[n]['weight'])

    return sum


def cost_fun(G, doms):
    cost_fun={}
    for node in G.nodes():
        if node in doms:
            cf=0
        else: 
            if G.node[node]['weight']==0:
                ww=0.1
            else:
                ww=G.node[node]['weight']
            cf=int(G.degree[node]/ww)
           # cf=int(G.node[node]['weight'])/G.degree[node]
            
            for u in G.neighbors(node):
                if u in doms:
                    s=1
                else: 
                    if G.node[node]['weight']==0:
                        ww=0.1
                    else:
                        ww=G.node[node]['weight']
                    
                    s=0
                    cf=cf+(1-s)*int(G.degree[u]/ww)
                    #cf=cf+(1-s)*int(G.node[u]['weight'])/G.degree[u]
        cost_fun[node]=cf
    return cost_fun

def is_dominated(G, doms, alpha):
    for u in G.nodes():    
        intersect = [ii for ii in doms if ii in G.neighbors(u)]                       
        rate=len(intersect)
        if rate < math.ceil(alpha*G.degree(u)):
            return  math.ceil(alpha*G.degree(u))-rate
    return 0

# main 
for alpha in [0.25, 0.5, 0.75]:
    G=nx.read_gexf("lp2.gexf")
     
    doms=[]
    costfun=cost_fun(G, doms)
    r=is_dominated(G, doms, alpha)
    while r>0:
        l=[u for u in G.nodes() if u not in doms]
       
        rel_costfun={k: costfun[k] for k in l} #calculate cost only for relevant nodes, not in D
        
        sorted_rel_costfun = sorted(rel_costfun.items(), key=operator.itemgetter(1), reverse=True)
        # add the first r nodes that maximise the cost function
        for i in range(r):
            doms.append(sorted_rel_costfun[i][0])
        costfun=cost_fun(G, doms)
        r=is_dominated(G, doms, alpha)
    print("Length of dominating set for alpha", alpha, "is", len(doms))
    print("Total sum of weights is", sum_klouts(G, doms))   
    
        