#This is a modified implementation of the algorithm for finding alpha-rate dominating sets. It has a random approach but at the same time it
#takes into account klouts. The
# output is an alpha rate dominating set with minimised klouts.

# Start with alpha_rate_B_preferential.py N, M, pro, alpha    where N is number of vertices,
# M number of edges added in each iteration (1 or 10)  for preferntial attachment, pro probability of triangles 0.8
# and alpha   0.25, 0.5 and 0.75
import time 
import operator
import scipy.optimize as sc
#from lpsolve55 import *
#from lp_maker import *
#import sys
import csv
import math
import networkx as nx
#import matplotlib.pyplot as plt
#import operator as op
#import functools
import random
#import logging
import os
#from random import randint
import numpy as np
import community


#Minimum degree for a graph
def min_delta(G, alpha):
    it = iter(G.nodes())
    first = next(it)
    min_deg= G.degree(first)
    for node in it:
        node_deg = G.degree(node)
        if node_deg < min_deg or min_deg==0:
            min_deg = node_deg
    print(min_deg)
    delta = math.floor(min_deg*(1-alpha)) + 1

    return delta

#Maximum degree for a graph
def max_delta(G):
    max_deg= 0
    for node in G.nodes():
        node_deg = G.degree(node)
        if node_deg > max_deg:
            max_deg = node_deg
        
    return max_deg

#Sum of klouts for the whole graph
def sum_klouts(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    sum = 0
    for n in nodes:
        sum+=int(args[0].node[n]['weight'])

    return sum

#Max of the klouts the whole graph
def max_klouts(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    max_klout = 0
    for n in nodes:
        if int(args[0].node[n]['weight']) > max_klout:
            max_klout = int(args[0].node[n]['weight'])

    return max_klout

#Min of the klouts the whole graph
def min_klouts(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    it = iter(nodes)
    first = next(it)
    min_klout = int(args[0].node[first]['weight'])
    for n in it:
        if int(args[0].node[n]['weight']) < min_klout:
            min_klout = int(args[0].node[n]['weight'])
       
    return min_klout




# Return sum of klout values in N[v].
def sum_klouts_neighbourhood(G, neighbours):
    sum = 0
    for nb in neighbours:
        sum+=int(G.node[nb]['weight'])

    return sum


#returns an alpha rate dominating set of G



def solve_LP(adj_matrix,objfun,constraint,bounds1):
    #print(np.count_nonzero(A))
    res=sc.linprog(objfun, A_ub=adj_matrix, b_ub=constraint, bounds=bounds1, method='interior-point' )
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



    #alpha=float(sys.argv[2])
#f1=sys.argv[1]
for i in range(1):
    #f1="/Users/dvg27/Documents/Documents/Python27/graphs/list_random_partition_network_p_in0.02p_out0.0001/0_rp.graphml" 
  
    G=nx.read_gexf("lp2.gexf")
    
    if (nx.number_connected_components(G)>1):
        print("G is disconnected")
        exit(1)
 
    index={}
    i=0
    for node in G.nodes():
        index[node]=i
        i=i+1
    
    partition = community.best_partition(G)
    size=float(len(set(partition.values())))
    print("Size_of_communities", size)
    
    for alpha in [0.25, 0.50, 0.75]:    
        time_taken=[]
        domset_total=[]
        for com in set(partition.values()) :
                objfun=[]
                constraint=[]
                #all_one=[]
                bounds1=[]
                index1={}
                list_nodes = [nodes for nodes in partition.keys()
                                            if partition[nodes] == com]
               # G1=nx.Graph(G.subgraph(list_nodes)) 
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
            #print(math.ceil(alpha*G.degree(node)))
                    bounds1.append((0,1))
                    for ngh in sorted(G1.neighbors(node)):
            #print index[node], index[ngh]
                        A[index1[node], index1[ngh]]=-1
                doms=[]
                min_sum = sum_klouts(G1)
                dom_set_array = {}
                sum_klouts_dom = []
                max_sum = 0
                viol_num=0
                average_sum = 0
                min_dom_size = G1.number_of_nodes()
                max_dom_size = 0
                avg_cardinality = 0
                number_of_runs=math.ceil(math.log(max_delta(G1), 2))
                #print(number_of_runs)
                #for i in range(0,100):
                t=time.time()
                lp_result=solve_LP(A,objfun,constraint,bounds1)
                t1=time.time() - t
                #print("Time taken for LP...", t1)
                for iii in range(0, 100):
                    #print("******", iii, "*******")
                    t=time.time()
                    doms=[]
                    viol_num=0
                    average_sum = 0
                    B=[]
                    viol_num=0
                    violation=1
                    nruns=0
                    while (violation and nruns<number_of_runs):
                        viol_num=0
                        r=random.uniform(0, 0.5)
                        #print r
                        for i in range(0,len(lp_result[1])):
                            if r<lp_result[1][i]:
                                B.append(i)
                
                        dom_set=list(set(B))
                        #print dom_set
                        #    for i in dom_set:
                        #    print 'node', G.nodes(data=True)[i]
                        for node in G1.nodes():
                                set_N=[]
                                for jj in G1.neighbors(node):
                                    set_N.append(index1[jj])
                                intersect = [ii for ii in dom_set if ii in set_N]
                                #print len(target)            
                               
                                rate=len(intersect)
                               
                                if  (len(set_N)!=0) and (rate < math.ceil(alpha*G1.degree(node))):
                                    #print "violation at ", index[node]
                                    #print rate, math.ceil(alpha*G.degree(node)), math.ceil(alpha*G.degree(node))-rate
                                    viol_num+=1
                                    break
                        #print "Viol_num", viol_num    
                        if viol_num==0:
                            violation=0
                        nruns+=1
                    #if(number_of_runs==flag) and (violation==1):
                        ##print("Violation!!")
           
                    len_dom_set = len(dom_set)
                    #print("Length of dom set is ", len(dom_set)) 
                    avg_cardinality+=len_dom_set
                    doms=[key for key, item in index1.items() if item in dom_set]
                    sum_klout = sum_klouts(G,doms)
                    dom_set_array[sum_klout]=dom_set
                    sum_klouts_dom.append(sum_klout)
                    if sum_klout < min_sum:
                        min_sum = sum_klout
                    if sum_klout > max_sum:
                            max_sum = sum_klout
                    if len_dom_set < min_dom_size:
                        min_dom_size = len_dom_set
                    if len_dom_set > max_dom_size:
                        max_dom_size = len_dom_set
                    average_sum+=sum_klout
                    t1=time.time() - t
                    #print("Computing randomised alpha rate:", t1)
                    time_taken.append(t1)
                avg_cardinality /= 100
                sum_klouts_dom.sort()
                #dom_set_array.sort()
                
                doms=[key for key, item in index1.items() if item in dom_set_array[min_sum]]   
                #print("Doms", len(doms))
                domset_total.extend(doms) 
        
        t=time.time()
        sum_klout = sum_klouts(G,domset_total)
        
        
        #print("Total lenght", len(domset_total), sum_klout)
        domset_total=sorted(list(set(domset_total)))
        for node in G.nodes():
            set_N=[]
            for jj in G.neighbors(node):
                set_N.append(index[jj])
            intersect = [ii for ii in domset_total if ii in set_N]
                
                #print intersect           
                       
            rate=len(intersect)
            #if (rate>1):
                #print intersect
            r= math.ceil(alpha*G.degree(node))-rate
            if  (len(set_N)!=0) and (r>0):
                #print("violation at ", index[node])
                #print(rate, math.ceil(alpha*G.degree(node)), math.ceil(alpha*G.degree(node))-rate, len(domset_total))
                viol_num+=1
              
                patchweights={}
                diff = [ii for ii in set_N if ii not in domset_total]  
                #print("Diff", len(diff))
                for node, edata in G.nodes(data=True):
                    if(node in diff):
                        #node, edata['weight']
                        patchweights[node]=edata['weight']
                ii=0
                xs=sorted(patchweights.items(), key=operator.itemgetter(1))
                for key, item in xs:
                    if (ii>=int(r)):
                        break 
                    domset_total.append(key)
                    #print key, patchweights[key]
                    
                    ii+=1
                #print(len(domset_total))
                #exit(1);
        #print("Violations", viol_num)
        viol_num=0
        #print("Patching finished! Time taken",float(len(time_taken))) 
        for node in G.nodes():
            set_N=[]
            for jj in G.neighbors(node):
                set_N.append(index[jj])
            intersect = [ii for ii in domset_total if ii in set_N]
                
                #print intersect           
                       
            rate=len(intersect)
            #if (rate>1):
                #print intersect
                 
            if  (len(set_N)!=0) and (rate < math.ceil(alpha*G.degree(node))):
               # print("violation at ", index[node])
                #print(rate, math.ceil(alpha*G.degree(node)), math.ceil(alpha*G.degree(node))-rate)
                viol_num+=1
                #break
        #print("Violations", viol_num, len(domset_total))
        t1=time.time() - t
        new_file = "B_comms.csv"
        script_dir = os.path.dirname(os.path.abspath(__file__))
        #print(script_dir)
        new_file = True
        path = os.path.join(script_dir, "AlgRRcomms.csv")
        if os.path.isfile(path):
            new_file = False
        if new_file:
            w = csv.writer(open(path, 'w'))
            #w = csv.writer(open(path, 'w', encoding="utf8"))
            w.writerow(["File", "Threshold", "DominatingSet-Size", "DominatingSet-Weights","AverageCardinalityDominatingSet", "AverageSumOFWeightsDomSet","MinCardinality", "MaxCardinality", 
                       "MinSumOfWeights", "MaxSumOfWeights","Average time"])
        else:
            w = csv.writer(open(path, 'a'))
            #w = csv.writer(open(path, 'a', encoding="utf8"))
        w.writerow(['file', alpha, len(domset_total), sum_klouts(G,domset_total), avg_cardinality, average_sum, min_dom_size, max_dom_size, min_sum, max_sum,t1+sum(time_taken)/float(len(time_taken))])
        print("Length of dominating set for alpha, ", alpha, " is ", len(domset_total))
        print("Total sum of weights is", sum_klouts(G, domset_total))      
           
           
        
        #plt.plot(sum_klouts_dom, dom_set_array, color='b', lw=1.5, zorder=1)
        #plt.scatter(sum_klouts_dom, dom_set_array,color='r', zorder=2)
        #print([min_sum-100, max_sum+100, min_dom_size-100, max_dom_size+100])
        #plt.suptitle('Algorithm B Results ')
        #plt.axis([min_sum-100, max_sum+100, min_dom_size-50, max_dom_size+50])
        #plt.xlabel("Sum of weights in dominating set")
        #plt.ylabel("Cardinality of dominating set")
        #print(os.path.splitext(sys.argv[1])[0])
        #plt.savefig(os.path.splitext(sys.argv[1])[0]+"_pro_" +str(pro)+"_alpha_"+str(alpha)+"PrefAlgB-communties.png")
        
        #plt.show()
        #try: # draw
        #    pos=nx.spring_layout(G,iterations=10)
        #    nx.draw(G,pos,node_size=0,alpha=0.4,edge_color='r',font_size=6)
        #    plt.savefig("twitternetwork.png")
        #    plt.show()
        #except: # matplotlib not available
        #    pass
