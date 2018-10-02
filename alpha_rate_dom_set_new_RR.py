#This is an implementation of the original randomised algorithm for finding alpha-rate dominating sets.
# run with n, m, alpha
import time
import csv
import math
import networkx as nx
import random
import os
import numpy as np
import scipy.optimize as sc
#LOG_FILENAME = 'AlgRR.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)


def min_delta(G, alpha):
    it = iter(G.nodes())
    first = next(it)
    min_deg= G.degree(first)
    for node in it:
        node_deg = G.degree(node)
        if node_deg < min_deg or min_deg==0:
            min_deg = node_deg
    #print(min_deg)
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
        
        sum+=int(args[0].node[str(n)]['weight'])

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
    print(first)
    min_klout = int(args[0].node[first]['weight'])
    for n in it:
        if int(args[0].node[n]['weight']) < min_klout:
            min_klout = int(args[0].node[n]['weight'])
       
    return min_klout



# N[v] return list of nodes that have an edge from vertex v.
# Again the vertices that have no outgoing edges areignored.
def vertex_neighbourhood(G, node):
    neighb = G.neighbors(node)
    neighb.append(node)

    return neighb



# Return sum of klout values in N[v].
def sum_klouts_neighbourhood(neighbours):
    sum = 0 	
    for n in neighbours:
        sum+=n['weight']    
    return sum

#returns an alpha rate dominating set of G



def solve_LP(adj_matrix,objfun,constraint,bounds1):
    print(np.count_nonzero(A))
    res=sc.linprog(objfun, A_ub=adj_matrix, b_ub=constraint, bounds=bounds1, method='interior-point' )
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
   


#alpha=float(sys.argv[2])
#f1=sys.argv[1]
#G=nx.dense_gnm_random_graph(n, m, seed=None)
for i in range(1):
    #f1="/Users/nv902387/Documents/Python27/graphs/list_random_networks//"+str(i)+"_random_network.gml" 
    G=nx.read_gexf("lp0.gexf")

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
    
    for alpha in [0.25, 0.5, 0.75]:
   
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
        
        min_sum = sum_klouts(G)
        dom_set_array = []
        sum_klouts_dom = []
        max_sum = 0
        average_sum = 0
        min_dom_size = G.number_of_nodes()
        max_dom_size = 0
        avg_cardinality = 0
        time_taken=[]
        t=time.time()
        lp_result=solve_LP(A, objfun,constraint,bounds1)
        tlp=time.time() - t
        t1=0
        print("Time taken for LP...", tlp, "objfun", lp_result[0])
       # for j in range(0, n):
        #    if lp_result[1][j]>0:
         #       print(j,lp_result[1][j])
        for i in range(0,100):
            #print("******", i, "*******")
            t=time.time()
            dom_set=[]
            B=[]
            viol_num=0
            violation=1
            flag=0
            number_of_runs=math.ceil(math.log(max_delta(G), 2))
            #print(number_of_runs)
            while (violation and flag<number_of_runs):
                viol_num=0
                r=random.uniform(0, 0.5)
                for i in range(0,n):
                    if r<lp_result[1][i]:
                        B.append(i)
                dom_set=list(set(B))
                for node in G.nodes():
                    set_N=[]
                    for jj in G.neighbors(node):
                        set_N.append(index[jj])            
                   
                    intersect = [ii for ii in dom_set if ii in set_N]
                     
                    rate=len(intersect)
                
                    if  (len(set_N)!=0) and (rate < alpha*G.degree(node)):
                        #print("Violation at ", index[node])
                        #print(rate, math.ceil(alpha*G.degree(node)), math.ceil(alpha*G.degree(node))-rate)
                        #print "violation!!!!"
                        viol_num+=1
                        break
                if viol_num==0:
                    violation=0
                flag=flag+1
            #if(number_of_runs==flag) and (violation==1):
                #print "Violation!!"
            
            #print("Computing randomised alpha rate:", t1)
            time_taken.append(t1)
            len_dom_set = len(dom_set)
            avg_cardinality+=len_dom_set
            #dom_set_array.append(len_dom_set)
            sum_klout = sum_klouts(G,dom_set)
            #sum_klouts_dom.append(sum_klout)
            if sum_klout < min_sum:
                min_sum = sum_klout
            if sum_klout > max_sum:
                    max_sum = sum_klout
            if len_dom_set < min_dom_size:
                min_dom_size = len_dom_set
            if len_dom_set > max_dom_size:
                max_dom_size = len_dom_set
            average_sum+=sum_klout
            t1=time.time()-t
        #    logging.info("\tthe size of alpha rate dominating set is: %d", len_dom_set)
        #    logging.info("\tmaximum of klouts in dominating set is: %d ",max_klouts(G,dom_set))
        #    logging.info("\tminimum of klouts in dominating set is: %d ",min_klouts(G,dom_set))
        #    logging.info("\tsum of klouts in dominating set is: %d ",sum_klout)
        #    logging.info("\taverage of klouts in dominating set is: %d ",sum_klout/len_dom_set)
        #    logging.info("\r\n")
        #
        #    logging.info("\r\n")
            for node in G.nodes():
                set_N=[]
                for jj in G.neighbors(node):
                    set_N.append(index[jj])
                intersect = [ii for ii in dom_set if ii in set_N]
                    
                    #print intersect           
                           
                rate=len(intersect)
                #if (rate>1):
                    #print intersect
                     
                if  (len(set_N)!=0) and (rate < math.ceil(alpha*G.degree(node))):
                    #print("violation at ", index[node])
                    #print(rate, math.ceil(alpha*G.degree(node)), math.ceil(alpha*G.degree(node))-rate)
                    viol_num+=1
                    #break
        print("Violations", viol_num, len(dom_set))
        
        print("Time taken", float(len(time_taken)))
        average_sum /=100.
        avg_cardinality /= 100.
        #sum_klouts_dom.sort()
        #dom_set_array.sort()
        
        new_file = "AlgRR.csv"
        script_dir = os.path.dirname(os.path.abspath(__file__))
        print(script_dir)
        new_file = True
        path = os.path.join(script_dir, "AlgRR.csv")
        if os.path.isfile(path):
            new_file = False
        if new_file:
            w = csv.writer(open(path, 'w'))
            #w = csv.writer(open(path, 'w', encoding="utf8"))
            w.writerow(["File", "Threshold", "AverageCardinalityDominatingSet", "AverageSumOFWeightsDomSet","MinCardinality", "MaxCardinality", 
                       "MinSumOfWeights", "MaxSumOfWeights","Average time", "LP time"])
        else:
            w = csv.writer(open(path, 'a'))
            #w = csv.writer(open(path, 'a', encoding="utf8"))
        w.writerow(["file", alpha, avg_cardinality, average_sum, min_dom_size, max_dom_size, min_sum, max_sum,sum(time_taken)/float(len(time_taken)), tlp])
        print("Length of dominating set for alpha, ", alpha, " is ", len(dom_set))
        print("Total sum of weights is", sum_klouts(G, dom_set))   
        #plt.plot(sum_klouts_dom, dom_set_array, color='b', lw=1.5, zorder=1)
        #plt.scatter(sum_klouts_dom, dom_set_array,color='r', zorder=2)
        #print([min_sum-100, max_sum+100, min_dom_size-100, max_dom_size+100])
        #plt.suptitle('Algorithm B Results ')
        #plt.axis([min_sum-100, max_sum+100, min_dom_size-50, max_dom_size+50])
        #plt.xlabel("Sum of weights in dominating set")
        #plt.ylabel("Cardinality of dominating set")
        #print(os.path.splitext(sys.argv[1])[0])
        #plt.savefig(os.path.splitext(sys.argv[1])[0]+"_alpha_"+str(alpha)+"RandAlgB.png")
        
    
    
