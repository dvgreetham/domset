#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:10:55 2018

@author: dvg27
"""
""" linprog  Solve standard form LP via ADMM translated from matlab 
    from https://web.stanford.edu/~boyd/papers/admm/linprog/linprog.html
%
% [x, history] = linprog(c, A, b, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize     c'*x
%   subject to   Ax = b, x >= 0
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
"""
    
import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix
import math
import networkx as nx    
import time
import os
import csv

class hist(object):
    objval = []
    r_norm = 0
    s_norm= 0
    eps_pri=0
    eps_dual=0
    
    def __init__(self, objval, r_norm, s_norm, eps_pri, eps_dual):
       self.objval = objval
       self.r_norm  = r_norm
       self.s_norm= s_norm
       self.eps_pri= eps_pri
       self.eps_dual= eps_dual
    
def sum_weights(*args):
    if len(args)==1:
        nodes=args[0].nodes()
    elif len(args)==2:
        nodes=args[1]
    sum = 0
    for n in nodes:
        sum+=float(args[0].node[n]['weight'])

    return sum
def ADMM(c, A, b, rho, alpha):

    m, n = A.shape
  #  A_t_A = A.T.dot(A)
   # w, v = np.linalg.eig(A_t_A)
    
    MAX_ITER = 1000
    ABSTOL   = 1e-4
    RELTOL   = 1e-2

    x = np.zeros([m, 1])
    z = np.zeros([m, 1])
    u = np.zeros([m, 1])
    
    history=hist(np.zeros((n, 1)), 0, 0, 0, 0)
    
    for k in range(MAX_ITER):
        if k%50==0:
            print(k)
        cc=np.array(c).reshape(n,1)
        zz=rho*(z-u)-cc
        bb=np.array(b).reshape(n,1)

        c1=np.vstack((zz,bb))
        #print(c1.shape, c1.size)
        a1=np.concatenate((rho*np.identity(n), A.T), axis=1)
        #print(a1.shape, a1.size)
        b1=np.concatenate((A, np.zeros((m, m))), axis=1)
       # print(b1.shape, b1.size)
        d1 = np.vstack((a1, b1)) 
        #print(np.linalg.cond(d1))
        d2= csr_matrix(d1)
        tmp=lsqr(d2, c1)[0]
        tmp=tmp.reshape(2*n,1)
#        tmp=spsolve(d2,c1)
#       
        #tmp=np.linalg.lstsq(d1,c1, rcond=None)[0]
        #tmp=np.linalg.solve(d1, c1)
        x = tmp[0:n]
       
        # z-update with relaxation
        zold = z
        x_hat = alpha*x + (1 - alpha)*zold
        z=x_hat + u
        z[z<0]=0
        
    
        u = u + (x_hat - z)
    
        # diagnostics, reporting, termination checks
    
        history.objval  = cc.T.dot(x)
    
        history.r_norm  = np.linalg.norm(x - z)
        history.s_norm  = np.linalg.norm(-rho*(z - zold))
    
        history.eps_pri = math.sqrt(n)*ABSTOL + RELTOL*max(np.linalg.norm(x), np.linalg.norm(-z))
        history.eps_dual= math.sqrt(n)*ABSTOL + RELTOL*np.linalg.norm(rho*u)
    
        
       
        if (history.r_norm < history.eps_pri and history.s_norm < history.eps_dual):
             break;
        
       
    return x,history.objval 

def are_dominated(G, doms, alpha):
    viol=0
    nodes=[]
    for u in G.nodes():    
        intersect = [ii for ii in G.neighbors(u) if ii in doms]                       
        rate=len(intersect)
        if rate < math.ceil(alpha*G.degree(u)):
            viol+=math.ceil(alpha*G.degree(u))-rate
            nodes.append(u)
    return viol, nodes

def test2(alpha, file, path):
    #G=nx.read_gexf("lp3.gexf")
    G=nx.read_gexf(file)
    t=time.time()
    
    n=G.order()
    A=np.zeros((n, n))
    objfun=[]
    constraint=[]
    bounds1=[]
    index={}
    
    i=0
    for node in G.nodes():
        index[node]=i
        i=i+1
    for node, edata in G.nodes(data=True):
        objfun.append(edata['weight'])
        #print(edata)
        constraint.append(math.ceil(alpha*G.degree(node)))
        #print(math.ceil(alpha*G.degree(node)))
        bounds1.append((0,1))
        for ngh in sorted(G.neighbors(node)):
        #print index[node], index[ngh]
            A[index[node], index[ngh]]=1
  
 
    [z, objval]=ADMM(objfun, A, constraint, rho=1, alpha=1)
    #print(objval, z)
    doms=set()
    for i in range(n):
          r=np.random.uniform(0, 0.5)
          if r>=z[i]:
              z[i]=0
          else:
              doms.add(i)

    domskeys=set()
    for key, value in index.items():
        if value in doms:
            domskeys.add(key)
    [viol, nodes]=are_dominated(G, domskeys, alpha)
 
    if viol>0:
        violation=1
        l=[u for u in G.nodes(data=True) if u[0] not in domskeys]
    
        weights=sorted(l, key=lambda x: x[1]['weight'])
     
  
        while (violation):
            [viol, nds]=are_dominated(G, domskeys, alpha)
          
   
            if viol==0: 
                violation=0
            else:
                avgrate=math.ceil(viol/len(nds))
                for ii in range(avgrate):
                    domskeys.add(weights.pop()[0])    
        print("Violations", are_dominated(G, domskeys, alpha))            
    t1=time.time() - t
    print("Time taken ...", t1)
    print("Length of dominating set  for alpha ", alpha, " is", len(domskeys))
    print("Total sum of weights is", sum_weights(G, domskeys))  
    w = csv.writer(open(path, 'a'))
    w.writerow([file, alpha, len(domskeys), sum_weights(G,domskeys), t1])      
    #print(domskeys)            
if __name__ == "__main__":
   
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, "ADMM.csv")
    print(path)
    if os.path.isfile(path)==False:
        w = csv.writer(open(path, 'w'))
        w.writerow(["File", "Threshold", "DomSet-Size", "DomSet-Weights","Time taken"])
#    for i in range(10):
#        graph="pp5-100-0.18-0.0001-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:
#            test2(alpha, graph, path)   
#    for i in range(10):
#        graph="plc-500-10-0.8-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:
#            test2(alpha, graph, path)   
#    for i in range(10):
#        graph="er-500-5000-"+str(i)+".gexf"
#        for alpha in [0.25, 0.5, 0.75]:
#            test2(alpha, graph, path)  "dsjc250-5-rw.gexf","r250-1-rw.gexf", "fpsol2-i-3-rw.gexf" 
    for  graph in ["fb1-rw.gexf", "fb2-rw.gexf", "bitcoinalpha.gexf"]:
        for alpha in [0.25, 0.5, 0.75]:
            test2(alpha, graph, path)  
   