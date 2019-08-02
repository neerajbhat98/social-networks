import networkx as nx
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import csv

def CreateGraph():
  g=nx.Graph()
  file=open('CA-CondMat.txt',newline='\n')
  #creating a graph from the relationship given in the condmat file
  for word in file:
    data=word.split()
    u=data[0]
    v=data[1]
    g.add_edge(u,v)    
  #total no of nodes and edges formed in the network
  g.remove_edges_from(g.selfloop_edges())
  print("Number of Nodes =" ,g.number_of_nodes())  
  print("Number of Edges =" ,g.number_of_edges())    
  return g

def beta_threshold(g):
    k=0
    n=g.number_of_nodes()
    k_square=0
    for node in g.nodes():
        k=k+g.degree(node)
        k_square=k_square+pow(g.degree(node),2)
    return k/k_square

def Coreness(g,k_shell):
  coreness={}
  for node in g.nodes():
    sum=0
    for neighbor in g.neighbors(node):
        sum=sum + k_shell[neighbor]   
    #sum=sum+d[node]
    coreness[node]=sum
 
  return coreness
 
def beta_threshold(g):
    k=0
    n=g.number_of_nodes()
    k_square=0
    for node in g.nodes():
        k=k+g.degree(node)
        k_square=k_square+pow(g.degree(node),2)
    return k/k_square    

def Extended_Coreness(g,coreness):
  extended_coreness={}
  for node in g.nodes():
    sum=0
    for neighbor in g.neighbors(node):
        sum=sum+coreness[neighbor]
    extended_coreness[node]=sum
  
  return extended_coreness
        

def Eigen_Vector_Centrality(g):  
  eigen_vector_centrality=nx.eigenvector_centrality(g)
  return eigen_vector_centrality


#calculating the hybrid centrality of each and every node in the network
def Hybrid_Centrality(g,Centrality1,Centrality2):
  hybrid_centrality={}
  for node in g.nodes():
      hybrid_centrality[node]=Centrality1[node]*Centrality2[node]
  return hybrid_centrality
 
    
def Degree_Centrality(g):
    degree_centrality=nx.degree_centrality(G)
    return degree_centrality  

    
 
def Page_Rank(g):
    page_rank=nx.pagerank(g)
    return page_rank

#finding out spreaders
def FindSpreaders(g,centrality,no_of_spreaders):
  spreaders_list=centrality   
  spreaders=[]    
  while no_of_spreaders>0:
    no_of_spreaders=no_of_spreaders-1
    maxvalue=-1;
    maxkey=0
    for k,v in spreaders_list.items():
        if v >=maxvalue:
            maxvalue=v
            maxkey=k
    spreaders.append(maxkey)
    for neighbor in g.neighbors(maxkey): 
          spreaders_list[neighbor]=-1
    spreaders_list[maxkey]=-1
  return spreaders


def H_Index(g):
    h_index={}
    for node in g.nodes():
       max_h_index=g.degree(node) 
       result=-1
       for i in range (0,max_h_index+1):
           total_nodes=0
           for neighbor in g.neighbors(node):
               if(g.degree(neighbor)>=i):
                   total_nodes=1+total_nodes;
           if(total_nodes>=i):
               result=max(result,i)
       h_index[node]=result      
    return h_index


def Find_Spreaders(g,centrality,no_of_spreaders):
   spreaders_list=centrality   
   spreaders=[]    
   while no_of_spreaders>0:
      no_of_spreaders=no_of_spreaders-1 
      matches = sorted(spreaders_list.items(), key=lambda kv: kv[1], reverse=True)
      key=(matches[0][0])
      #print(key)
      spreaders.append(matches[0][0])
      del spreaders_list[key]    
   return spreaders

def NonNeighboringSpreaderFilter(g,centrality,no_of_spreaders):
  spreaders_list=centrality   
  spreaders=[]    
  while no_of_spreaders>0:
    no_of_spreaders=no_of_spreaders-1
    maxvalue=-1;
    maxkey=0
    for k,v in spreaders_list.items():
        if v >=maxvalue:
            maxvalue=v
            maxkey=k
    spreaders.append(maxkey)
    for neighbor in g.neighbors(maxkey): 
          spreaders_list[neighbor]=-1
    spreaders_list[maxkey]=-1
  return spreaders

def Find_Spreaders_arka_sir(g,centralities,no_of_spreaders):
   master_centrality=dict()
   avg_centrality=dict()
   sd_centrality=dict()
   all_list=dict()
   count=0
   for name,current_centrality in centralities.items():
       avg_centrality[name]=np.array(list(current_centrality.values())).mean()
       sd_centrality[name]=np.array(list(current_centrality.values())).std()
   spreaders=list()
   flag=0
   for node in g.nodes():
       flag=0
       x=0
       for name,current_centrality in centralities.items():
           x+=(current_centrality[node]-avg_centrality[name])/sd_centrality[name]
           if current_centrality[node]<(avg_centrality[name]):
              flag=1
       if flag==0:
           count+=1
           master_centrality.update({node:x})
           #master_centrality.update({node:centralities['page_rank'][node]})
       all_list.update({node:x})
   #spreaders=sorted(master_centrality, key=lambda key: master_centrality[key], reverse=True)[0:no_of_spreaders]
  # spreaders= NonNeighboringSpreaderFilter(g,master_centrality,no_of_spreaders)
   #spreaders= AntiNonNeighboringSpreaderFilter(g,master_centrality,no_of_spreaders,all_list)
  
   return master_centrality

 

def Cluster_Coefficient(g):
    cluster_rank=nx.clustering(g)
    return cluster_rank
 
def SPL(G,spreaders):
  total=0
  n=len(spreaders)
  for i in range (0,n):
    for j in range (i+1,n):
      if nx.has_path(G,spreaders[i],spreaders[j]):
        total=total+nx.shortest_path_length(G,spreaders[i],spreaders[j])
  total=total/(n*n-n)  
  total=2*total
  return total

def SIR_Rank(g,beta,simulations,centrality):
    rank_dict=dict()
    for node in g.nodes():
        spreaders=list()
        spreaders.append(node)
        ftc=SIRModel(g,spreaders,beta,simulations)
        rank_dict.update({node:ftc})
    ranks=sorted(rank_dict, key=lambda key: rank_dict[key], reverse=True)
    return ranks

def Spreaders_Rank(centrality):
     ranks=sorted(centrality, key=lambda key: centrality[key], reverse=True)  
     return ranks

def kendall_tau(spreaders_rank,sir_rank):
    spreaders_ranking=dict()
    sirs_ranking=dict()
    x1=list()
    x2=list()
    index=0
    #create hashmap of rankings in spreaders
    for node in spreaders_rank:
        spreaders_ranking.update({node:index})
        index+=1
    index=0
    #create hashmap of rankings in sir
    for node in sir_rank:
        sirs_ranking.update({node:index})
        index+=1
    for node in g.nodes():
        x1.append(spreaders_ranking[node])
        x2.append(sirs_ranking[node])    
    tau, p_value = stats.kendalltau(x1, x2)
    return tau

#implementing SIR model
def SIRModel(g,spreaders,beta,simulations):
 ftc=0
 s=simulations
 while simulations>0 :
   simulations=simulations-1
   #print(simulations)
   infected=spreaders
  
   status={}
   for i in g.nodes():
      status[i]=0;
   for i in infected:
     status[i]=1;
  
   n=g.number_of_nodes() 
   infected_nodes=len(infected)
   recovered_nodes=0
   time_stamp=0
   infected=spreaders
  # print(infected)
   while(len(infected)>0):
     susceptible_to_infected=[]
     time_stamp=time_stamp+1
     #print("time=",time_stamp)
     #print("infected=",infected)
     for i in infected:
        susceptible=[]
        status[i]=2
        for neighbor in g.neighbors(i):
            if(status[neighbor]==0):
                susceptible.append(neighbor)
        #print("susceptible=",susceptible)        
        total_susceptible=len(susceptible)
        #print("total no of susceptible nodes are=",total_susceptible)
        no_of_susceptible_to_infected=round(beta*total_susceptible)
        #print('after calculating probability=', no_of_susceptible_to_infected)
        while no_of_susceptible_to_infected>0:
            random_index=random.randint(0,total_susceptible-1)
            if susceptible[random_index] not in susceptible_to_infected:
             susceptible_to_infected.append(susceptible[random_index])
             status[susceptible[random_index]]=1
             no_of_susceptible_to_infected=no_of_susceptible_to_infected-1
             #print("infected to be =",susceptible[random_index])
     infected_nodes=len(susceptible_to_infected)
     recovered_nodes=len(infected)
    # print("infected :",infected_nodes)
    # print("recovered:",recovered_nodes)
     ftc=ftc+(recovered_nodes)/n 
     infected=susceptible_to_infected  
  
 return ftc /s

if __name__=="__main__":
    
    
    fraction=0.06
    simulations=5 
    
   # CASE1: IMPROVED CORENESS * EIGEN VECTOR CENTRALITY
    g=CreateGraph() 
    spreaders=int(g.number_of_nodes()*fraction)
    G=g.copy()
    betath=beta_threshold(g)
    beta=betath
    k_shell=nx.core_number(g)
    coreness=Coreness(G,k_shell)
    extended_coreness=Extended_Coreness(G,coreness)
    eigen_vector_centrality=Eigen_Vector_Centrality(G)
    hybrid_centrality=Hybrid_Centrality(G,coreness,eigen_vector_centrality)   
    sir_rank=SIR_Rank(g,beta,simulations,hybrid_centrality)
    spreaders_rank=Spreaders_Rank(hybrid_centrality)
    kd1=kendall_tau(spreaders_rank,sir_rank)
    print(kd1)
   #CASE2: EIGEN VECTOR CENTRALITY
    eigen_vector_centrality=Eigen_Vector_Centrality(G)
    sir_rank=SIR_Rank(g,beta,simulations,eigen_vector_centrality)
    spreaders_rank=Spreaders_Rank(eigen_vector_centrality)
    kd2=kendall_tau(spreaders_rank,sir_rank)
    print(kd2)
    
    #CASE3: DEGREE CENTRALITY
    degree_centrality=Degree_Centrality(G)
    spreaders_list_degree=FindSpreaders(G,degree_centrality,spreaders)
    
    #CASE4 K-SHELL CENTRALITY 
    sir_rank=SIR_Rank(g,beta,simulations,k_shell)
    spreaders_rank=Spreaders_Rank(k_shell)
    kd4=kendall_tau(spreaders_rank,sir_rank)
    print(kd4)    
              
   # CASE 5: EXTENDED CORENESS * H-INDEX   
    h_index=H_Index(G) 
    hybrid_centrality=Hybrid_Centrality(G,extended_coreness,h_index)   
    sir_rank=SIR_Rank(g,beta,simulations,hybrid_centrality)
    spreaders_rank=Spreaders_Rank(hybrid_centrality)
    kd5=kendall_tau(spreaders_rank,sir_rank)
    print(kd5)
 
        
    #CASE 6: EXTENDED_CORENESS
    sir_rank=SIR_Rank(g,beta,simulations,extended_coreness)
    spreaders_rank=Spreaders_Rank(extended_coreness)
    kd6=kendall_tau(spreaders_rank,sir_rank)
    print(kd6)
    
    
    #case 7 : arka sir wala
    eigen=eigen_vector_centrality
    degree=degree_centrality
    page=nx.pagerank(G)
    cluster=nx.clustering(G)
    centralities2=dict()
    centralities2['page_rank']=page
    centralities2['degree']=degree
    centralities2['eig']=eigen
    centralities2['cluster_rank']=cluster
    master_centrality=Find_Spreaders_arka_sir(g,centralities2,spreaders)
    sir_rank=SIR_Rank(g,beta,simulations,master_centrality)
    spreaders_rank=Spreaders_Rank(master_centrality)
    kd7=kendall_tau(spreaders_rank,sir_rank)
    print(kd7)
 
 