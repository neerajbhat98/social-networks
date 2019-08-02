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
   spreaders= NonNeighboringSpreaderFilter(g,master_centrality,no_of_spreaders)
   #spreaders= AntiNonNeighboringSpreaderFilter(g,master_centrality,no_of_spreaders,all_list)
  
   return spreaders

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
 
#implementing SIR model
def SIRModel(g,spreaders,beta,simulations):
 infected_scale=np.array(np.zeros(50))
 ftc=0
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
   infected_scale[time_stamp]=infected_scale[time_stamp]+(infected_nodes+recovered_nodes)/n
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
     infected_scale[time_stamp]=infected_scale[time_stamp]+(infected_nodes+recovered_nodes)/n
     infected=susceptible_to_infected  
  
 return infected_scale,ftc 

if __name__=="__main__":
    
    fraction=0.06
    simulations=100 
    beta=0.06
    spreaders=int(g.number_of_nodes()*fraction)
    row=['Centrality','Beta','fraction_of_Spreaders','Ftc']
    with open('CondMatresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
   # CASE1: IMPROVED CORENESS * EIGEN VECTOR CENTRALITY
    g=CreateGraph()
    G=g.copy()
    betath=beta_threshold(g)
    k_shell=nx.core_number(g)
    coreness=Coreness(G,k_shell)
    extended_coreness=Extended_Coreness(G,coreness)
    eigen_vector_centrality=Eigen_Vector_Centrality(G)
    hybrid_centrality=Hybrid_Centrality(G,coreness,eigen_vector_centrality)   
    spreaders_list_coreness_eigen=FindSpreaders(G,hybrid_centrality,spreaders)
    infected_scale1,ftc1=SIRModel(G, spreaders_list_coreness_eigen,beta,simulations)
    for i in range(0,50):
        infected_scale1[i]=infected_scale1[i]/simulations
    ftc1=ftc1/simulations    
    time1=[]
    for i in range (0,50):                                                                                  
        time1.append(5*i)
  
    row=["ICC",beta,fraction,ftc1]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)    
    print("ICC done")
    
   #CASE2: EIGEN VECTOR CENTRALITY
    eigen_vector_centrality=Eigen_Vector_Centrality(G)
    spreaders_list_eigen=FindSpreaders(G,eigen_vector_centrality,spreaders)
    infected_scale2,ftc2=SIRModel(G,spreaders_list_eigen,beta,simulations)
    for i in range(0,50):
        infected_scale2[i]=infected_scale2[i]/simulations
    time2=[]
    ftc2=ftc2/simulations
    for i in range (0,50):
        time2.append(5*i)   
    row=["Eigen",beta,fraction,ftc2]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print("Eigen done")
    
    #CASE3: DEGREE CENTRALITY
    degree_centrality=Degree_Centrality(G)
    spreaders_list_degree=FindSpreaders(G,degree_centrality,spreaders)
    infected_scale3,ftc3=SIRModel(G,spreaders_list_degree,beta,simulations)
    for i in range(0,50):
        infected_scale3[i]=infected_scale3[i]/simulations
    time3=[]
    ftc3=ftc3/simulations
    for i in range (0,50):
        time3.append(5*i)   
    row=["Degree",beta,fraction,ftc3]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print("degree done")  
    
    #CASE4 K-SHELL CENTRALITY 
    spreaders_list_k_shell=FindSpreaders(G,k_shell,spreaders)
    infected_scale4,ftc4=SIRModel(G,spreaders_list_k_shell,beta,simulations)
    for i in range(0,50):
        infected_scale4[i]=infected_scale4[i]/simulations
    time4=[]
    ftc4=ftc4/simulations
    for i in range (0,50):
        time4.append(5*i) 
    row=["K_Shell",beta,fraction,ftc4]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print("k-shell done")
        
              
   # CASE 5: EXTENDED CORENESS * H-INDEX   
    h_index=H_Index(G) 
    hybrid_centrality=Hybrid_Centrality(G,extended_coreness,h_index)   
    spreaders_list_extended_coreness_h_index=FindSpreaders(G,hybrid_centrality,spreaders)
    infected_scale5,ftc5=SIRModel(G,spreaders_list_extended_coreness_h_index,beta,simulations)
    for i in range(0,50):
        infected_scale5[i]=infected_scale5[i]/simulations
    time5=[]
    ftc5=ftc5/simulations
    for i in range (0,50):
        time5.append(5*i)  
    row=["Improved Hybrid",beta,fraction,ftc5]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print("improved hybrid done")
        
    #CASE 6: EXTENDED_CORENESS
    spreaders_list_extended_coreness=FindSpreaders(G,extended_coreness,spreaders)
    infected_scale6,ftc6=SIRModel(G,spreaders_list_extended_coreness,beta,simulations)
    for i in range(0,50):
        infected_scale6[i]=infected_scale6[i]/simulations
    time6=[]
    ftc6=ftc6/simulations
    for i in range (0,50):
        time6.append(5*i)
    row=["Extended",beta,fraction,ftc6]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print('extended coreness done')  
  
    #CASE 7 : ARKA SIR
    eigen=eigen_vector_centrality
    degree=degree_centrality
    page=nx.pagerank(G)
    cluster=nx.clustering(G)
    centralities2=dict()
    centralities2['page_rank']=page
    centralities2['degree']=degree
    centralities2['eig']=eigen
    centralities2['cluster_rank']=cluster
    spreaders_list=Find_Spreaders_arka_sir(g,centralities2,spreaders)
    infected_scale7,ftc7=SIRModel(G,spreaders_list,beta,simulations)
    for i in range(0,50):
        infected_scale7[i]=infected_scale7[i]/simulations
    time7=[]
    ftc7=ftc7/simulations
    for i in range (0,50):
        time7.append(5*i)
    row=["Arka sir base",beta,fraction,ftc7]
    with open('CondMAtresults.csv','a') as f:
      w=csv.writer(f)
      w.writerow(row)
    print('arka sir done')  
  
 
    
    xnew1 = np.linspace(min(time1[:-30]),max(time1[:-30]),50)
    power_smooth1 = spline(time1[:-30],infected_scale1[:-30],xnew1)
    
    xnew2 = np.linspace(min(time2[:-30]),max(time2[:-30]),50)
    power_smooth2 = spline(time2[:-30],infected_scale2[:-30],xnew2)
     
    xnew3 = np.linspace(min(time3[:-30]),max(time3[:-30]),50)
    power_smooth3 = spline(time3[:-30],infected_scale3[:-30],xnew3)
     
    xnew4 = np.linspace(min(time4[:-30]),max(time4[:-30]),50)
    power_smooth4 = spline(time4[:-30],infected_scale4[:-30],xnew4)
    
    xnew5 = np.linspace(min(time5[:-30]),max(time5[:-30]),50)
    power_smooth5 = spline(time5[:-30],infected_scale5[:-30],xnew5)

    xnew6 = np.linspace(min(time6[:-30]),max(time6[:-30]),50)
    power_smooth6 = spline(time6[:-30],infected_scale6[:-30],xnew6)    

    
    xnew7 = np.linspace(min(time7[:-30]),max(time7[:-30]),50)
    power_smooth7 = spline(time7[:-30],infected_scale7[:-30],xnew7)       
    fig=plt.figure() 
    
    plt.plot(xnew5[:-35],power_smooth5[:-35],label='Improved Hybrid',marker='*')
    plt.plot(xnew1[:-35],power_smooth1[:-35],label='Hybrid',marker='v')  
    plt.plot(xnew7[:-35],power_smooth7[:-35],label='page+cluster+degree+eigen',marker='v')  
    plt.plot(xnew6[:-35],power_smooth6[:-35],label='Extended Coreness',marker='.')
     
    plt.plot(xnew4[:-35],power_smooth4[:-35],label='K-shell',marker='X')
    plt.plot(xnew2[:-35],power_smooth2[:-35],label='Eigen',marker='o')  
     
   # plt.plot(xnew3[:-35],power_smooth3[:-35],label='degree')
    
    
    plt.ylabel("Infected Scale")
    plt.xlabel("Time")
    plt.legend()
    
    
   
    fig.savefig("Condmat")  ###please provide the file name
    