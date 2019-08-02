import networkx as nx
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import csv

def CreateGraph():
  g=nx.DiGraph()
  file=open('Wiki-Vote.txt',newline='\n')
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
 
    
     
 
def Page_Rank(g):
    page_rank=nx.pagerank(g)
    return page_rank


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

def H_Index(g):
    h_index={}
    for node in g.nodes():
       max_h_index=g.out_degree(node) 
       result=-1
       for i in range (0,max_h_index+1):
           total_nodes=0
           for neighbor in g.neighbors(node):
               if(g.out_degree(neighbor)>=i):
                   total_nodes=1+total_nodes;
           if(total_nodes>=i):
               result=max(result,i)
       h_index[node]=result      
    return h_index


def SPL(G,spreaders):
  total=0
  n=len(spreaders)
  pair=0
  for i in range (0,n):
    for j in range (i+1,n):
      if nx.has_path(G,spreaders[i],spreaders[j]):
        pair=pair+1  
        total=total+nx.shortest_path_length(G,spreaders[i],spreaders[j])
  return total

def Out_Degree_Centrality(G):
    out_degree={};
    for node in G.nodes():
        out_degree[node]=G.out_degree(node)
    return out_degree    
        
def gravity(g,mass,r ):
    grav = {}
    for node in (g.nodes()):
        grav[node] = 0
        neighbour_nodes =list(g.neighbors(node))
        for neighbour in neighbour_nodes:
            if (nx.shortest_path_length(g,source= node, target=neighbour)) > r:
                break
            if(node not in grav):
                grav[node] = 0
            if (neighbour == node):
                break
            grav[node] += (mass[neighbour] * mass[node])/((nx.shortest_path_length(g,source= node, target=neighbour))**2)  
            if((nx.shortest_path_length(g,source= node, target=neighbour))) < r:
                for n in g.neighbors(neighbour):
                    if(n not in neighbour_nodes):
                        neighbour_nodes.append(n)
        neighbour_nodes = []
        
    return grav     
    
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
    g=CreateGraph()
    G=g.copy()
    spreaders=int(g.number_of_nodes()*fraction)
    betath=beta_threshold(g)
    
   # CASE1: IMPROVED CORENESS * EIGEN VECTOR CENTRALITY
    
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
    print("ICC done")    
    
  
    #Case 3: page rank:
    page_rank_centrality=Page_Rank(G)
    spreaders_list_page_rank=Find_Spreaders(G,page_rank_centrality,spreaders)
    infected_scale3,ftc3=SIRModel(G,spreaders_list_page_rank,beta,simulations)
    for i in range(0,50):
        infected_scale3[i]=infected_scale3[i]/simulations
    ftc3=ftc3/simulations    
    time3=[]
    for i in range (0,50):
        time3.append(5*i)    
    print("Page Rank done")    
      
    #case 4:extended* h index    
    h_index=H_Index(G) 
    hybrid_centrality=Hybrid_Centrality(G,extended_coreness,h_index)   
    spreaders_list_extended_coreness_h_index=FindSpreaders(G,hybrid_centrality,spreaders)
    infected_scale4,ftc4=SIRModel(G,spreaders_list_extended_coreness_h_index,beta,simulations)
    for i in range(0,50):
        infected_scale4[i]=infected_scale4[i]/simulations
    ftc4=ftc4/simulations    
    time4=[]
    for i in range (0,50):
        time4.append(5*i)   
    print('improved hybrid done')  
        
   #case 5: k shell
    spreaders_list_k_shell=Find_Spreaders(G,k_shell,spreaders)
    infected_scale5,ftc5=SIRModel(G,spreaders_list_k_shell,beta,simulations)
    for i in range(0,50):
        infected_scale5[i]=infected_scale5[i]/simulations
    ftc5=ftc5/simulations
    time5=[]
    for i in range (0,50):
        time5.append(5*i)    
    print('k shell done')  
        
    #case 6: extended coreness
    spreaders_list_extended_coreness=FindSpreaders(G,extended_coreness,spreaders)
    infected_scale6,ftc6=SIRModel(G,spreaders_list_extended_coreness,beta,simulations)
    for i in range(0,50):
        infected_scale6[i]=infected_scale6[i]/simulations
    time6=[]
    ftc6=ftc6/simulations
    for i in range (0,50):
        time6.append(5*i)
    
    print('extended done')   
      
    
   #case 8 : summing
    eigen=Eigen_Vector_Centrality(g)
    degree= Out_Degree_Centrality(g)
    page=nx.pagerank(g)
    n_degree=dict()
    mean=np.array(list(degree.values())).mean()
    sd=np.array(list(degree.values())).std()
    for node in g.nodes():
        n_degree[node]=(degree[node]-mean)/sd
    
    n_eigen=dict()
    mean=np.array(list(eigen.values())).mean()
    sd=np.array(list(eigen.values())).std()
    for node in g.nodes():
        n_eigen[node]=(eigen[node]-mean)/sd
    
    n_page=dict()
    mean=np.array(list(page.values())).mean()
    sd=np.array(list(page.values())).std()
    for node in g.nodes():
        n_page[node]=(page[node]-mean)/sd
    
    master=dict()
    for node in g.nodes():
        master[node]=n_degree[node]+n_eigen[node]+n_page[node]       
    spreaders_list_sum=Find_Spreaders(G,master,spreaders)
    infected_scale8,ftc8=SIRModel(G,spreaders_list_sum,beta,simulations)
    for i in range(0,50):
        infected_scale8[i]=infected_scale8[i]/simulations
    time8=[]
    ftc8=ftc8/simulations
    for i in range (0,50):
        time8.append(5*i)
    
    print('sum done')  
    
    
   #gravity
    k_shell = nx.core_number(g) 
    grav_k_shell = gravity(g,k_shell,r=2)
    spreaders_list_grav_k_shell=Find_Spreaders(g,grav_k_shell,spreaders)
    infected_scale9,ftc9=SIRModel(g, spreaders_list_grav_k_shell,beta,simulations)
    for i in range(0,50):
        infected_scale9[i]=infected_scale9[i]/simulations
    ftc9=ftc9/simulations    
    time9=[]
    for i in range (0,50):                                                                                  
        time9.append(5*i) 
    print("grav_k_shell done")
    
    
    xnew1 = np.linspace(min(time1[:-30]),max(time1[:-30]),50)
    power_smooth1 = spline(time1[:-30],infected_scale1[:-30],xnew1)
        
     

    xnew3 = np.linspace(min(time3[:-30]),max(time3[:-30]),50)
    power_smooth3 = spline(time3[:-30],infected_scale3[:-30],xnew3)

    xnew4 = np.linspace(min(time4[:-30]),max(time4[:-30]),50)
    power_smooth4 = spline(time4[:-30],infected_scale4[:-30],xnew4)
    
    
    xnew5 = np.linspace(min(time5[:-30]),max(time5[:-30]),50)
    power_smooth5 = spline(time5[:-30],infected_scale5[:-30],xnew5)
     
     
    xnew7 = np.linspace(min(time7[:-30]),max(time7[:-30]),50)
    power_smooth7 = spline(time7[:-30],infected_scale7[:-30],xnew7)
    
    xnew8 = np.linspace(min(time8[:-30]),max(time8[:-30]),50)
    power_smooth8 = spline(time8[:-30],infected_scale8[:-30],xnew8)
    
    
    xnew9 = np.linspace(min(time9[:-30]),max(time9[:-30]),50)
    power_smooth9 = spline(time9[:-30],infected_scale9[:-30],xnew9)
    
    
    plt.plot(xnew4[:-35],power_smooth4[:-35],label='IH',marker='*')
    plt.plot(xnew1[:-35],power_smooth1[:-35],label='Hybrid',marker='v')  
  #  plt.plot(xnew7[:-35],power_smooth7[:-35],label='Page+Cluster+Degree+Eigen',marker='v')  
    plt.plot(xnew3[:-35],power_smooth3[:-35],label='PR',marker='X')
    plt.plot(xnew6[:-35],power_smooth6[:-35],label='Extended',marker='X')
    plt.plot(xnew5[:-35],power_smooth5[:-35],label='K-shell',marker='o')
##    plt.plot(xnew2[:-35],power_smooth2[:-35],label='Eigen',marker='.')
    plt.plot(xnew8[:-35],power_smooth8[:-35],label='SSD',marker='.')
    plt.plot(xnew9[:-35],power_smooth9[:-35],label='Gravity',marker='.')
    
   
    
    plt.ylabel("Infected Scale F($_t$)")
    plt.xlabel("Time")
    plt.legend()
    