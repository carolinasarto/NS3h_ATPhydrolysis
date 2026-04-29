"""
Conserved Water Site Analysis Script
------------------------------------

This script analyzes conserved water binding sites across dengue virus serotypes
(DEN1–DEN4) in both apo and RNA-bound systems.

Workflow:
- Loads water cluster data (watcent.pdb, watclust.dat) for each serotype and atom type
  (ETA_C1_DATA, ETA_O1_DATA).
- Removes duplicate water sites closer than 2.5 Å.
- Identifies conserved sites across serotypes if they are within 4.5 Å.
- Compares multiple serotype pairings (D1–D3 vs D2–D4, D1–D2 vs D3–D4, D1–D4 vs D2–D3)
  to detect sites conserved across all four.
- Calculates midpoints for conserved sites.
- Outputs counts, coordinates, and new PDB files (prot_<atom>_<system>.pdb) with
  conserved waters appended to reference protein structures.

Usage:
1. Ensure directory structure includes den1/, den2/, den3/, den4/ with apo/ and RNA/
   subfolders, each containing all_aligned/ETA_C1_DATA/ and ETA_O1_DATA/ with
   watcent.pdb and watclust.dat.
2. Provide reference protein PDB files (ref_apo.pdb, ref_RNA.pdb).
3. Run the script:
   python conserved_sites.py

"""


import numpy as np
import pandas as pd
import os
import subprocess
import shutil
import sys
import math
import itertools
import matplotlib.pyplot as plt


analisis='all_aligned/'
#cutoff_dist=4.0
cutoff_repetidos=2.5
#cutoff_repetidos_rna=4.5
#cutoff_repetidos_apo=3.0
cutoff_dist_lista=[4.5] #np.arange(1.0,3.5,0.5)

def distancia (seriea,a,serieb,b):
  dx=seriea['x'].iloc[a]-serieb['x'].iloc[b]
  dy=seriea['y'].iloc[a]-serieb['y'].iloc[b]
  dz=seriea['z'].iloc[a]-serieb['z'].iloc[b]
  d=math.sqrt(dx*dx+dy*dy+dz*dz)
  return d

def punto_medio (seriea,a,serieb,b):
  equis=(seriea['x'].iloc[a]+serieb['x'].iloc[b])/2
  ye=(seriea['y'].iloc[a]+serieb['y'].iloc[b])/2
  zeta=(seriea['z'].iloc[a]+serieb['z'].iloc[b])/2
  pm=[equis,ye,zeta]
  return pm

def leer_linea (un_archivo, linea):
    with open(un_archivo, 'r') as file:
        for current_line_number, line in enumerate(file, start=0):
            if current_line_number == linea:
                return  line
                break                

dir_base='/home/usuario/ns3/inhibidores/md_mixsv/'

for cutoff_dist in cutoff_dist_lista:
    print(cutoff_dist)
    conservados=[[],[],[],[]]

    z=0
    for system in ['apo/','RNA/']:

        ## Hago busqueda separada de sitios de C y de O
        for atom in ['ETA_C1_DATA','ETA_O1_DATA']:
            
##             print(system+' '+atom+'----------1-----------Primero chequeo repetidos y guardos los datos de cada sistema')
  
             dfs=[]
             dfs_watcent=[]
  
             ## Hay 4 serotipos de dengue den1, den2, den3, den4
             for i in range (1,5):
                 
                 watcent_data = dir_base+"den%d/" % i+system + analisis + atom+'/watcent.pdb'
                 print(watcent_data)
                 watclust_data = dir_base+"den%d/" % i+system + analisis + atom+'/watclust.dat'
                 
                 df_wat_cent = pd.read_csv(watcent_data, delimiter='\s+')
  
                 df_info = pd.read_csv(watclust_data, delimiter='\s+')
                 df_coord = pd.read_csv(watcent_data, delimiter='\s+',usecols=[2,3,5,6,7,8],names=['atom','solvent','number','x','y','z'])
                 
             ## Eliminar los sitios repetidos
  
                 todos_distintos=True
                 c=0
                 loc_repetidos=[]
                 while todos_distintos and c < len(df_coord):
                     for k in range (len(df_coord)-1):
                         for j in range (k+1, len(df_coord)):
  
                             dist_sites = distancia(df_coord,k,df_coord,j)
                             if dist_sites < cutoff_repetidos:
                                 print(system+" Den"+str(i)+" tiene sitios a menos de "+str(cutoff_repetidos)+' en '+atom+": "+str(k+1)+" y "+str(j+1)) 
                                 loc_repetidos.append(k) 
                                 todos_distintos=False
                                 
                             else:
                                 c+=1
  
                 if todos_distintos == False:
                     for dataframe in [df_wat_cent, df_info, df_coord]:
                         dataframe = dataframe.drop(index=loc_repetidos, inplace=True)
                 else:
                      print("") #Son todos sitios distintos para Den"+str(i)) 
  
                 print(system+' '+atom+' D'+str(i)+' '+str(len(df_wat_cent)))
                 
                 dfs_watcent.append(df_wat_cent)
                 df_concat = pd.concat([df_info,df_coord],axis=1)
                 dfs.append(df_concat)
                 
     #############################################################################################################################      
             print('--------------2---------------------Comparo D1 y D3 primero D2 y D4 despues')
             ## Busqueda de sitios conservados entre den1 y den3
  
             sites_d1_d3=[]
  
             for i in range(len(dfs[0])):
                 for j in range(len(dfs[2])):
                     dist_sites = distancia(dfs[0],i,dfs[2],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[0],i,dfs[2],j)
  
                         info_sites=[]
                         info_sites.append(dfs[0]['number'].iloc[i])
                         info_sites.append(dfs[2]['number'].iloc[j])
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
                 
                         sites_d1_d3.append(info_sites)
             df_sites_d1_d3 = pd.DataFrame(sites_d1_d3, columns=['Site d1', 'Site d3', 'Distance', 'x', 'y', 'z'])
             print('D1-D3: '+str(len(df_sites_d1_d3)))
  
             ## Busqueda de sitios conservados entre den1 y den2
             
             sites_d4_d2=[]
  
             for i in range(len(dfs[3])):
                 for j in range(len(dfs[1])):
                     dist_sites = distancia(dfs[3],i,dfs[1],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[3],i,dfs[1],j)
  
                         info_sites=[]
                         info_sites.append(dfs[3]['number'].iloc[i])
                         info_sites.append(dfs[1]['number'].iloc[j])
                         info_sites.append(dist_sites)
  
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d4_d2.append(info_sites)
  
             df_sites_d4_d2 = pd.DataFrame(sites_d4_d2, columns=['Site d4', 'Site d2', 'Distance', 'x', 'y', 'z'])
             
             print('D4-D2: '+str(len(df_sites_d4_d2)))
  
             ## Busqueda de sitios en los 4 serotipos
  
             sites_d1_d2_d3_d4=[]
  
             for i in range(len(sites_d4_d2)):
                 for j in range(len(sites_d1_d3)):
                     dist_sites = distancia(df_sites_d4_d2,i,df_sites_d1_d3,j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (df_sites_d4_d2,i,df_sites_d1_d3,j)
  
                         info_sites=[]
                         info_sites.append(df_sites_d4_d2['Site d4'].iloc[i])
                         info_sites.append(df_sites_d4_d2['Site d2'].iloc[i])
                         info_sites.append(df_sites_d1_d3['Site d1'].iloc[j])
                         info_sites.append(df_sites_d1_d3['Site d3'].iloc[j])
                         
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d1_d2_d3_d4.append(info_sites)
  
  
             df_sites_d1_d2_d3_d4 = pd.DataFrame(sites_d1_d2_d3_d4, columns=['Site d4', 'Site d2', 'Site d1', 'Site d3', 'Distance', 'x', 'y', 'z'])
             
            ## Eliminar los sitios repetidos

             todos_distintos=True
             c=0
             loc_repetidos=[]
             while todos_distintos and c < len(df_sites_d1_d2_d3_d4):
                 for k in range (len(df_sites_d1_d2_d3_d4)-1):
                     for j in range (k+1, len(df_sites_d1_d2_d3_d4)):

                         dist_sites = distancia(df_sites_d1_d2_d3_d4,k,df_sites_d1_d2_d3_d4,j)
                         if dist_sites < cutoff_repetidos:
#                             print(system+" Den"+str(i)+" tiene sitios a menos de "+str(cutoff_repetidos)+' en '+atom+": "+str(k+1)+" y "+str(j+1))
                             loc_repetidos.append(k)
                             todos_distintos=False

                         else:
                             c+=1

             if todos_distintos == False:
                 for dataframe in [df_sites_d1_d2_d3_d4]:
                     dataframe = dataframe.drop(index=loc_repetidos, inplace=True)
             else:
                  print("Son todos sitios distintos") 

             conservados[z].append(len(df_sites_d1_d2_d3_d4))
             print('D1-D2-D3-D4: '+str(len(df_sites_d1_d2_d3_d4)))
             print(df_sites_d1_d2_d3_d4)
  
     ###############################################################################################33
             print('--------------3---------------------Comparo D1 y D2 primero D3 y D4 despues')
             ## Busqueda de sitios conservados entre den1 y den2
  
             sites_d1_d2=[]
  
             for i in range(len(dfs[0])):
                 for j in range(len(dfs[1])):
                     dist_sites = distancia(dfs[0],i,dfs[1],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[0],i,dfs[1],j)
  
                         info_sites=[]
                         info_sites.append(dfs[0]['number'].iloc[i])
                         info_sites.append(dfs[1]['number'].iloc[j])
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d1_d2.append(info_sites)
  
             df_sites_d1_d2 = pd.DataFrame(sites_d1_d2, columns=['Site d1', 'Site d2', 'Distance', 'x', 'y', 'z'])
             print('D1-D2: '+str(len(df_sites_d1_d2)))
  
             ## Busqueda de sitios conservados entre den3 y den4
  
             sites_d3_d4=[]
  
             for i in range(len(dfs[2])):
                 for j in range(len(dfs[3])):
                     dist_sites = distancia(dfs[2],i,dfs[3],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[2],i,dfs[3],j)
  
                         info_sites=[]
                         info_sites.append(dfs[2]['number'].iloc[i])
                         info_sites.append(dfs[3]['number'].iloc[j])
                         info_sites.append(dist_sites)
  
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d3_d4.append(info_sites)
  
             df_sites_d3_d4 = pd.DataFrame(sites_d3_d4, columns=['Site d3', 'Site d4', 'Distance', 'x', 'y', 'z'])
  
             print('D3-D4: '+str(len(df_sites_d3_d4)))
  
             ## Busqueda de sitios en los 4 serotipos
  
             sites_d1_d2_d3_d4=[]
  
             for i in range(len(sites_d1_d2)):
                 for j in range(len(sites_d3_d4)):
                     dist_sites = distancia(df_sites_d1_d2,i,df_sites_d3_d4,j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (df_sites_d1_d2,i,df_sites_d3_d4,j)
  
                         info_sites=[]
                         info_sites.append(df_sites_d1_d2['Site d1'].iloc[i])
                         info_sites.append(df_sites_d1_d2['Site d2'].iloc[i])
                         info_sites.append(df_sites_d3_d4['Site d3'].iloc[j])
                         info_sites.append(df_sites_d3_d4['Site d4'].iloc[j])
  
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d1_d2_d3_d4.append(info_sites)
  
  
             df_sites_d1_d2_d3_d4 = pd.DataFrame(sites_d1_d2_d3_d4, columns=['Site d1', 'Site d2', 'Site d3', 'Site d4', 'Distance', 'x', 'y', 'z'])
  
            ## Eliminar los sitios repetidos

             todos_distintos=True
             c=0
             loc_repetidos=[]
             while todos_distintos and c < len(df_sites_d1_d2_d3_d4):
                 for k in range (len(df_sites_d1_d2_d3_d4)-1):
                     for j in range (k+1, len(df_sites_d1_d2_d3_d4)):

                         dist_sites = distancia(df_sites_d1_d2_d3_d4,k,df_sites_d1_d2_d3_d4,j)
                         if dist_sites < cutoff_repetidos:
#                             print(system+" Den"+str(i)+" tiene sitios a menos de "+str(cutoff_repetidos)+' en '+atom+": "+str(k+1)+" y "+str(j+1))
                             loc_repetidos.append(k)
                             todos_distintos=False

                         else:
                             c+=1

             if todos_distintos == False:
                 for dataframe in [df_sites_d1_d2_d3_d4]:
                     dataframe = dataframe.drop(index=loc_repetidos, inplace=True)
             else:
                  print("Son todos sitios distintos")

             conservados[z].append(len(df_sites_d1_d2_d3_d4))
             print('D1-D2-D3-D4: '+str(len(df_sites_d1_d2_d3_d4)))
             print(df_sites_d1_d2_d3_d4)
  
     ###############################################################################################33
#             print('--------------4---------------------Comparo D1 y D4 primero D2 y D3 despues')
  
             ## Busqueda de sitios conservados entre den1 y den4
  
             sites_d1_d4=[]
  
             for i in range(len(dfs[0])):
                 for j in range(len(dfs[3])):
                     dist_sites = distancia(dfs[0],i,dfs[3],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[0],i,dfs[3],j)
  
                         info_sites=[]
                         info_sites.append(dfs[0]['number'].iloc[i])
                         info_sites.append(dfs[3]['number'].iloc[j])
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d1_d4.append(info_sites)
  
             df_sites_d1_d4 = pd.DataFrame(sites_d1_d4, columns=['Site d1', 'Site d4', 'Distance', 'x', 'y', 'z'])
#             print('D1-D4: '+str(len(df_sites_d1_d4)))
  
             ## Busqueda de sitios conservados entre den2 y den3
  
             sites_d2_d3=[]
  
             for i in range(len(dfs[1])):
                 for j in range(len(dfs[2])):
                     dist_sites = distancia(dfs[1],i,dfs[2],j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (dfs[1],i,dfs[2],j)
  
                         info_sites=[]
                         info_sites.append(dfs[1]['number'].iloc[i])
                         info_sites.append(dfs[2]['number'].iloc[j])
                         info_sites.append(dist_sites)
  
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d2_d3.append(info_sites)
  
             df_sites_d2_d3 = pd.DataFrame(sites_d2_d3, columns=['Site d2', 'Site d3', 'Distance', 'x', 'y', 'z'])
  
             print('D2-D3: '+str(len(df_sites_d2_d3)))
  
             ## Busqueda de sitios en los 4 serotipos
  
             sites_d1_d2_d3_d4=[]
  
             for i in range(len(sites_d1_d4)):
                 for j in range(len(sites_d2_d3)):
                     dist_sites = distancia(df_sites_d1_d4,i,df_sites_d2_d3,j)
                     if dist_sites < cutoff_dist:
  
                         punto_medio_site = punto_medio (df_sites_d1_d4,i,df_sites_d2_d3,j)
  
                         info_sites=[]
                         info_sites.append(df_sites_d1_d4['Site d1'].iloc[i])
                         info_sites.append(df_sites_d1_d4['Site d4'].iloc[i])
                         info_sites.append(df_sites_d2_d3['Site d2'].iloc[j])
                         info_sites.append(df_sites_d2_d3['Site d3'].iloc[j])
  
                         info_sites.append(dist_sites)
                         for k in range(3):
                             info_sites.append(punto_medio_site[k])
  
                         sites_d1_d2_d3_d4.append(info_sites)
  
  
             df_sites_d1_d2_d3_d4 = pd.DataFrame(sites_d1_d2_d3_d4, columns=['Site d1', 'Site d4', 'Site d2', 'Site d3', 'Distance', 'x', 'y', 'z'])
             
            ## Eliminar los sitios repetidos

             todos_distintos=True
             c=0
             loc_repetidos=[]
             while todos_distintos and c < len(df_sites_d1_d2_d3_d4):
                 for k in range (len(df_sites_d1_d2_d3_d4)-1):
                     for j in range (k+1, len(df_sites_d1_d2_d3_d4)):

                         dist_sites = distancia(df_sites_d1_d2_d3_d4,k,df_sites_d1_d2_d3_d4,j)
                         if dist_sites < cutoff_repetidos:
#                             print(system+" Den"+str(i)+" tiene sitios a menos de "+str(cutoff_repetidos)+' en '+atom+": "+str(k+1)+" y "+str(j+1))
                             loc_repetidos.append(k)
                             todos_distintos=False

                         else:
                             c+=1

             if todos_distintos == False:
                 for dataframe in [df_sites_d1_d2_d3_d4]:
                     dataframe = dataframe.drop(index=loc_repetidos, inplace=True)
             else:
                  print("Son todos sitios distintos")
            
             conservados[z].append(len(df_sites_d1_d2_d3_d4))
             print('D1-D2-D3-D4: '+str(len(df_sites_d1_d2_d3_d4)))
             print(df_sites_d1_d2_d3_d4)
  
             z+=1

             # Define the input and output file names
             input_file = 'ref_'+system[:-1]+'.pdb'  # File to copy content from
             output_file = 'prot_'+atom+'_'+system[:-1]+'.pdb'  # File to write the new content to


             # Open the input file to read and the output file to write
             # Read the input file
             with open(input_file, 'r') as file:
                 lines = file.readlines()

                 # Remove the last line que es el END del pdb
                 lines = lines[:-1]
                 for i in range (len(df_sites_d1_d2_d3_d4)):
                     for serotipos in ['Site d1']: #, 'Site d2', 'Site d3', 'Site d4']:
                         d_sites = df_sites_d1_d2_d3_d4[serotipos].tolist()
                         lines.append(leer_linea('/home/usuario/ns3/inhibidores/md_mixsv/den'+serotipos[-1]+'/'+system+'all_aligned/'+atom+'/watcent.pdb',d_sites[i]-1))
                         lines.append("TER\n")

                 lines.append("END")

             # Write to the output file
             with open(output_file, 'w') as file:
                 file.writelines(lines)
             file.close()


     ###############################################################################################
    print(conservados)
  
  
