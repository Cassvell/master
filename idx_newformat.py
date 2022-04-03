'''
En esta rutina tiene como propósito tomar un archivo de datos de índices 
geomagnéticos con una ventana de tiempo de varios días o anual, para 
fragmentar el archivo original, dividiendolo en varios archivos que abarcan un
día en sus respectivas ventanas de tiempo. Estos archivos pueden ser de tres 
tipos: dst, kp y sym
'''
################################################################################
#importación de librerías
################################################################################
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import re
import os
################################################################################
#introducimos la ventana de tiempo que abarca el archivo original y que está 
#indicada en el nombre del mismo
idate = input("format code:\n (yyyy-mm-dd): ")
fdate = input("format code:\n (yyyy-mm-dd): ")

code_stat = 'D' #código que indica el estado de los archivos, D, P o Q
                #(Definitivos, Provisionales o Rápidos)

file_type=input("select file type:\n(dst, kp or sym): ") #selección del tipo de
#archivo de acuerdo al índice geomagnético de interés

path='/home/c-isaac/Escritorio/proyecto/tormentas-list/rutidl/'+file_type+'/'

code_name = 0   #código que indica el tipo de índice geomagnético en el nombre 
                #del archivo
                
head = 0        #número de líneas en el encabezado
sample = 0      #tiempo de muestreo

#revisión de si existe el archivo que se busca procesar, en función de la 
#ventana de tiempo y del tipo de índice de interés seleccionado. 
if file_type == 'dst':
    if os.path.isfile(path+'Dst_'+idate+'_'+fdate+'_D'+'.dat'):
        code_stat = 'D'
        file_type = 'dst'
        code_name = 'Dst'        
        
    elif os.path.isfile(path+'Dst_'+idate+'_'+fdate+'_P'+'.dat'):
        code_stat = 'P'
        file_type = 'dst'
        code_name = 'Dst'
        head    = 24   
         
    else:
        code_stat = 'Q'
        file_type = 'dst'
        code_name = 'Dst'
        head    = 24     
        
elif file_type == 'kp':        
    if os.path.isfile(path+'Kp_'+idate+'_'+fdate+'_D'+'.dat'):
        code_stat = 'D'
        file_type = 'kp'
        code_name = 'Kp'  
        head    = 35      
         
    elif os.path.isfile(path+'Kp_'+idate+'_'+fdate+'_P'+'.dat'):
        code_stat = 'P'
        file_type = 'kp'
        code_name = 'Kp'
        head    = 35         
    else:
        code_stat = 'Q'
        file_type = 'kp'
        code_name = 'Kp' 
        head    = 35
    
elif file_type == 'sym':
    sample=input('chose sample rate: \n h or m: ')
    
    if os.path.isfile(path+'ASY_'+idate+'_'+fdate+'h_D'+'.dat'):
        if sample=='h':    
            code_stat = 'h_D'
        else:
            code_stat = 'm_D'
                        
        file_type = 'sym'
        code_name = 'ASY'
                
    elif os.path.isfile(path+'ASY_'+idate+'_'+fdate+'h_P'+'.dat'):
        if sample=='h':    
            code_stat = 'h_P'
        else:
            code_stat = 'm_P'
        file_type = 'sym'
        code_name = 'ASY'
        head    = 24
        
    else:
        if sample=='h':    
            code_stat = 'h_Q'
        else:
            code_stat = 'm_Q'    
        file_type = 'sym'
        code_name = 'ASY'
        head    = 24 
                         
else:
    print('file does not exist.\n Try again with another date or another name')        
################################################################################
################################################################################
################################################################################   
#generación del DataFrame

if file_type != 'sym':
    df = pd.read_csv(path+code_name+'_'+idate+'_'+fdate+'_'+code_stat+'.dat', \
                     header=head, delim_whitespace=True)

else:
    df = pd.read_csv(path+code_name+'_'+idate+'_'+fdate+code_stat+'.dat', \
                     header=head, delim_whitespace=True)
                     
df = df.drop(columns=['|'])
################################################################################
################################################################################
#fragmentación del archivo original en varios archivos con un día de ventana de
#tiempo
################################################################################
################################################################################

enddata = 0 #hora del del día en que debe terminar la ventana de tiempo

idx = 0     #vector de tiempo a generar

step = 0    #indica la tasa de muestreo, dependiendo del índice geomagnético 
            #seleccionado
            
if file_type == 'dst':
    enddata = fdate+ ' 23:00:00'
    idx = pd.date_range(start = pd.Timestamp(idate), \
                        end = pd.Timestamp(enddata), freq='H')
    step=24
    fhour=23                    
elif file_type == 'kp':   
    enddata = fdate+ ' 21:00:00'                     
    idx = pd.date_range(start = pd.Timestamp(idate), \
                        end = pd.Timestamp(enddata), freq='3H')
    step=8
    fhour=7
    
elif file_type == 'sym':
    if sample =='h':
        enddata = fdate+ ' 23:00:00'
        idx = pd.date_range(start = pd.Timestamp(idate), \
                            end = pd.Timestamp(enddata), freq='H')
        step=24
        fhour=23                                 
    else:
        enddata = fdate+ ' 23:59:00'
        idx = pd.date_range(start = pd.Timestamp(idate), \
                        end = pd.Timestamp(enddata), freq='T')
        step=1440
        fhour=1439                         
else:
    print('try with another sample rate')      


df = df.set_index(idx)
df = df.drop(columns=['DATE', 'TIME'])
df = df.reset_index()
                  
#print(df)                 
#en este bucle, se generan los archivos por cada día y se guardan en la dir 
#indicada.
for i in range(0,len(idx),step):
    mask = (df['index'] >= idx[i]) & (df['index'] <= idx[i+fhour])
    df_new = df[mask]
    date = str(idx[i])
    date = date[0:10]
    
    if file_type == 'dst' or file_type == 'kp':
        name_new = file_type+'_'+date+'.txt'
        new_path = '/home/c-isaac/Escritorio/proyecto/tormentas-list/rutidl/'+\
        file_type+'/daily/'
        
        df_new.to_csv(new_path+name_new, sep= ' ', index=False)  

    else:
        name_new = file_type+'_'+date+code_stat+'.txt'
        new_path = '/home/c-isaac/Escritorio/proyecto/tormentas-list/rutidl/'+\
        file_type+'/daily/'
        
        df_new.to_csv(new_path+name_new, sep= ' ', index=False)        









