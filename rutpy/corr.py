import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc

from PIL import Image

df    = pd.read_csv('/home/c-isaac/Escritorio/proyecto/tormentas-list/rutidl/'+\
                    'output/TGM_corr_V3.txt', header=0, sep='\s+')


plt.title('Tasa de cambio en la correlación de Dst y '+\
           r'$\mathrm{\Delta H}$ [nT]')

plt.plot(df['nT'], df['TGM(General)'], label='TGM(General)', color='k', \
         linewidth=3.0)
         
plt.plot(df['nT'], df['TGM1'], label='TGM1', color='k', linewidth=1.0)
plt.plot(df['nT'], df['TGM2'], label='TGM2', color='r', linewidth=1.0)
plt.plot(df['nT'], df['TGM3'], label='TGM3', color='y', linewidth=1.0)
plt.plot(df['nT'], df['TGM4'], label='TGM4', color='m', linewidth=1.0)
plt.plot(df['nT'], df['TGM5'], label='TGM5', color='g', linewidth=1.0)
plt.plot(df['nT'], df['TGM6'], label='TGM6', color='b', linewidth=1.0)
plt.ylabel('Correlacion')
plt.xlabel('Limite para Dst (nT)')

plt.legend()
plt.grid()
plt.show()
