#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 18:22:33 2021

@author: ccastellanos
preprocesado:
- lectura del archivo
- generación de una serie de tiempo
- codificar una forma de leer todos los archivos dentro de un directorio específico y que
se apliquen las mismas rutinas
- Iterar un comando de entrada para cambiar entre archivos csv del directorio.

"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
import pandas as pd

df = pd.read_csv('data/2017-05-28.csv', delim_whitespace=True, header=59)

#PARAMERTROS QUE INTERVIENEN EN LA GEOEFECTIVIDAD

#campo magnético
B_escal = df.iloc[:, 8]
B_vec = df.iloc[:, 9]
B_lat = df.iloc[:, 10]
B_z = df.iloc[:, 16]			#coord GSM para estudio de TGM

#viento solar
vs_n = df.iloc[:, 23]
vs_v = df.iloc[:, 24]
vs_p = df.iloc[:,35]

#plasma parameter
beta = df.iloc[:,36]
alf_mach = df.iloc[:,37]
ms_mach = df.iloc[:,38]

#indices
dst = df.iloc[:,42]
kp = df.iloc[:,40]
#convertimos el DOY y el año para obtener la fecha para cada evento

from datetime import datetime, date, timedelta


df["combined"] = df.iloc[:, 0]*1000 + df.iloc[:,1] #la primer columna es el año, la segunda es el DOY. 
df["dt"] = pd.to_datetime(df["combined"], format = "%Y%j") #obtenemos una columna para la fecha

dt = df.iloc[:, 58]

df['hr'] = pd.to_timedelta(df.iloc[:,2], unit='h')
hr = df.iloc[:,59]


#comvertimos la columna de horas a datetime

df['Date'] = datetimes  = dt + hr
Date = df.iloc[:,60]

#plt.figure()
fig, (ax1, ax2, ax3,  ax4, ax5, ax6)= plt.subplots(6)
fig.suptitle('Parámetros durante la tormenta')

fig.set_figheight(18)
fig.set_figwidth(15)

ax1.plot(Date, B_escal)
ax1.set_ylabel('|B| [nT]')
ax1.grid()

ax2.plot(Date, B_vec, 'g')
ax2.set_ylabel('B [nT]')
ax2.grid()


ax3.plot(Date, B_lat, 'r')
ax3.set_ylabel('B [°]')
ax3.grid()


ax4.plot(Date, vs_n, 'c')
ax4.set_ylabel('[cm⁻³]')
ax4.grid()


ax5.plot(Date, vs_v, 'm')
ax5.set_ylabel('[kms⁻¹]')
ax5.grid()


ax6.plot(Date, vs_p, 'k')
ax6.set_ylabel('[nPa]')
ax6.grid()


#ax7.plot(Date, beta, 'r--')
#ax7.set_ylabel('Beta')
#ax7.grid()


#ax8.plot(Date, ms_mach)
#ax8.set_ylabel('[Mach]')
#ax8.grid()


#ax9.plot(Date, dst)
#ax9.set_ylabel('[nT]')
#ax9.grid()


plt.show()









