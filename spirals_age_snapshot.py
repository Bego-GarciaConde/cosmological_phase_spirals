#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:27:54 2021

@author: bego
"""

import numpy as np
import scipy.io
from scipy import stats
import matplotlib
from matplotlib import colors as mcol
from matplotlib import animation, rc
import matplotlib.pylab as plt
#from IPython.display import HTML
import matplotlib.gridspec as gridspec

import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter
import cmasher as cmr
import pandas as pd
import os
import sys
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

    
nbOfColours=257
base = matplotlib.cm.get_cmap('plasma_r')
mycolorlist = base(np.linspace(0, 1, nbOfColours))
pixelstofade=2
mycolorlist[0]=[1,1,1,1]
incrementsR = 1.*(1 - mycolorlist[pixelstofade][0])/pixelstofade
incrementsG = 1.*(1 - mycolorlist[pixelstofade][1])/pixelstofade
incrementsB = 1.*(1 - mycolorlist[pixelstofade][2])/pixelstofade
print (incrementsR, incrementsG, incrementsB)
for p in range(pixelstofade):
	n = pixelstofade-p
	mycolorlist[p][0] = mycolorlist[pixelstofade][0] + n*incrementsR
	mycolorlist[p][1] = mycolorlist[pixelstofade][1] + n*incrementsG
	mycolorlist[p][2] = mycolorlist[pixelstofade][2] + n*incrementsB
templatecmap = matplotlib.cm.get_cmap('hot')
mycolormap = templatecmap.from_list('mycustomcolormap', mycolorlist, nbOfColours)


datos_edades = pd.read_csv("edades.csv", sep = ",",index_col = 0)



name = 750
nn = 1.2
#radios = [5,8,11,14,17,20]
#edades = [0, 1000, 3000, 5000, 7000, 13000]
edades = [0,1,3,5,7]
radianes = 0.2618
posicion = 6
radios = [10, 12]
rangex=[-2.5,2.5]
rangey=[-80,80]
binsx=50
binsy=50
       
aspect=(rangex[1]-rangex[0])/(rangey[1]-rangey[0])
deltax=(rangex[1]-rangex[0])/binsx
deltay=(rangey[1]-rangey[0])/binsy

fig, ax = plt.subplots( len(edades)-1,3,
                       sharex=True, sharey=True,figsize = (6,10))
fig.subplots_adjust(hspace=0.25, wspace=0.015)
# axes are in a two-dimensional array, indexed by [row, col]
lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
file="/media/temp/bego/snapshots/%s_stars_Rvir.csv" %name
dfA= pd.read_csv(file)

X =  np.array(dfA["X_re"])
Y =  np.array(dfA["Y_re"])
VX =  np.array(dfA["VX_re"])
VY =  np.array(dfA["VY_re"])


dfA["Phi_re"] = np.mod( np.arctan2(Y, X), 2*np.pi)
 #        print(df["Phi_re"])
dfA["Vphi_re"] = (-X*VY + Y*VX)/np.sqrt(X**2 + Y**2)
dfA ["Vr_re"] = (X*VX + Y*VY)/np.sqrt(X**2 + Y**2)

lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
for j in range(len(edades)-1):

    df = dfA[(dfA['Age']> edades[j]*1000 )& (dfA['Age']< edades[j+1]*1000 ) ].copy()
    dfRs_r= df[(df['Phi_re']>posicion*2*radianes-radianes)&(df['Phi_re']<posicion*2*radianes +radianes)].copy()

    df_red1= dfRs_r[(dfRs_r['R']>radios[0])&(dfRs_r['R']<radios[1])].copy()
    for i in range(3):
        if i == 0:
            z=ax[j,i].hist2d(df_red1['Z'], df_red1['VZ'],bins=[binsx,binsy], norm=mcolors.PowerNorm(nn),range=[rangex,rangey], cmap=mycolormap)#, cmap='jet'
            ax[j,i].set_title("%s-%s Gyr" %(edades[j], edades[j+1]), fontsize = 12)

        if i == 1:    
            Vphi_resta = df_red1["Vphi_re"] - np.mean(df_red1["Vphi_re"])
            df_red1["Vphi_resta"] = Vphi_resta   
            stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vphi_resta'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
            im=np.flip(stat2.statistic.T*1.,0)
            im1=ax[j,i].imshow(im,cmap='cmr.guppy', extent=[rangex[0],rangex[1],rangey[0],rangey[1]],aspect=aspect,vmin=-15,vmax=15)

                
        if i == 2:
            stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vr_re'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
            im=np.flip(stat2.statistic.T*1.,0)
            im1=ax[j,i].imshow(im,cmap='seismic', extent=[rangex[0],rangex[1],rangey[0],rangey[1]],aspect=aspect,vmin=-30,vmax=30)

  
fig.subplots_adjust(right=0.98)

fig.text(0.5, 0.07, "Z (kpc)", va='center', ha='center', fontsize=20)
#fig.text(0.5, 0.95, "Stars of 3-5 Gyr", va='center', ha='center', fontsize=28)
fig.text(0., 0.5, "$V_{Z} (km/s)$", va='center', ha='center', rotation='vertical', fontsize=20)

plt.savefig("spirals_%s_age_%s.png"%(name, posicion), format='png', dpi=200, bbox_inches='tight')
#plt.show()
