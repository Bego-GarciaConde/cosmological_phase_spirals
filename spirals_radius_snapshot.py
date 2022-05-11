#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:21:52 2021

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




name = 820
nn = 0.8
#radios = [5,8,11,14,17,20]
#edades = [0, 1000, 3000, 5000, 7000, 13000]
edades = [0,5]
radios = [4,6,8,10,12,14,16,18]
radianes = 0.2618
posicion =11
#posicion = 10
rangex=[-2.8,2.8]
rangey=[-90,90]
binsx=35
binsy=35
       
aspect=(rangex[1]-rangex[0])/(rangey[1]-rangey[0])
deltax=(rangex[1]-rangex[0])/binsx
deltay=(rangey[1]-rangey[0])/binsy

fig, ax = plt.subplots(3,len(radios)-1,
                       sharex=True, sharey=True,figsize = (12,6))
fig.subplots_adjust(hspace=0, wspace=0)
path_csv = "/media/temp/bego/snapshots_resim/"
path_datos = "/home/bego/GARROTXA/datos_GARROTXA_resim/"
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
plt.rc("font", family = "serif")
# axes are in a two-dimensional array, indexed by [row, col]
lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
file=path_csv + "%s_stars_Rvir.csv" %name
dfA= pd.read_csv(file)

#X =  np.array(dfA["X"])
#Y =  np.array(dfA["Y"])
#VX =  np.array(dfA["VX"])
#VY =  np.array(dfA["VY"])


#dfA["Phi_re"] = np.mod( np.arctan2(Y, X), 2*np.pi)
 #        print(df["Phi_re"])
#dfA["Vphi"] = (-X*VY + Y*VX)/np.sqrt(X**2 + Y**2)
#dfA ["Vr_re"] = (X*VX + Y*VY)/np.sqrt(X**2 + Y**2)

lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]

for j in range(len(radios)-1):

    df = dfA[(dfA['R']> radios[j])& (dfA['R']< radios[j+1]) ].copy()
    dfRs_r= df[(df['Phi']>posicion*2*radianes-radianes)&(df['Phi']<posicion*2*radianes +radianes)].copy()

    df_red1= dfRs_r[(dfRs_r['Age']>edades[0]*1000 )&(dfRs_r['Age']<edades[1]*1000 )].copy()
    for i in range(3):
#ax = plt.Subplot(fig, inner[j,i])
        if i == 0:
#            print(df_red1["Z"])
#            print(df_red1["VZ"])

            ax[i,j].hist2d(df_red1['Z'], df_red1['VZ'],bins=[binsx,binsy], norm=mcolors.PowerNorm(nn),range=[rangex,rangey], cmap=mycolormap, vmin = 0.01, vmax = 60)#, cmap='jet'
           # ax[i,j].scatter(df_red1['Z'], df_red1['VZ'],
       #                         marker='o', s = 2, color = "black", alpha = 0.1 )
            ax[i,j].set_title("%s - %s  kpc" %(radios[j], radios[j+1]), fontsize = 12)

        if i == 1:    
            Vphi_resta = df_red1["Vphi"] - np.mean(df_red1["Vphi"])
            df_red1["Vphi_resta"] = Vphi_resta   
            stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vphi_resta'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
            im=np.flip(stat2.statistic.T*1.,0)
          #  ax[i,j].scatter(df_red1['Z'],df_red1['VZ'], marker='o', s = 2, c = df_red1['Vphi_resta'], cmap='cmr.guppy', vmin = -15, vmax = 15, alpha = 0.6)
            im1=ax[i,j].imshow(im,cmap='cmr.guppy', extent=[rangex[0],rangex[1],rangey[0],rangey[1]],aspect=aspect,vmin=-10,vmax=10)

                
        if i == 2:
            stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vr'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
            im=np.flip(stat2.statistic.T*1.,0)
            im1=ax[i,j].imshow(im,cmap='seismic', extent=[rangex[0],rangex[1],rangey[0],rangey[1]],aspect=aspect,vmin=-45,vmax=45)
          # ax[i,j].scatter(df_red1['Z'],df_red1['VZ'], marker='o', s =2, c = df_red1['Vr_re'], cmap='seismic', vmin = -40, vmax = 40, alpha = 0.6 )
#        if i == 3:
#            stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Age'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
 #           im=np.flip(stat2.statistic.T*1.,0)
 #           im1=ax[i,j].imshow(im,cmap='jet', extent=[rangex[0],rangex[1],rangey[0],rangey[1]],aspect=aspect)
     
    #    ax[i,j].set_ylim(-80, 80)
    #    ax[i,j].set_xlim(-2.7, 2.7)
fig.subplots_adjust(right=0.98)
#plt.title("%.2f Gyr" %(lookback), fontsize = 12)

fig.text(0.5, 0.04, "Z [kpc]", va='center', ha='center', fontsize=16)
#fig.text(0.5, 0.95, "Stars of 3-5 Gyr", va='center', ha='center', fontsize=28)
fig.text(0.05, 0.5, "$\mathrm{V_{Z}}$ [km/s]", va='center', ha='center', rotation='vertical', fontsize=16)
fig.text(0.5, 0.96, "Lookback time: %.2f Gyr"%lookback , va='center', ha='center', fontsize=20)
plt.savefig("/home/bego/GARROTXA/Figuras/radios_%s_ytRS.png"%name, format='png', dpi=100, bbox_inches='tight')
#plt.show()
