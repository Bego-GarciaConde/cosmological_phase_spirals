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
#datos_edades = pd.read_csv("tablas/edades.csv", sep = ",",index_col = 0)

radios = [10,12]
etiqueta = "all8-10"
radianes = 0.2618
edades = [0, 5000]
#tiempo_referencia = lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
rangex=[-2.5,2.5]
rangey=[-90,90]
binsx=35
binsy=35
#binsx = 45
#binsy= 45
aspect=(rangex[1]-rangex[0])/(rangey[1]-rangey[0])
#print(aspect)
deltax=(rangex[1]-rangex[0])/binsx
deltay=(rangey[1]-rangey[0])/binsy



nn = 0.8

x = [-2.5 ,-1, 0,1,2.5]
y = [80, 40, 0, -40, -80]
nbins_reg_d = [0,30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
path_csv = "/media/temp/bego/snapshots_resim/"
#path_csv = "/media/temp1/bego/snapshots/"
#path_csv = "/home/estudiantsBCN/snapshots/"
path_datos = "/home/bego/GARROTXA/datos_GARROTXA_resim/"
#path_datos = "/home/bego/GARROTXA/datos_GARROTXA/"
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)


snapshots_analysis= [ 660,666, 670, 674, 679,  684,
690,  694,  700, 706, 711, 716,  720,726, 731, 736,  740,  746, 751,
  756, 761,764, 770, 776, 780,786, 790,794, 798, 805, 810, 816, 820, 826,830, 836,
 840,  844, 850, 855, 860, 864, 870,  875,  881, 886, 890,  894, 900, 904,910, 915,  921,
 927]
# snapshots_analysis = [600, 602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 
# 629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
# 664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
# 690, 692, 694, 698, 704, 706, 708,711, 712,714, 716, 718, 720, 
# 722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
#  755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778, 780, 
# 782, 784, 786, 788, 790, 792, 794, 797, 798, 802, 805, 806, 808, 810, 812, 814, 816,
#  818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
# 853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 888,
# 890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 
# 927, 929, 930, 932, 934, 937, 939, 941,942, 944, 946, 948, 950, 952, 954,956, 
# 958, 961, 963, 965, 966, 968, 970, 972, 974, 976, 979,
#  980, 982, 984, 986, 989, 990, 993, 994, 996, 999] 
fig = plt.figure(figsize=(70,1.75*len(snapshots_analysis)))
outer = gridspec.GridSpec(1, 3, wspace=0.05, hspace=0.03)

#
for k in range(3):


    inner = gridspec.GridSpecFromSubplotSpec(nrows = len(snapshots_analysis),
                                             ncols = 12,
                    subplot_spec=outer[k], wspace=0.05, hspace=0.005)

    for j in range(len (snapshots_analysis)):
        name = snapshots_analysis[j]
        print(name)
        lookback = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
        file=path_csv +"%s_stars_Rvir.csv" %name
        dfA= pd.read_csv(file)
        dfB = dfA[(dfA['Z']< 2.7)& (dfA['Z']> -2.7 ) ].copy()
       # dfB["Vphi"] = -dfB["Vphi"]
        df = dfB[(dfB['Age']> edades[0] )& (dfB['Age']< edades[1] ) ].copy()
        df["Vphi"] = -df["Vphi"]

        df= df[(df['R']>radios[0])&(df['R']<radios[1])].copy()

        for i in range(12):
            ax = plt.Subplot(fig, inner[j,i])
            if i == 0: 
                dfRs_r1= df[(df['Phi']>2*(np.pi)-radianes)].copy()
                dfRs_r2=df[(df['Phi']<  radianes)].copy()
                frames = [dfRs_r1, dfRs_r2]
                df_red1= pd.concat(frames)
            else:
        #        df_red1= df[(df['Phi_re']>i*2*radianes-radianes)&(df['Phi_re']<i*2*radianes +radianes)].copy()
#    #     
               df_red1= df[(df['Phi']>i*2*radianes-radianes)&(df['Phi']<i*2*radianes +radianes)].copy()

            Vphi_resta = df_red1["Vphi"] - np.mean(df_red1["Vphi"])
            df_red1["Vphi_resta"] = Vphi_resta


            if k ==0:
                    z=ax.hist2d(df_red1["Z"], df_red1["VZ"], bins = [binsx, binsy],
                       norm = mcolors.PowerNorm(nn), range=[rangex, rangey],cmap = mycolormap, vmin = 0.01, vmax = 60)

                    if i == 0:
                        if j % 5==0:
                            ax.text(-14, 3, "%.2f Gyr" %lookback, fontsize = 54)
      
            if k == 1:
                    stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vphi_resta'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
                    im=np.flip(stat2.statistic.T*1.,0)
                    vphi=ax.imshow(im,cmap="cmr.guppy_r", extent=[rangex[0],rangex[1],rangey[0],rangey[1]],
                                   aspect=aspect,vmin=-10,vmax=10)

        
            if k ==2:
                    stat2=stats.binned_statistic_2d(df_red1['Z'],df_red1['VZ'],df_red1['Vr'], statistic='median', bins=(binsx,binsy), range=[rangex,rangey])
                    im=np.flip(stat2.statistic.T*1.,0)
                    vr=ax.imshow(im,cmap="seismic", extent=[rangex[0],rangex[1],rangey[0],rangey[1]],
                                   aspect=aspect,vmin=-45,vmax=45)

            ax.axis("off")

          
            fig.add_subplot(ax)

cbar_ax_dens = fig.add_axes([0.13, 0.895, 0.235, 0.005])
cbar_dens = fig.colorbar(z[3], cax=cbar_ax_dens, orientation = "horizontal")
#cbar_dens = fig.colorbar(z, cax=cbar_ax_dens, orientation = "horizontal")

cbar_dens.ax.tick_params(labelsize= 35, top= True,bottom= False,
                    labeltop=True,  labelbottom= False)
cbar_dens.set_label(label="Counts", size = 54, labelpad = -135)
cbar_dens.ax.tick_params(labelsize= 35)




cbar_ax_vphi = fig.add_axes([0.395, 0.895, 0.235, 0.005])
cbar_vphi = fig.colorbar(vphi, cax=cbar_ax_vphi, orientation = "horizontal")          
cbar_vphi.ax.tick_params(labelsize= 35, top= True,bottom= False,
                   labeltop=True,  labelbottom= False)
#cbar_vphi.ax.tick_params(labelsize= 35, top= False,bottom= True,
 #                   labeltop=False,  labelbottom= True)

cbar_vphi.set_label(label="$\mathrm{V_{\phi}}$ [km/s]", size = 54,  labelpad = -135)
cbar_vphi.ax.tick_params(labelsize= 35)


#outer[1].set_xlabel("Z (kpc)", fontsize = 30)
cbar_ax_vr = fig.add_axes([0.66, 0.895, 0.235, 0.005])
cbar_vr = fig.colorbar(vr, cax=cbar_ax_vr, orientation = "horizontal")          
cbar_vr.ax.tick_params(labelsize= 35, top= True,bottom= False,
                   labeltop=True,  labelbottom= False)

#cbar_vr.ax.tick_params(labelsize= 35, top= False,bottom= True,
 #                   labeltop=False,  labelbottom= True)
cbar_vr.set_label(label="$\mathrm{V_{R}} $[km/s]", size = 54 , labelpad = -135)
#fig.subplots_adjust(left=0.94)
cbar_vr.ax.tick_params(labelsize= 35)
fig.text(0.045, 0.5, "$\mathrm{V_{Z}}$ [km/s]", va='center', ha='center',
             rotation='vertical', fontsize=54)
fig.text(0.51, 0.105, "Z [kpc]", va='center', ha='center', fontsize=54)
plt.savefig("spirals_0-5_10-12_ytRS_fig2_corrected.png", format='png', dpi=100, bbox_inches='tight')
#fig.tight_layout()
