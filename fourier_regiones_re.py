
#IMPORTAMOS LIBRERIAS


import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RectBivariateSpline
import pandas as pd
import matplotlib.colors as mcolors
import os
from scipy import stats
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from os.path import expanduser
home = expanduser("~")
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.gridspec
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg 
import matplotlib.pyplot as plt 
#from PIL import Image
from matplotlib import rcParams
import logging
logging.basicConfig(level=logging.ERROR)

#TABLA CON LISTA DE SNAPSHOTS, EDADES, LOOKBACK
#datos_edades=pd.read_csv("/home/bego/GARROTXA/datos_GARROTXA_resim/edades.csv", sep = ",",index_col = 0)
#datos_edades = pd.read_csv("/home/bego/GARROTXA/datos_GARROTXA_resim/edades.csv", sep = ",",index_col = 0)
path_datos = "datos_GARROTXA_resim/"
path_csv = "/media/temp/bego/snapshots_resim/"
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)


def fourier(X,Y, peso=None, modo =1, maxmode=6):
    R = np.sqrt(X**2+Y**2)
    limit_p = np.percentile(R, 90) 
    nbins = 15
    threshold_noise = 1.1
    upper_bin = round(13*nbins/15) #We will take only the intermediate bins for the estimator
    bottom_bin = round(2*nbins/15)
    AA =  np.zeros((7,nbins))
    armangle =  np.zeros((7,nbins))
    steparm = limit_p/nbins
    radarm = np.arange(0,limit_p +steparm, steparm)

    A = [0,0,0,0,0,0,0]
    B = [0,0,0,0,0,0,0]
    dd = np.sqrt(X**2 + Y**2)
    indr=np.digitize(dd,radarm)-1
    nparticles = np.zeros(nbins)

    for m in range(0,maxmode+1):
       #Iterating over radial bins
       for i in range(nbins):
            X_i=X[indr==i] 
            Y_i=Y[indr==i]

            A[m] = 0
            B[m] = 0
            a = []
            b = []
            a = np.arctan2(Y_i, X_i)
            b= np.arcsin(Y_i/np.sqrt(Y_i**2 + X_i**2))
            if peso is None:
                if m == 0:
                    A[m] = np.sum(np.cos(m*a))
                    B[m] = np.sum(np.sin(m*a))
                else :
                    A[m] = np.sum(2*np.cos(m*a))
                    B[m] = np.sum(2*np.sin(m*a))

            else:
                peso_i=peso[indr==i]
                if m ==0:
                    A[m] = np.sum(peso_i*np.cos(m*a))
                    B[m] = np.sum(peso_i*np.sin(m*a))
                else :
                    A[m] = np.sum(2*peso_i*np.cos(m*a))
                    B[m] = np.sum(2*peso_i*np.sin(m*a))
            
            AA[m,i] = np.sqrt(A[m]**2+ B[m]**2)
            if m > 0:
                armangle[m,i] = np.arctan2(B[m],A[m])
            elif m == 0:
                armangle[m,i] = 0
                nparticles[i]= len(a)


   # modo1 = AA[1,:]/nparticles
   # modo2 = AA[2,:]/nparticles
   # modo3 = AA[3,:]/nparticles
   # modo4 = AA[4,:]/nparticles
   # modo5 = AA[5,:]/nparticles
   # modo6 = AA[6,:]/nparticles
    print("Angle ", np.std(armangle[modo,bottom_bin:upper_bin]),np.pi/6)
    if np.std(armangle[modo,bottom_bin:upper_bin])< np.pi/6:
        mode_amplitude = 0
    else:
        # amplitude_modo1 = np.mean(AA[1,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        # amplitude_modo2 = np.mean(AA[2,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        # amplitude_modo3 = np.mean(AA[3,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        # amplitude_modo4 = np.mean(AA[4,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        # amplitude_modo5 = np.mean(AA[5,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        # amplitude_modo6 = np.mean(AA[6,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin])
        amplitude_modo1 = np.mean(AA[1,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        amplitude_modo2 = np.mean(AA[2,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        amplitude_modo3 = np.mean(AA[3,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        amplitude_modo4 = np.mean(AA[4,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        amplitude_modo5 = np.mean(AA[5,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        amplitude_modo6 = np.mean(AA[6,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])
        if (amplitude_modo1 > threshold_noise*amplitude_modo6)& (amplitude_modo1 > threshold_noise*amplitude_modo5):
            if (amplitude_modo1 > threshold_noise*amplitude_modo3)& (amplitude_modo1 > threshold_noise*amplitude_modo4):
                #mode_amplitude = np.mean(AA[modo,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin])  
               # mode_amplitude = np.mean(AA[modo,bottom_bin:upper_bin])/np.mean(AA[0,bottom_bin:upper_bin]) 
              #  mode_amplitude = np.mean(amplitude_modo1)
                mode_amplitude = amplitude_modo1
               # np.mean(AA[modo,bottom_bin:upper_bin]/nparticles[bottom_bin:upper_bin]) 
            else:
                mode_amplitude = 0
        else:
            mode_amplitude = 0


    return mode_amplitude


def amplitudes_into_table(weight, amplitude_data): #weight: dens vphi or vr
    fourier_tabla = pd.read_csv( path_datos +f"fourier_regiones_{weight}_{etiqueta}.csv" , index_col = 0, sep = ",")
    new_row ={'Snapshot':name, "Lookback":lb , 
              "R1":amplitude_data[0], "R2":amplitude_data[1],
              "R3":amplitude_data[2], "R4":amplitude_data[3],
              "R5":amplitude_data[4], "R6":amplitude_data[5],
              "R7":amplitude_data[6], "R8":amplitude_data[7],
              "R9":amplitude_data[8], "R10":amplitude_data[9],
              "R11":amplitude_data[10], "R12":amplitude_data[11]
             }

    fourier_tabla = fourier_tabla.append(new_row, ignore_index = True)
    fourier_sorted = fourier_tabla.sort_values(by=["Snapshot"], ascending = True)
    fourier_sorted.to_csv(path_datos +f"fourier_regiones_{weight}_{etiqueta}.csv", sep = ",")


#LISTA DE SNAPSHOTS PARA ANALIZAR 

# snapshots_analysis = [ 602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626,
#  629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660,
#  660, 662, 664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
#  690, 692, 694, 698, 704, 706, 708,711, 712,714, 716, 718, 720,
#  722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
#    755, 756, 758, 761,763, 764, 766,
#   768, 770, 772, 774, 776, 778, 780,
# 782, 784, 786, 788, 790, 792, 794, 797, 798,
# 801, 802, 805, 806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828,
# 830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
# 853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 886, 888,
# 890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 927, 
# 929, 930, 932, 934, 937, 939, 941,942, 944, 946, 948, 950, 952, 954,956, 958, 961, 963, 
# 965, 966, 
# 968, 970, 972, 974, 976, 979, 980, 982, 984, 989, 990, 993, 994, 996, 999]
snapshots_analysis= [886]

iniciar = 0 #inicia una tabla nueva
edades = [0,5000] #Myr
radianes = 0.2618 
#peso = "dens"
radios = [10,12]
etiqueta = f"{radios[0]}-{radios[1]}_0-5Gyr_ytRS_mean0"

if iniciar == 1:
    dato_vacio = []
    df_vacio = pd.DataFrame(dato_vacio, columns=["Snapshot", "Lookback",
                                                 "R1", "R2", "R3", "R4",
                                                "R5", "R6", "R7", "R8"
                                                ,"R9", "R10", "R11", "R12"])
    df_vacio.to_csv(path_datos + f"fourier_regiones_Vphi_{etiqueta}.csv", sep = ",")
    df_vacio.to_csv(path_datos + f"fourier_regiones_Vr_{etiqueta}.csv", sep = ",")
    df_vacio.to_csv(path_datos + f"fourier_regiones_dens_{etiqueta}.csv", sep = ",")
 
#    df_vacio2 = pd.DataFrame(dato_vacio, columns=["Snapshot", "Npart"])
#    df_vacio2.to_csv(path_datos + "Npart.csv"  , sep = ",")

for k in range(len(snapshots_analysis)):
    name = snapshots_analysis[k]
    print(name)
    #Lookback 
    lb = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
    #Load snapshot table in csv
    dfA= pd.read_csv(path_csv +"%s_stars_Rvir.csv" %name)
    dfA["Vphi"] = -dfA["Vphi"]
    df = dfA[(dfA['Age']>edades[0] )& (dfA['Age']< edades[1] ) ].copy()
    df_radio = df[(df['R']>radios[0] )& (df['R']< radios[1] ) ].copy()

    amplitude_modes_vphi = np.zeros(12, dtype=np.float32)
    amplitude_modes_vr = np.zeros(12, dtype=np.float32)
    amplitude_modes_dens = np.zeros(12, dtype=np.float32)
    particulas = []
    #ANALIZAMOS CADA REGION
    for j in range(12):
        print("region", j)
        if j == 0:  
            dfRs_r1= df_radio[(df_radio['Phi']>2*(np.pi)-radianes)].copy()
            dfRs_r2=df_radio[(df_radio['Phi']< radianes)].copy()
            frames = [dfRs_r1, dfRs_r2]
            dfRs_r= pd.concat(frames)

        else:
            dfRs_r= df_radio[(df_radio['Phi']>j*2*radianes-radianes)&(df_radio['Phi']<j*2*radianes +radianes)].copy()
            
        Y = np.array(dfRs_r["VZ"]*np.std(dfRs_r["Z"])/np.std(dfRs_r["VZ"])) #Rescale Vz
        Vphi = np.array(dfRs_r["Vphi"]- np.mean(dfRs_r["Vphi"]))
        X = np.array(dfRs_r["Z"])
        Vr = np.array(dfRs_r["Vr"])
        #   Mass = np.array(dfRs_r["Mass"])
        fourier_resultado_vphi = fourier(X,Y, peso = Vphi)
        amplitude_modes_vphi[j]=round(fourier_resultado_vphi,3) 

        fourier_resultado_vr = fourier(X,Y, peso = Vr)
        amplitude_modes_vr[j]=round(fourier_resultado_vr,3)

        fourier_resultado_dens = fourier(X,Y)
        amplitude_modes_dens[j]=round(fourier_resultado_dens,3)


    amplitudes_into_table("dens", amplitude_modes_dens)
    amplitudes_into_table("Vphi", amplitude_modes_vphi)
    amplitudes_into_table("Vr", amplitude_modes_vr)




