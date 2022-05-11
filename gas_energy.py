import yt
import math
from yt import YTArray

from yt.utilities.cosmology import Cosmology

co = Cosmology(hubble_constant=0.7, omega_matter=0.3, 
               omega_lambda=0.7, omega_curvature=0.0)

import numpy as np 
from yt.units import G 
import matplotlib.pyplot as plt
import os
import array


import matplotlib
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
import matplotlib.gridspec

from matplotlib import rcParams

path_csv = "/media/temp/bego/snapshots_resim/"
path_datos = "/home/bego/GARROTXA/datos_GARROTXA_resim/"

#snapshots_analysis = [600, 602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626,
#629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660,
#660, 662, 664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
#690, 692, 694, 698, 700, 702, 704, 706, 708,711, 712,714, 716, 718, 720,
#722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
# 755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778, 780,
#782, 784, 786, 788, 790, 792, 794, 797, 798,
#801, 802, 805, 806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828,
#830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
#853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 886, 888,
#890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 927, 929, 930, 932, 934, 937, 939, 941,942, 944, 946, 948, 950, 952, 954,956, 958, 961, 963, 965]

snapshots_analysis = [# 966, 968, 970, 972, 974, 976, 979, 980, 982]
 984, 989, 990, 993, 994, 996, 999]
#snapshots_analysis = [599, 600,602, 604, 606, 608, 610, 615,  617, 618,
#619, 620,622, 625, 629, 631, 633, 635, 640, 645, 650, 655, 660,
# 665, 670, 675, 680, 683, 686, 690, 693, 696, 701, 702, 703, 704, 706, 707, 708,
# 709,  712, 714, 716, 717,  718, 719, 721, 723, 724, 725, 726, 727,
# 728, 729, 731, 733, 734, 737, 738, 739, 740, 741, 742, 743, 744, 746, 747, 748,
# 749, 750, 751, 752, 753, 754, 756, 757, 758, 759, 760, 761, 762, 764, 766, 767,
# 768, 769, 772, 773, 774, 775, 777, 778, 779, 782, 784, 785, 786, 787, 788, 789,
# 791, 792, 793, 794, 796, 797, 798, 799, 800, 804, 807, 810, 813, 815, 817, 820,
# 823, 825, 827, 830, 831, 832, 833, 834, 835, 836, 837, 839, 840, 841, 842, 843,
# 844, 846, 847, 848, 849, 850, 851, 852, 854, 855, 856, 858, 859, 861, 862, 863,
# 864, 866, 867, 868, 869, 872, 874, 876, 878, 880, 882, 884, 886, 888, 890, 893,
# 895, 897, 900, 902, 903, 905, 906, 908, 911, 913, 914, 916, 918, 920, 922, 923,
# 925, 927, 930, 935, 937, 940, 941, 945, 947, 949, 951, 957, 959, 960, 962, 963,
# 965, 967, 970, 972, 975, 977, 979, 981,983, 984, 986, 988, 990, 992, 994, 996]
datos_edades=pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
#if not os.path.isfile("gas_energy.csv"):
iniciar = 0
if iniciar ==1:
    dato_vacio = []
    df_vacio = pd.DataFrame(dato_vacio, columns=['Snapshot', "Age", "Lookback", "Etot", "Epot", "Ekin",
             "Mass_in", "Mass_out", "Lz"])
    df_vacio.to_csv(path_datos + "gas_energy.csv", sep = ",")



#lookback = []
#age = []
#Etot = []
#Ekin = []
#Epot = []
#Lz_tot = []
#mass_in = []
#mass_out = []
for i in range (len(snapshots_analysis)):
    name = snapshots_analysis[i]
    print(name)
    archivo = path_csv + "Gas_%s.csv" %name
    df= pd.read_csv(archivo)
    lb = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
    edad = datos_edades.loc[datos_edades['Snapshot'] == name, 'Age'].iloc[0]
#    lookback.append(lb)i
    centro = np.loadtxt(path_datos + "center_{}.txt".format(name))
    R200 = centro[3]
    
#    age.append(edad)
    G = 4.3e-6
    masa_particula = np.array(df["Mass"])
    Mass_acu =  np.array(df["Masa_acumulada"])
    Vx = np.array(df["VX"])
    Vy = np.array(df["VY"])
    Vz = np.array(df["VZ"])
    
    X = np.array(df["X"])
    Y = np.array(df["Y"])
    Z = np.array(df["Z"])
    
    Vtot = np.sqrt (Vx**2 + Vy**2 + Vz**2)
    r = np.array(df["Distance"])
    E_pot = -G*Mass_acu*masa_particula/r
    E_k = 0.5*masa_particula*Vtot**2
    E_tot = E_pot +E_k
#    print(E_tot)

    radio = np.array(df["R"])
    vel_phi = np.array(df["Vphi"])
    vel_z =  np.array(df["VZ"])
    vel_r = np.array(df["Vr"])
    vel_radial = (X*Vx + Y*Vy +  Z*Vz)/np.sqrt(X**2 + Y**2 + Z**2)

    Lz = radio *np.sqrt(vel_phi**2 + vel_z**2) 

    df['Etot'] = E_tot
    df['Lz'] = Lz
    df["Vradial"] = vel_radial
    df["Epot"] = E_pot
    df["Ekin"] = E_k
    #Inflows and outflows in a shell
#    dist_shell = np.max (df["R"])
    df_shell = df[(df['Distance']> 0.15*R200)&(df['Distance']<0.15*R200+ 0.5)].copy()
    
    gas_in = df_shell[(df_shell['Vradial']< -50)].copy()
    gas_out =  df_shell[(df_shell['Vradial']> 50)].copy()
    
    gas_in0 = np.sum(gas_in["Mass"])
#    mass_in.append(gas_in0)
    gas_out0 = np.sum(gas_out["Mass"])
#    mass_out.append(gas_out0)
    
    
    #Energía total en disco 
    df_disk1 = df[(df['Z']< 2.5) & (df['Z']> -2.5)].copy()
    df_disk2 = df_disk1 [(df_disk1 ['R']< 0.15*R200)].copy()
    energia_tot = np.sum(df_disk2["Etot"])/np.sum(df_disk2["Mass"])
#    Etot.append(energia_tot)
    
    energia_kin = np.sum(df_disk2["Ekin"])/np.sum(df_disk2["Mass"])
#    Ekin.append(energia_kin)
    
    energia_pot = np.sum(df_disk2["Epot"])/np.sum(df_disk2["Mass"])
#    Epot.append(energia_pot)
    
    angular_momentum= np.sum(df_disk2["Lz"])/np.sum(df_disk2["Mass"])
#    Lz_tot.append(angular_momentum)
    
    energy = pd.read_csv(path_datos + "gas_energy.csv", index_col = 0, sep = ",")
    new_row ={'Snapshot':name,"Age":edad,   "Lookback": lb,
             "Etot": energia_tot, "Epot" :energia_pot , "Ekin" : energia_kin,
             "Mass_in" :gas_in0,
              "Mass_out": gas_out0, "Lz" :angular_momentum} 
    energy = energy.append(new_row, ignore_index = True)
    energy.to_csv(path_datos + "gas_energy.csv", sep = ",")


    del df
    del df_disk1
    del df_disk2
    del df_shell
    
    
#datos = {'Snapshot':snapshots_analysis,"Age":age,   "Lookback": lookback,
 #            "Etot": Etot, "Epot" :Epot , "Ekin" : Ekin, "Mass_in" :mass_in, 
 #             "Mass_out": mass_out, "Lz" :Lz_tot}


#tabla= pd.DataFrame(data=datos)
#buf = "gas_energy.csv"  
#tabla.to_csv(buf,  sep=",")
