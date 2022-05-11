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

#snapshots_analysis= [600, 602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626,
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


snapshots_analysis = [#965, 966, 968, 970, 972, 974, 976, 979, 980, 982]
982,  984, 989, 990, 993, 994, 996, 999]
path_datos = "/home/bego/GARROTXA/datos_GARROTXA_resim/"
path_csv = "/media/temp/bego/snapshots_resim/"
datos_edades = pd.read_csv(path_datos + "edades.csv", index_col = 0, sep = ",")


iniciar = 0
#if  os.path.isfile("SFH_eficicencia.csv")== False:
if iniciar == 1:
    dato_vacio = []
    df_vacio = pd.DataFrame(dato_vacio, columns=['Snapshot', 'Age', "Lookback",
          "Eficiencia","Masa", "Gas_disponible",  
        "Eficiencia_bulge","Masa_bulge", "Gas_disponible_bulge",  
        "Eficiencia_disk","Masa_disk", "Gas_disponible_disk", "Nparticulas"])
    df_vacio.to_csv(path_datos + "SFH_eficiencia.csv", sep = ",")
#else:
 #   datos_SFH = pd.read_csv("SFH_eficiencia.csv", index_col = 0, sep = ",")


age= []
lookback = []
eficiencia_SFH = []
masa_creada_SFH = []
gas_disponible_SFH = []
lista = []
for i in range(1,len(snapshots_analysis)):
    name = snapshots_analysis[i]
#    datos_SFH = pd.read_csv("SFH_eficiencia.csv", index_col = 0, sep = ",")
    print(name)
    lista.append(name)
    age1 = datos_edades.loc[datos_edades['Snapshot'] == name, 'Age'].iloc[0]
    print(i)
    name0 = snapshots_analysis[i-1]
    age0 = datos_edades.loc[datos_edades['Snapshot'] == name0, 'Age'].iloc[0]
    age.append(age0)
    lb = datos_edades.loc[datos_edades['Snapshot'] == name, 'Lookback'].iloc[0]
    lookback.append(lb)

    df= pd.read_csv(path_csv + "%s_stars_Rvir.csv" %name, sep = ",")
    centro = np.loadtxt(path_datos + "center_{}.txt".format(name))
    R200 = centro[3]

    df_r= df[(df['R']<0.15*R200)].copy()
    df_red= df_r[(df_r["Z"]>-2.7)&(df_r['Z']<2.7)].copy()
    edad_snapshots= (age1 - age0)*1000
    print("Edad entre snapshots (Gyr)")
    print(edad_snapshots)

    df_age_i = df_red[(df_red["Age"]< edad_snapshots)].copy()
    df_age_bulge = df_age_i[(df_age_i["R"]< 4)].copy()
    df_age_disk = df_age_i[(df_age_i["R"]>= 4)].copy()
    masa_creada= np.sum(df_age_i["Mass"])
    masa_creada_bulge = np.sum(df_age_bulge["Mass"])
    masa_creada_disk = np.sum(df_age_disk["Mass"])
    nparticulas = len(df_age_i["Mass"])
    masa_creada_SFH.append(masa_creada)
    del df
    df_gas = pd.read_csv(path_csv +"Gas_%s.csv" %name, sep = ",")
    df_gas_disco = df_gas[(df_gas["R"]< 0.25*R200)&(np.abs(df_gas['Z'])<2.7)].copy()
    df_gas_SF = df_gas_disco[(df_gas_disco["nH"]> 1)&(df_gas_disco['Temperature']<9000)].copy()
    del df_gas
    df_gas_bulge = df_gas_SF[(df_gas_SF["R"]< 4)].copy()
    df_gas_disk = df_gas_SF[(df_gas_SF["R"]>= 4)].copy()
    gas_disponible_bulge = np.sum(df_gas_bulge["Mass"])
    gas_disponible_disk = np.sum(df_gas_disk["Mass"])
    gas_disponible = np.sum(df_gas_SF["Mass"])
    gas_disponible_SFH.append(gas_disponible)
    eficiencia = masa_creada/(masa_creada + gas_disponible)
    eficiencia_bulge = masa_creada_bulge/(masa_creada_bulge + gas_disponible_bulge)
    eficiencia_disk = masa_creada_disk/(masa_creada_disk + gas_disponible_disk)
    eficiencia_SFH.append(eficiencia)
    del df_gas_SF
    del df_r 
    datos_SFH = pd.read_csv(path_datos + "SFH_eficiencia.csv", index_col = 0, sep = ",")
    new_row = {'Snapshot':name, 'Age': age1, "Lookback": lb,
           "Eficiencia":eficiencia, "Masa": masa_creada, "Gas_disponible": gas_disponible,
            "Eficiencia_bulge":eficiencia_bulge,"Masa_bulge": masa_creada_bulge, "Gas_disponible_bulge": gas_disponible_bulge,
            "Eficiencia_disk":eficiencia_disk,"Masa_disk": masa_creada_disk, "Gas_disponible_disk": gas_disponible_disk,
            "Nparticulas":nparticulas }#
    datos_SFH = datos_SFH.append(new_row, ignore_index = True)
    datos_SFH.to_csv(path_datos + "SFH_eficiencia.csv", sep = ",")

