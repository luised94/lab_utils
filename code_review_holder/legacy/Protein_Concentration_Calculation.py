# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:19:33 2020

@author: liusm
"""

import pandas as pd

from scipy.constants import Avogadro
import os
os.chdir("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/")

path_to_Median_Abundance_Data = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/2020-06-04_Protein_Median_Abundance.xlsx"

df = pd.read_excel(path_to_Median_Abundance_Data)
df = df[df.iloc[:,4] != "None"]


# type(df.iloc[1,4]) Check type of column. Shows that it is imported as character. Next line converts it to 
df.iloc[:,[3,4]] = df.iloc[:,[3,4]].astype('int64')
df["Median_Abundance"].mean()


#BIONUMBERS used to calculate certain statistics 
#ID: 110465 G1 phase 10-50µm^3: S phase 20-60µm^3: G2 phase 40-80µm^3 M phase 60-100µm^3 µm^3
#ID: 108315 Size and composition of yeast cells 70um^3 haploid, 120um^3 diploid
#ID: 101396 Volume of Nucleus 3.3 µm^3
#ID: 104708 Fraction of nuclear volume out of cell volume 7 % volume is 4.9um^3
#ID  110468 Fraction of nucleolus out of nucleus volume  20 %
#ID  106356 Volume of nucleus 10.3 μm^3 Range: ±3.7 μm^3
#ID  107660 Nuclear volume of yeast spore cell  0.28 μm^3
#Volume of nucleus ranges from 3.3um^3 to 13um^3
cell_Size = (70 * pow(10, -12)) * pow(10, -3) #in liters
nucleus_Size = (4.9 * pow(10, -12)) * pow(10,-3)#in liters
#Calculate concentration of proteins in whole cell (not excluding volumes) and nucleus
#Not specific to proteins actual subcellular localization thats why both are calculated
df["Concentration_in_Cell"] = ((df.iloc[:,4] / Avogadro) / cell_Size ) * pow(10, 9)#nanomolar
df["Concentration_in_Nucleus"] = ((df.iloc[:,4] / Avogadro) / nucleus_Size ) * pow(10, 9) #nanomolar
df.to_excel('2020-06-04_Protein_Concentration.xlsx', index=False)
df["Concentration_in_Cell"].mean()
df["Concentration_in_Nucleus"].mean()

#Following outputs concentration for given proteins by three letters
df[df['Protein_Name'].str.contains('Orc')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Mcm')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Cdc')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Dpb')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Psf')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Sld')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Rfa')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Rpb')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Toa')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Spt')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')
df[df['Protein_Name'].str.contains('Dbf')].iloc[:,[0,3,7]].sort_values(by='Protein_Name')


df.plot(kind='density',xlim =[-5000,5000])
