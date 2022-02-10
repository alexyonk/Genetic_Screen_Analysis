#                   Genetic Screen Analysis v2.0
#                   *Created by AY on 9/21/2021*
#                   *Last Updated on 10/14/2021*
#     *For any issues or bugs, please contact alex.yonk2@gmail.com*
#*This script was designed to import excel data, calculate log2 scores, mean, std, z score, 
#and assess the statistical significance of changes in gene expression

#Install and import these libraries
from tkinter import filedialog
import pandas as pd
import numpy as np
from scipy import stats

#Allows the user to select the file and read in the excel data as a pandas dataframe
file_path = filedialog.askopenfilename()
data = pd.read_excel(file_path)
file_path = file_path[:-5]
Genes = list(set(data['Strain']))
KeyList = []
UpRegList = []
DownRegList = []
FullUpList = []
FullDownList = []

#Split the dataframe by strain, experiment, and replicate
Split = dict(iter(data.groupby(['Strain','Experiment','Replicate#'])))

#Calculate normalization, zscore, normal distribution, and inverse normal distribution
for i in range(len(Split)):
    median = np.log2(list(Split.values())[i].loc[:,'MEDIAN'])
    Z = stats.zscore(median)
    ss = stats.norm.sf(Z)
    ss1 = 1 - ss
    KeyList.append([median, Z, ss, ss1])

#Append the previous calculated values into each corresponding dataframe with the appropriate header
for i in range(len(Split)):
    list(Split.values())[i]['Normalized'] = KeyList[i][0]
    list(Split.values())[i]['Zscore'] = KeyList[i][1]
    list(Split.values())[i]['NormDist'] = KeyList[i][2]
    list(Split.values())[i]['InvNormDist'] = KeyList[i][3]


for i in range(len(Split)):
    for ii in range(len(list(Split.values())[i])):
        if list(Split.values())[i].iloc[ii,7] < 0.05:
            UpReg = list(Split.values())[i].iloc[ii]
            UpRegList.append(UpReg)
        if list(Split.values())[i].iloc[ii,8] < 0.05:
            DownReg = list(Split.values())[i].iloc[ii]
            DownRegList.append(DownReg)

Up = pd.concat(UpRegList, axis = 1)
Up = Up.sort_values(by=['Strain','Experiment','Replicate#'], axis = 1, ascending = True)
Down = pd.concat(DownRegList, axis = 1)
Down = Down.sort_values(by=['Strain','Experiment','Replicate#'], axis = 1, ascending = True)

writer = pd.ExcelWriter(file_path + '_Analysis.xlsx',engine = 'xlsxwriter')
Up.to_excel(writer, sheet_name = 'UpReg')
Down.to_excel(writer,sheet_name = 'DownReg')
for n,Split in Split.items():
    Split.to_excel(writer,sheet_name = str(n))
writer.save()