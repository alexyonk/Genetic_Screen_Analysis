#                   Genetic Screen Analysis v3.0
#                   *Created by AY on 9/21/2021*
#                   *Last Updated on 12/2/2021*
#     *For any issues or bugs, please contact alex.yonk2@gmail.com*
#*This script was designed to import standardized excel data, calculate log2 scores, mean, std, z score, 
#and assess the statistical significance of changes in gene expression

#Install and import these libraries
from tkinter import filedialog
import pandas as pd
import numpy as np
from scipy import stats
import copy

#Allows the user to select the file and read in the excel data as a pandas dataframe
file_path = filedialog.askopenfilename()
data = pd.read_excel(file_path)
file_path = file_path[:-5]

#Set up lists for appending via for loops
Strains = list(set(data['Strain']))
Days = list(set(data['Experiment']))
UpRegList = []
DownRegList = []
FullUpList = []
FullDownList = []
NormList = []
EVlist = []

#Split the dataframe by strain, experiment, and replicate
Normalization = dict(iter(data.groupby(['Strain','Experiment','Replicate#'])))
ZDist = copy.deepcopy(Normalization)

#Calculate normalization, zscore, normal distribution, and inverse normal distribution
#Create new dataframes within a dictionary for each strain, experiment, and replicate
#Append calculated values into the appropriate created dataframes that includes:
#Log2 transformation as Normalized
#Median Absolute Deviation as MAD
#Total median as Medi
#Z score as Zscore
#Statistical significance in the positive direction as NormDist
#Statistical significance in the negative direction as InvNormDist
for i in range(len(Normalization)):
    Log = np.log2(list(Normalization.values())[i].loc[:,'MEDIAN'])
    Median = np.median(list(Normalization.values())[i].loc[:,'MEDIAN'])
    MAD = stats.median_abs_deviation(list(Normalization.values())[i].loc[:,'MEDIAN'],axis = None)
    Z = stats.zscore(Log)
    ss = stats.norm.sf(Z)
    ss1 = 1 - ss
    list(Normalization.values())[i]['Normalized'] = Log
    list(ZDist.values())[i]['Normalized'] = Log
    list(ZDist.values())[i]['MAD'] = MAD
    list(ZDist.values())[i]['Medi'] = Median
    list(ZDist.values())[i]['Zscore'] = Z
    list(ZDist.values())[i]['NormDist'] = ss
    list(ZDist.values())[i]['InvNormDist'] = ss1

#Grab all EV information for each strain, experiment, replicate dataframe into a new list
for i in range(len(Normalization)):
    for ii in range(len(list(Normalization.values())[i])):
        EV = list(Normalization.values())[i][(list(Normalization.values())[i]['Bacterial  Clone'] == 'Empty Vector')]
    EVlist.append(EV)
EVMean = copy.deepcopy(EVlist)

#Calculate the EV average from each group
for i in range(len(EVlist)):
    EVMean[i] = np.nanmean(EVlist[i].loc[:,'Normalized'])

#Calculate the log2 fold change for each strain, experiment, replicate
#Append the fold change back into the normalization dictionary
for i in range(len(Normalization)):
    FC = (list(Normalization.values())[i].loc[:,'Normalized'] / EVMean[i])
    list(Normalization.values())[i]['Fold Change'] = FC

#Using total median, MAD, and indivdiual values, calculate and append the Robust Z score values back into each list
for i in range(len(ZDist)):
    Res = 0.6745 * (list(ZDist.values())[i].loc[:,'MEDIAN'] - list(ZDist.values())[i].loc[:,'Medi']) / list(ZDist.values())[i].loc[:,'MAD']
    list(ZDist.values())[i]['RZscore'] = Res

writer = pd.ExcelWriter(file_path + '_Analysis.xlsx',engine = 'xlsxwriter')
for n,Normalization in Normalization.items():
    Normalization.to_excel(writer,sheet_name = str(n))
writer.save()






'''
Results = pd.concat(ZSplit.values())
NSplit = dict(iter(ans.groupby(['Experiment','Replicate#'])))









for i in range(len(ZSplit)):
    for ii in range(len(list(ZSplit.values())[i])):
        if list(ZSplit.values())[i].iloc[ii,7] < 0.05:
            UpReg = list(ZSplit.values())[i].iloc[ii]
            UpRegList.append(UpReg)
        if list(ZSplit.values())[i].iloc[ii,8] < 0.05:
            DownReg = list(ZSplit.values())[i].iloc[ii]
            DownRegList.append(DownReg)

Up = pd.concat(UpRegList, axis = 1)
Up = Up.sort_values(by=['Strain','Experiment','Replicate#'], axis = 1, ascending = True)
Down = pd.concat(DownRegList, axis = 1)
Down = Down.sort_values(by=['Strain','Experiment','Replicate#'], axis = 1, ascending = True)

writer = pd.ExcelWriter(file_path + '_Analysis.xlsx',engine = 'xlsxwriter')
Up.to_excel(writer, sheet_name = 'UpReg')
Down.to_excel(writer,sheet_name = 'DownReg')
for n,ZSplit in ZSplit.items():
    ZSplit.to_excel(writer,sheet_name = str(n))
writer.save()

'''