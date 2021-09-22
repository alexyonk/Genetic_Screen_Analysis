#                   Genetic Screen Analysis v1.0
#                   *Created by AY on 9/21/2021*
#                   *Last Updated on 9/22/2021*
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

#Perform a log2 transformation on the median data and convert the output numpy array into a dataframe
#Assign the dataframe to the last column
median = data.loc[:,'MEDIAN']
a = np.log2(median)
a = pd.DataFrame(a)
data = data.assign(Normalized=a.values)

#Calculate the Z score of the log2 transformation separated by each experimental value
#Assign the Z score values to the last column
Z = data[['Normalized','Experiment']].groupby('Experiment').Normalized.transform(lambda x: stats.zscore(x,ddof = 1))
data = data.assign(Zscore = Z.values)

#Calculate the survival function of the zscores to determine statistical significance
#Calculate the opposite survival function score (for negative values)
#Assign both columns to the end of the dataframe
#NOTE = df is normdist (for positive values), while df2 is 1-normdist (for negative values)
ss = stats.norm.sf(data['Zscore'])
df = pd.DataFrame(ss,columns=['Normdist'])
df2 = 1 - df['Normdist']
data = data.assign(Normdist = df.values)
data = data.assign(InverseNormdist = df2.values)

#Output the file as a CSV, skipping the index column
data.to_csv(file_path + '_Analysis.csv', index = False)