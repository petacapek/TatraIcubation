#About the repository:
This repository contains data from a long-term incubation experiment conducted in České Budějovice. The code analyses respiration data from root-free and root-amended soils. 

#About the soils: 
The incubation experiment used six different soils sampled in the Tatra Mountains along an elevation and vegetation gradient.

#About the experiment
Soils were incubated with or without 13C-and 15N-labelled root material under 6 different conditions simulating possible winter conditions. For details on the design and purpose of the incubation experiment, email me at petacapek@gmail.com or capekp00@jcu.cz. 

#About data
The analysed data are in the “TrLongerFull.csv” file. Raw xlsx files are in the SourceDataXLSX folder, but they are not publicly available. Please write to me if you want the data. To run the statistical analysis, you don't need the xlsx files. 

#About R scripts
The repository contains several R scripts. The statistical analysis is performed separately for root-free soils (RootFreeStat.R), SOC-derived (SOCDerivedStat.R) and root-derived (RootDerivedStat.R) respiration. These three files require data in CSV format (see above) and the R libraries listed in ‘packages.R’.  “DataPrep.R” file converts original data from xlsx files to “TrLongerFull.csv”. The file ‘FieldThermalSummary.R’ calculates various climate metrics in the Tatra Mountains. This file requires thermal data in one xlsx file. To get that file, write to me.  Lastly, the “OldScripts” folder contains two additional R scripts. These scripts are older, messy versions that are impossible for anyone other than the author to understand. I keep them as a backup. 



#About R scipts output

The R scripts generate several figures in PDF. The figures are exported from R to the “StatOut” folder, which is not in the repository. Make sure you create the folder on your laptop before you start reproducing the figures. 



ENJOY THE SCIENCE


