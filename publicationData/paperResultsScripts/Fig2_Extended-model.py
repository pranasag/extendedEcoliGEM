import cbmpy
import numpy as np
import os
import sys
import pandas as pd
import re

modelLoc = sys.argv[1]
growthMediumLoc = sys.argv[2]
scriptLoc = sys.argv[3]
proteomicsLoc = sys.argv[4]
resultsFolder = sys.argv[5]

model = cbmpy.CBRead.readSBML3FBC(modelLoc, scan_notes_gpr = False)
growthData = pd.read_csv(growthMediumLoc)
proteomicsData = pd.read_csv(proteomicsLoc)
resultsPath = '%s/%s' %(scriptLoc, resultsFolder)
if not (os.path.isdir(resultsPath)): os.mkdir(resultsPath)
os.chdir(resultsPath)

"""
Total protein volume constraint for E. coli
See the supplementary material of the paper for the derivation of the constraint
"""
protSum=float(0.62/0.34)
pID = 'UP000000625'
constraint = []
UniProtIDs = pd.read_csv('proteinMasses.txt', sep = '\t')

for entry in UniProtIDs.index: constraint.append([(7.3*pow(10,-4)*float(UniProtIDs['Mass'][entry].replace(',',''))), 'P_%s_synthesis' %(UniProtIDs['Entry'][entry])])

model.addUserConstraint(pid = None, fluxes = constraint, operator = '<=', rhs = protSum)

os.chdir(resultsPath)

"""
Here, we define the multiplier for the concentrations of nutrients in the growth medium. We will use this to perform glucose (and amino acid, for the supplemented MOPS variants) limitation simulations.
"""
multiplier = 1.0 #No changes in Glc abundance
for i in growthData['Reaction ID']:
	model.setReactionLowerBound(i, multiplier * growthData['Lower Bound'].loc[growthData['Reaction ID']==i].values[0])

fbaResult = cbmpy.CBCPLEX.cplx_analyzeModel(model)

fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model, pre_opt=True)
cbmpy.CBWrite.writeFVAdata(fva[0], fva[1], 'glcTitration_%s_%.2f.csv' %(os.path.split(growthMediumLoc)[1].replace('.csv', ''), j))
