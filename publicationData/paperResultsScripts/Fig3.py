import cbmpy
import numpy as np
import os
import sys
import pandas as pd

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
Metabolic constraints			
"""
for i in growthData['Reaction ID']:
	model.setReactionLowerBound(i, growthData['Lower Bound'].loc[growthData['Reaction ID']==i].values[0])

"""
Proteomic constraints
Should be commented out for the "w/o proteomic constraints" condition.
"""
for i in proteomicsData['Entry']:
	if (proteomicsData['conc_mmolgDW'].loc[proteomicsData['Entry']==i].values[0] == ('#VALUE!')): continue
	elif np.isnan(float(proteomicsData['conc_mmolgDW'].loc[proteomicsData['Entry']==i].values[0])): continue
	elif (proteomicsData['conc_mmolgDW'].loc[proteomicsData['Entry']==i].values[0] == 0.0): continue
	else: model.setReactionBounds('P_%s_synthesis' %(i), float(proteomicsData['0.9conc'].loc[proteomicsData['Entry']==i].values[0]), 1000.0)

"""
Total protein volume constraint for E. coli
See the supplementary material of the paper for the derivation of the constraint
"""
protSum=float(0.62/0.34)
pID = 'UP000000625'
constraint = []
UniProtIDs = pd.read_csv('proteinMasses.txt', sep = '\t')

for entry in UniProtIDs.index: constraint.append([(7.3*pow(10,-4)*float(UniProtIDs['Mass'][entry].replace(',', ''))), 'P_%s_synthesis' %(UniProtIDs['Entry'][entry])])

model.addUserConstraint(pid = None, fluxes = constraint, operator = '<=', rhs = protSum)

fbaResult = cbmpy.CBCPLEX.cplx_analyzeModel(model)
if np.isnan(float(fbaResult)): sys.exit(0) #terminate if infeasible

fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model, pre_opt=True)
cbmpy.CBWrite.writeFVAdata(fva[0], fva[1], os.path.split(growthMediumLoc)[1])
