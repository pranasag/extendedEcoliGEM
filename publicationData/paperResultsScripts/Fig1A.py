import cbmpy
import numpy as np
import os
import sys
import pandas as pd

model_loc = sys.argv[1]
dataset_loc = sys.argv[2]
location = sys.argv[3]
model = cbmpy.CBRead.readSBML3FBC(model_loc)
lb = pd.read_csv(dataset_loc)

"""
Metabolic constraints			
"""
for i in lb['Reaction ID']: model.setReactionLowerBound(i, lb['Lower Bound'].loc[lb['Reaction ID']==i].values[0])

fba_result = cbmpy.CBCPLEX.cplx_analyzeModel(model)
if np.isnan(float(fba_result)): sys.exit(0)

fva_piece = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model, pre_opt=True, default_on_fail=True)
os.chdir(str(location+'/mops_results/'))
cbmpy.CBWrite.writeFVAdata(fva_piece[0], fva_piece[1], '%s' %(os.path.split(dataset_loc)[1].replace('.csv.fvadata.csv', '')))
