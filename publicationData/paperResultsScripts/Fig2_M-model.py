import cbmpy
import numpy as np
import os
import re
import sys
import pandas as pd

model_loc = sys.argv[1]
dataset_loc = sys.argv[2]
location = sys.argv[3]
model = cbmpy.CBRead.readSBML3FBC(model_loc)
lb = pd.read_csv(dataset_loc)

"""
Glucose limitation simulation for the M-model. Since this is rather straightforward to implement, we can run many conditions in one script in a rather short time.
"""

"""
Metabolic constraints			
"""
for j in np.arange(0.05, 1.05, 0.05):
	for i in [x for x in list(lb['Reaction ID']) if re.search('__L_e|__D_e', x) is not None]: model.setReactionLowerBound(i, j * lb['Lower Bound'].loc[lb['Reaction ID']==i].values[0])
	for i in [x for x in list(lb['Reaction ID']) if re.search('__L_e|__D_e', x) is None]: model.setReactionLowerBound(i, lb['Lower Bound'].loc[lb['Reaction ID']==i].values[0])
	
	fba_result = cbmpy.CBCPLEX.cplx_analyzeModel(model)
	if np.isnan(float(fba_result)): sys.exit(0)
	
	fva_piece = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model, pre_opt=True, default_on_fail=True)
	os.chdir(str(location+'/glcTitrationResults/'))
	cbmpy.CBWrite.writeFVAdata(fva_piece[0], fva_piece[1], 'glcTitration_%s_%.2f' %(os.path.split(dataset_loc)[1].replace('.csv', ''), j))
