#!/bin/python
"""
Pranas Grigaitis, 2019-2020
Dept. Modelling of Biological Processes, BioQuant/COS Heidelberg, Heidelberg University &
Systems Biology Lab, AIMMS, VU Amsterdam
pranas.grigaitis@bioquant.uni-heidelberg.de
p.grigaitis@vu.nl
"""

"""
Usage: python main.py <SBML_file>
Requirements:
- Python modules: pandas, numpy, cbmpy and the dependencies of the latter.
All the dependencies are covered when installing CBMPy. 
The following commands do the trick for Debian-based systems (requires pip):
# apt install python-pip
$ pip install --user numpy scipy sympy cbmpy python-libsbml matplotlib xlrd xlwt
- BLAST+ and its SwissProt protein sequence DB.
Installation for Debian-based systems:
# apt install ncbi-blast+
$ update_blastdb --decompress swissprot
"""
from datetime import date
import hashlib
import multiprocessing
import os
import subprocess
import sys
import warnings
"""
The following lines check whether the non-base Python packages are there. 
If the test is passed, the packages are imported. Upon failure, a warning is raised 
and sys.exit() exits the script with an error code 2.
"""
present_packages = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in present_packages.split()]
required_packages = ['numpy', 'cbmpy', 'pandas', 'requests', 'SOAPpy']
if len([x for x in installed_packages if x in required_packages]) != len(required_packages):
	warnings.warn('Missing dependencies: %s. Exiting.' %(', '.join([y for y in required_packages if y not in [x for x in installed_packages if x in required_packages]])))
	sys.exit(2)
else:
	import pandas as pd
	import numpy as np
	import cbmpy

"""
If imports are successful, we can check also for the dependencies outside Python.
In our case, this is the BLAST+ and SwissProt Protein DB. 
"""
try: bp_output = subprocess.check_output(['blastp', '-h'])
except subprocess.CalledProcessError as e:
	warnings.warn('Dependency missing: blastp. Please install NCBI BLAST+ (see documentation).')
	sys.exit(e.returncode)

bp_output = subprocess.check_output(['blastdbcmd', '-list', '.'])
if 'swissprot' not in bp_output:
	warnings.warn('Dependency missing: SwissProt protein DB. Please download the database to your local machine (see documentation).')
	sys.exit(2)

"""
If these tests are also passed, now we're ready to run the pipeline itself.
First parameters we need to define are the organism name, its proteome code and 
the login information to BRENDA (www.brenda-enzymes.org). Many of these will be used to 
collect the kinetic data for different proteins of the organism of interest, as well as 
their proteome composition.
One can find the Proteome ID at UniProt (https://www.uniprot.org/proteomes).
"""
organism = 'Escherichia coli'
proteomeID = 'UP000000625'
kineticsExportName = '%s_kinetic_values_%s.csv' %(organism.replace(' ', '_'), date.today())

# BRENDA parameters
brendaEmail = ''
brendaPassword = ''
brendaPassword = hashlib.sha256(brendaPassword).hexdigest()

"""
Now we can call the kinetics pipeline - since it uses BLAST+ and gathers info from
different kinetic data sources (so far we use BRENDA only), it will take a while: 
trying to pararelize the process as much as possible...
"""
kineticsParameters = [sys.executable, 'kcats.py', kineticsExportName, organism, proteomeID, brendaEmail, brendaPassword]
kineticsPipeline = subprocess.Popen(kineticsParameters)

"""
In parallel, we then work on the stoichiometric model itself. For this we will use CBMPy, 
using the function database written to facilitate conversion of the models. 
First, we load the model.
"""
import cbmpyFunctions as cbf
model = cbmpy.CBRead.loadModel(sys.argv[1])
model.name = os.path.basename(sys.argv[1]) + '_extended'

"""
Next step in the pipeline is to process all the gene-protein-reaction (GPR) associations which 
have an OR logical operator. These reactions should be split into two individual reactions in 
order to have correct representation of the isoenzymes and/or their complexes.
Then, we split reversible metabolic reactions into two (except exchange ones).
"""
model = cbf.splitGPRs(model)
model = cbf.splitReversibleRxns(model)

"""
We should check whether the model's working after any step and if so, write the SBML file.
We do this by running FBA and checking whether there's solution instead of a 'nan'.
"""
cbf.checkModelConsistency(model)

"""
The next step is required to integrate the proteins into the metabolic reactions as reactants.
First we use a default kcat of 65.9 1/s, later on the pipeline will update these according to the 
kinetic data.
"""
kcatDefault = 65.9

model = cbf.integrateProteinsToMetNetwork(model, kcatDefault, proteomeID)

kineticsPipeline.wait()
if kineticsPipeline.returncode != 0:
	raise Exception('Kinetics pipeline exited with a non-zero code %d!' %(kineticsPipeline.returncode))
model = cbf.assignKineticValues(model, kineticsExportName)
model = cbf.proteinTurnover(model, proteomeID, kineticsExportName, kcatDefault)

"""
Since protein turnover is considered to consume ~40% of the whole energy requirement for the 
cellular growth (Lahtvee et al. 2014), we do substract 40% from the non-growth related ATP requirement (reaction 'R_ATPM').
"""
for gene in model.getGeneIds(): model.deleteGene(gene)

turnoverCoefficient = 0.4

model.setReactionBounds('R_ATPM', model.getReactionLowerBound('R_ATPM') * (1 - turnoverCoefficient), 1000.0)
cbf.checkModelConsistency(model)
