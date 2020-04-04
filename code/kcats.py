#!/bin/python
"""
Script for fetching kcat values from the target organism and homologous enzymes for 
a proteome-constrained model

Pranas Grigaitis, 2019-2020
Dept. Modelling of Biological Processes, BioQuant/COS Heidelberg, Heidelberg University &
Systems Biology Lab, AIMMS, VU Amsterdam
pranas.grigaitis@bioquant.uni-heidelberg.de
p.grigaitis@vu.nl
"""
import hashlib
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
import requests
import subprocess
import sys
from SOAPpy import SOAPProxy

def BRENDAunifiedCall(email, password, endpointURL, client, EC, organism, homologs):
	BRENDAparameters = '%s,%s,%s*%s#'%(email, password, "ecNumber", EC)
	reply = None
	while reply is None: reply = client.getTurnoverNumber(BRENDAparameters)
	kcatValues = reply.split('!')
	kcatValues = [x for x in kcatValues if re.search('[Ww]ild-type|[Ww]ild type|[Ww][Tt]', x)]
	
        kcatOrganism = [x for x in kcatValues if re.search(organism, x)]
	kcatHomologs = []
	for h in homologs: kcatHomologs = kcatHomologs + [x for x in kcatValues if re.search(h, x)]
	
	turnoverValues = ProcessBRENDAvalues(kcatOrganism)
	print(turnoverValues)	
	if len(turnoverValues) > 0: 
		return turnoverValues
	else: 
		turnoverValues = ProcessBRENDAvalues(kcatHomologs)
		if len(turnoverValues) > 0: return turnoverValues
		else: 
			turnoverValues = ProcessBRENDAvalues(kcatValues)
			return turnoverValues

def ProcessBRENDAvalues(listOfValues):
	turnoverValues = []
	listOfValues = ''.join(listOfValues).split('#')
	for kcat in listOfValues:
		if 'turnoverNumber*' in kcat: turnoverValues.append(float(kcat.replace('turnoverNumber*', '')))
	turnoverValues = [x for x in turnoverValues if x > 0.0]
	
	return turnoverValues

# BRENDA parameters
email = sys.argv[4] 
password = str(sys.argv[5])
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client = SOAPProxy(endpointURL)

"""
Acquiring the proteome of interest:
Please visit UniProt and collect both the proteome ID from https://www.uniprot.org/proteomes/, as well as the species name in single and double quotations (required for the REST API of SABIO-RK).
Please check that you set to download the *reference* proteome and not the pan-proteome.
We will acquire both FASTA of the proteome (for BLAST against SwissProt DB), as well as the list UniProt identifiers (for matching UniProts with gene names and EC numbers).

"""
filename = sys.argv[1]
organism = sys.argv[2]
proteomeID = sys.argv[3]



proteome = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/%s_83333.fasta.gz" %(proteomeID)

os.system("wget -Nq %s" %(proteome))
os.system("gunzip -f %s" %(proteome.split('/')[-1]))

# Get current enzyme nomenclature from EXPASY ENZYME
os.system('wget -Nq ftp://ftp.expasy.org/databases/enzyme/enzyme.dat')

with open('enzyme.dat', 'r') as f:
	enzymeDB = f.readlines()
	f.close()

"""
Part I: BLASTing for homologous enzymes against SwissProt DB
proteome.fasta > proteome.blast.out
(Takes some 30 minutes on a regular laptop)

We will need this data for (a) determination of kinetic values, as well as (b) filling in the EC values from UniProt that were not annotated for the original organism.

Here, BLAST flag -num_threads is set to (#logical cores - 1) on your machine. By default num_threads = 1, which makes it rather inefficient. 
"""
blastSubprocess = os.system('blastp -db swissprot -query %s -out %s -num_threads %d -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scomnames qcovs"' %((proteome.split('/')[-1]).replace('.gz', ''), (proteome.split('/')[-1]).replace('.gz', '.blast.out'), (multiprocessing.cpu_count() - 1)))

blastResult = pd.read_csv(proteome.split('/')[-1].replace('.gz', '.blast.out'), names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'scomnames', 'qcovs'], sep = '\t')

"""
Filtering by pident > 70 and qcovs > 70 to yield only very homologous entries (Aloy et al., 2003)
"""
blastResult = blastResult[(blastResult['pident'] > 70.0) & (blastResult['qcovs'] > 70.0)]

"""
Part II: identification of EC numbers for both the original proteome and homologous proteins
"""
# Get UniProtIDs from the current proteome: 
UniProtIdentifiers = pd.read_csv("https://www.uniprot.org/uniprot/?query=proteome:%s&format=tab&columns=id,genes(OLN),ec" %(proteomeID), sep = '\t')

# Match UniProts with EC numbers.
ECnumbers = [x.replace('\n', '') for x in enzymeDB if (x[:2] == 'ID') | (x[:2] == 'DR')]
ECentries = {}
for e in range(0, len(ECnumbers)):
	if ECnumbers[e][:2] == 'ID': 
		key = ECnumbers[e]
		entries = []
		if e == len(ECnumbers) - 1: 
			ECentries.update({key.replace('ID   ', ''):''.join(entries)})
			continue
		e = e + 1
		while not ECnumbers[e][:2] == 'ID':
			entries.append(ECnumbers[e])
			e = e + 1
		ECentries.update({''.join(entries).replace('DR   ', ''):key.replace('ID   ', '')})

ECUniProts = ECentries.keys()

for identifier in range(0, len(UniProtIdentifiers)):
	homologs = list(blastResult['sseqid'][blastResult['qseqid'].str.contains(UniProtIdentifiers['Entry'][identifier])])
	ECs = []
	for h in range(0, len(homologs)): 
		keys = [x for x in ECUniProts if homologs[h].split('|sp|')[1].split('.')[0] in x]
		for k in keys: ECs.append(ECentries[k])
	UniProtIdentifiers['EC number'][identifier] = '; '.join(pd.Series(ECs).unique())

UniProtIdentifiers['kcat'] = np.nan

"""
Part III: ECs with kcat values missing are probed at BRENDA (Jeske et al., 2019). Parsing routine constructed with help of instructions on BRENDA website.
"""
for entry in range(0, len(UniProtIdentifiers)):
	turnoverValues = []
	if UniProtIdentifiers['EC number'][entry] is np.nan: continue
	if not np.isnan(UniProtIdentifiers['kcat'][entry]): continue
	homologs = blastResult[blastResult['qseqid'].str.contains(UniProtIdentifiers['Entry'][entry])]
	homologsList = list(homologs['scomnames'])
	for h in homologsList:
		if len(h.split(' ')) >2: homologsList.append(' '.join(h.split(' ')[:2]))
	homologsList = list(pd.Series(homologsList).unique())
	for EC in UniProtIdentifiers['EC number'][entry].split('; '):
		turnoverValues = turnoverValues + BRENDAunifiedCall(email, password, endpointURL, client, EC, organism, homologsList)
	if len(turnoverValues) == 0: continue
	UniProtIdentifiers.loc[entry, 'kcat'] = np.median(turnoverValues)

"""
Part IV: estimation of median kcats over different classes of enzymes. We take already known kcat values from the list we have now and compute a median for a class instead of having a global dummy value.
"""
UniProtIdentifiers.to_csv(filename, index=False)
for ECclass in range(1, 8):
	classEnzymes = UniProtIdentifiers[UniProtIdentifiers['EC number'].str.contains('^%d.' %(ECclass), na = False)]
	classMedian = np.median(classEnzymes['kcat'].dropna())
	for identifier in classEnzymes.index:
		if np.isnan(UniProtIdentifiers['kcat'][identifier]): UniProtIdentifiers.loc[identifier, 'kcat'] = classMedian

UniProtIdentifiers['kcat(h-1)'] = UniProtIdentifiers['kcat'] * 3600.0
UniProtIdentifiers.to_csv(filename, index=False)
