#!/usr/python
import cbmpy
import copy as cp
import numpy as np
import pandas as pd
import sys
import warnings

def addAALigaseReactions(m, aT, aL, kcatDf, kcatDef):
	"""
	{
	input: CBMPy model instance, table of the amino acids in the model and the UniProt info
	on the amino acid-tRNA ligases;
	output: updated CBMPy model instance.
	}
	Creates amino acid-tRNA ligase reactions.
	"""
	kineticValues = pd.read_csv(kcatDf)
	atpSpecies = [['M_atp_c', -1.0], ['M_amp_c', +1.0], ['M_ppi_c', +1.0]]
        for a in aT.keys():
		rxnName = '%s_tRNA_ligase' %(aT[a][2].replace('L-',''))
		m.createReaction(rxnName, reversible = False, silent = True)
		for coef in [3,4]:
                        if m.getSpecies(aT[a][coef]) is None:
                                m.createSpecies(aT[a][coef], boundary = False, compartment = 'c')
                for coef in [1,3]: m.createReactionReagent(rxnName, aT[a][coef], -1.0, silent = True)
		for atp in atpSpecies: m.createReactionReagent(rxnName, atp[0], atp[1], silent = True)
		m.createReactionReagent(rxnName, aT[a][4], +1.0, silent = True)
		ligases = aL['Entry'][aL['Protein names'].str.contains(aT[a][2].replace('L-',''))]
		for ligase in ligases:
			kcatValue = kineticValues['kcat'].loc[kineticValues['Entry'] == ligase].values[0]
			if np.isnan(kcatValue): kcatValue = kcatDef
			m = addProteinToRxn(m, rxnName, ligase, kcatValue)
	
	return m

def addChaperoneRxns(m, protName, chap, proteinMW, kcat):
	"""
	{
	input: CBMPy model instance, protein UniProtID, chaperone identity, 
	protein molecular weight and the kcat of the chaperone complex;
	output: updated CBMPy model instance. 
	}
	Defines chaperone-mediated protein folding.
	"""
	foldedSpecies = 'M_%s_c' %(chap)
	usedSpecies = 'M_%s_used_c' %(chap)
	
	kcat = 1.0/(kcat * 3600.0)
	proteinMW = float(proteinMW.replace(',', ''))
	
	m.createReaction('P_%s_folding_%s' %(protName, chap), reversible = False, silent = True)
	foldingSpecies = [[foldedSpecies, -1.0 * kcat], [usedSpecies, +1.0 * kcat], ['M_%s_uf_c' %(protName), -1.0], ['M_%s_c' %(protName), +1.0]]
	for species in foldingSpecies: m.createReactionReagent('P_%s_folding_%s' %(protName, chap), species[0], species[1], silent = True) 

	atpSpecies = [['M_atp_c', -14.0], ['M_h2o_c', -14.0], ['M_adp_c', +14.0], ['M_pi_c', +14.0]]
	for a in atpSpecies: m.createReactionReagent('P_%s_folding_%s' %(protName, chap), a[0], a[1], silent = True)
	
	if (proteinMW > 70000) & (chap == 'GroE'): m.setReactionBounds('P_%s_folding_%s' %(protName, chap), 0.0, 0.0)
	
	return m
	
def addDegradationRxns(m, protName, degradationSystem, protSeq, kcat, aT):
	"""
	{
	input: CBMPy model instance, protein UniProtID, protein degradation system identity, 
	protein sequence,  kcat of the degradation system and the amino 
	acid identity table;
	output: updated CBMPy model instance. 
	}
	Defines protein degradation systems.
	"""
	foldedSpecies = 'M_%s_c' %(degradationSystem)
	usedSpecies = 'M_%s_used_c' %(degradationSystem)
	
	kcat = 1.0/(kcat * 3600.0)
	
	m.createReaction('P_%s_degradation_%s' %(protName, degradationSystem), reversible = False, silent = True)
	m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), 'M_%s_used_c' %(protName), -1.0, silent = True) 
	m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), 'M_h2o_c', -1.0 * len(protSeq), silent = True) 
	m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), foldedSpecies, -1.0 * kcat, silent = True)
	m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), usedSpecies, +1.0 * kcat, silent = True)

	atpSpecies = [['M_atp_c', -6.0], ['M_h2o_c', -6.0], ['M_adp_c', +6.0], ['M_pi_c', +6.0]]
	for a in atpSpecies: 
		try: 
			m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), a[0], a[1], silent = True)
		except: 
			m.getReaction('P_%s_degradation_%s' %(protName, degradationSystem)).setStoichCoefficient(a[0], m.getReaction('P_%s_degradation_%s' %(protName, degradationSystem)).getReagentWithSpeciesRef(a[0]).getCoefficient() + a[1])
	
	for aa in aT.keys():
                if protSeq.count(aa) == 0: continue
		m.createReactionReagent('P_%s_degradation_%s' %(protName, degradationSystem), aT[aa][1], +1.0 * protSeq.count(aa), silent = True)
	
	return m

def addProteinSynthesisRxn(m, protName, protSeq, protLength, kcat, aT): 
	"""
	{
	input: CBMPy model instance, protein UniProtID, protein synthesis, 
	protein sequence length, kcat of the ribosome (aa/s) and the amino 
	acid identity table;
	output: updated CBMPy model instance. 
	}
	Defines protein synthesis process in ribosomes.
	"""
	foldedSpecies = 'M_Ribosome_c'
	usedSpecies = 'M_Ribosome_used_c'
	kcat = float(protLength)/(kcat*3600.0)
	atpSpecies = [['M_atp_c', -4.0], ['M_h2o_c', -4.0], ['M_adp_c', +4.0], ['M_pi_c', +4.0]]
	
	m.createReaction('P_%s_synthesis' %(protName), reversible = False, silent = True)
	m.setReactionBounds('P_%s_synthesis' %(protName), 1e-9, 1e3)
	m.createReactionReagent('P_%s_synthesis' %(protName), foldedSpecies, -1.0 * kcat, silent = True)
	m.createReactionReagent('P_%s_synthesis' %(protName), usedSpecies, +1.0 * kcat, silent = True)
	m.createReactionReagent('P_%s_synthesis' %(protName), 'M_%s_uf_c' %(protName), +1.0, silent = True)
	
	for a in atpSpecies: 
		m.createReactionReagent('P_%s_synthesis' %(protName), a[0], a[1]*protLength, silent = True)
	
	for aa in aT.keys():
                if protSeq.count(aa) == 0: continue
		m.createReactionReagent('P_%s_synthesis' %(protName), aT[aa][4], -1.0 * protSeq.count(aa), silent = True)
		m.createReactionReagent('P_%s_synthesis' %(protName), aT[aa][3], +1.0 * protSeq.count(aa), silent = True)
	
	return m
			
def addProteinUsageAndDilution(m, protName):
	"""
	{
	input: CBMPy model instance and protein UniProtID;
	output: updated CBMPy model instance.
	}
	Adds protein usage and dilution reactions.
	"""	
	foldedSpecies = 'M_' + protName + '_c'
	usedSpecies = 'M_' + protName + '_used_c'
	
	m.createReaction('P_%s_usage' %(protName), reversible = False, silent = True)
	m.createReactionReagent('P_%s_usage' %(protName), foldedSpecies, -1.0, silent = True)
	m.createReactionReagent('P_%s_usage' %(protName), usedSpecies, +1.0, silent = True)
	
	return m

def addProteinToRxn(m, rxn, protName, kcat):
	"""
	{
	input: CBMPy model instance, metabolic reaction ID, protein name (UniProt identifier) and the default kcat; 
	output: updated CBMPy model instance with two new reagent species: protein (folded) and protein (used)
	}
	"""
	foldedSpecies = 'M_%s_c' %(protName)
	usedSpecies = 'M_%s_used_c' %(protName)
	
	m.createReactionReagent(rxn, foldedSpecies, (-1.0/(kcat * 3600.0)), silent = True)
	m.createReactionReagent(rxn, usedSpecies, (+1.0/(kcat * 3600.0)), silent = True)
	
	return m

def assignKineticValues(m, df):
	"""
	{input: CBMPy model instance and the name of the DataFrame of the kinetic values;
	output: updated CBMPy model instance}
	For any protein that has a kcat value obtained, this function finds the metabolic 
	reactions using this enzyme and assigns respective kinetic values.
	"""
	df = pd.read_csv(df)
	kineticValues = df[['Entry', 'kcat(h-1)']].to_dict()
	for entry in kineticValues['Entry'].keys():
		proteinEntity = kineticValues['Entry'][entry]
		kcatValue = kineticValues['kcat(h-1)'][entry]
		if not np.isnan(kcatValue):
			kcatValue = 1.0/kcatValue
			if (1.0/kcatValue) < 1e-6: kcatValue = 1e-6 
			foldedSpecies = 'M_' + proteinEntity + '_c'
			usedSpecies = 'M_' + proteinEntity + '_used_c'
			metabolicReactions = m.getSpecies(foldedSpecies).isReagentOf()
			for rxn in metabolicReactions: 
				m.getReaction(rxn).setStoichCoefficient(foldedSpecies, -1.0 * kcatValue)
				m.getReaction(rxn).setStoichCoefficient(usedSpecies, +1.0 * kcatValue)
	
	return m
	
def checkModelConsistency(m):
	"""
	{
	input: CBMPy model instance, 
	output: SBML3FBC v2 file of the model
	}
	If the model obtains an FBA solution that is not nan, writes model into a file.
	"""
	fbaSolution = cbmpy.CBCPLEX.cplx_analyzeModel(m)
	if np.isnan(fbaSolution):
		raise Exception('Error! The model is not computable anymore. Exiting.')
	else:
		cbmpy.CBWrite.writeSBML3FBCV2(m, m.name.replace('.xml','') + '.xml')
		return True

def cloneRxn(m, parentRxn, newName):
	"""
	{
	input: CBMPy model instance, the name of the reaction to be cloned, the 
	new name to be assigned; 
	output: an updated CBMPy model instance with the new reaction
	}
	Clones the reactions in parental model to avoid the mess that fbase.clone() does.
	"""
	newRxnID = (m.getReaction(parentRxn).id + '_%s' %(newName))
	newRxnName = (m.getReaction(parentRxn).name + ', %s' %(newName.replace('_', ' and ')))
	m.createReaction(newRxnID, name = newRxnName, reversible = m.getReaction(parentRxn).reversible, create_default_bounds = False, silent = True)
	m.setReactionBounds(newRxnID, m.getReaction(parentRxn).getLowerBound(), m.getReaction(parentRxn).getUpperBound())
	
	for reagent in m.getReaction(parentRxn).getStoichiometry():
		m.createReactionReagent(newRxnID, reagent[1], float(reagent[0]), silent = True)
	
	return m

def cplxFormation(m, subunits, cplx, stoich):
	"""
	{
	input: CBMPy model instance, a Pandas DataFrame with information on subunits, 
	the name of the complex, stoichiometry of the subunits and ATP demand to form the 
	complex.
	output: updated CBMPy model instance.
	}
	Defines the formation of macromolecular complexes (e.g. ribosome or proteases).
	"""	
	foldedSpecies = 'M_' + cplx + '_c'
	usedSpecies = 'M_' + cplx + '_used_c'
	for i in [foldedSpecies, usedSpecies]: 
		m.createSpecies(i, boundary = False, compartment = 'c')
	
	m.createReaction('P_%s_formation' %(cplx), reversible = False, silent = True)
	m.createReaction('P_%s_disassembly' %(cplx), reversible = False, silent = True)

	for s in subunits['Entry']:
		m.createReactionReagent('P_%s_formation' %(cplx), 'M_%s_c' %(s), (-1.0)*stoich, silent = True)
		m.createReactionReagent('P_%s_disassembly' %(cplx), 'M_%s_used_c' %(s), (+1.0)*stoich, silent = True)
	
	m.createReactionReagent('P_%s_formation' %(cplx), foldedSpecies, 1.0, silent = True)
	m.createReactionReagent('P_%s_disassembly' %(cplx), usedSpecies, -1.0, silent = True)
	
	
	return m

def createProteinSpecies(m, protName):
	"""
	{
	input: CBMPy model instance and the protein name (UniProt identifier); 
	output: updated CBMPy model instance with three new species: protein (folded), 
	protein (unfolded) and protein (used)
	}
	This function adds three new metabolites in order to add them to the subsequent reactions.
	We assume these proteins localize in the cytosol ("_c") and the following codes are used:
	"M_UniProt_c" for an active protein;
	"M_UniProt_uf_c" for the unfolded species;
	"M_UniProt_used_c" for the used species.
	"""
	foldedSpecies = 'M_' + protName + '_c'
	unfoldedSpecies = 'M_' + protName + '_uf_c'
	usedSpecies = 'M_' + protName + '_used_c'
	
	# Check whether species are not present, otherwise skip
	if (m.getSpecies(foldedSpecies) is None) & (m.getSpecies(unfoldedSpecies) is None) & (m.getSpecies(foldedSpecies) is None):
		for i in [foldedSpecies, unfoldedSpecies, usedSpecies]: m.createSpecies(i, boundary = False, compartment = 'c')
	elif (m.getSpecies(foldedSpecies) is None) | (m.getSpecies(unfoldedSpecies) is None) | (m.getSpecies(foldedSpecies) is None): 
		raise Exception('Error! Some of the protein species of UniProt ID %s are defined, some are not!' %(protName))
	
	return m
	
def getAminoAcidTable(m):
	"""
	{
	input: CBMPy model instance;
	output: a table of amino acids, their names in the model and these of 
	respective tRNA species (uncharged/charged)
	}
	Defines the amino acid table.
	"""
	AATable = {
	'A': ['ala'], 'C': ['cys'], 'D': ['asp'], 'E': ['glu'], 'F': ['phe'], 'G': ['gly'],
	'H': ['his'], 'I': ['ile'], 'K': ['lys'], 'L': ['leu'], 'M': ['met'], 'N': ['asn'],
	'P': ['pro'], 'Q': ['gln'], 'R': ['arg'], 'S': ['ser'], 'T': ['thr'], 'V': ['val'],
	'W': ['trp'], 'Y': ['tyr']
	}
	
	for a in AATable.keys():
		aaSpecies = m.getSpecies('M_%s__L_c' %(AATable[a][0]))
		if aaSpecies is None: aaSpecies = m.getSpecies('M_%s_c' %(AATable[a][0]))
		AATable[a].append(aaSpecies.id)
		AATable[a].append(aaSpecies.name)
		AATable[a].append('M_trna%s_c' %(AATable[a][0]))
		AATable[a].append('M_%strna_c' %(AATable[a][0]))
	
	return AATable

def integrateProteinsToMetNetwork(m, kcat, pID):
	"""
	{
	input: CBMPy model instance, the default kcat (in 1/s) and the proteome ID for 
	the organism (for fetching the gene-protein mappings from UniProt); 
	output: updated CBMPy model instance
	}
	Determines the proteins required for the metabolic reactions and adds them according to
	the following scheme (for each E_i and "used species" E*_i):
	A + B + 1/kcat E_i -> C + D + 1/kcat E*_i
	"""
	# For this step, we need protein annotations from the UniProtKB:
	UniProtIDs = pd.read_csv("https://www.uniprot.org/uniprot/?query=proteome:%s&format=tab&columns=id,genes,genes(OLN),ec" %(pID), sep = '\t')
	
	UniProtGenes = ' '.join(list(UniProtIDs['Gene names'])).split(' ')
	for proteinID in UniProtIDs['Entry']: m = createProteinSpecies(m, proteinID)
	rxnList = m.getReactionIds('_')
	for rxn in rxnList: 
		rxnNameSplit = rxn.split('_')
		for gene in rxnNameSplit:
			if gene not in UniProtGenes: continue
			else:
				proteinID = UniProtIDs['Entry'].loc[UniProtIDs['Gene names'].str.contains(gene)].values[0]
				m = addProteinToRxn(m, rxn, proteinID, kcat)	
	
	return m

def proteinTurnover(m, pID, kcatLocation, kcatDefault):
	"""
	{
	input: CBMPy model instance, the proteome ID, the location of the file containing 
	kcat values (for amino acid-tRNA synthethases), the default kcat (in case these are not
	available),
	output: updated CBMPy model instance
	}
	In this main function, several things are determined in parallel, which describe how 
	the protein turnover works in the model we are augmenting. These, in sequence, are: 
	- Establishing (dis)assembly of macromolecular processes;
	- Amino acid-tRNA ligases;
	- For each protein species:
	- |: Ribosomal protein synthesis
	- |  Chaperone-mediated protein folding
	- |  Protein degradation by Lon (E. coli only) or Clp proteases
	- |  "Void turnover" (folded -> used) of proteins :| 
	"""
	UniProtIDs = pd.read_csv("https://www.uniprot.org/uniprot/?query=proteome:%s&format=tab&columns=id,genes,genes(OLN),protein%snames,length,mass,sequence" %(pID, '%20'), sep = '\t')
	
	"""
	First, we define the subsets of entities, which comprise macromolecular 
	complexes (e.g. ribosomes, Lon protease) and add them one by one 
	"""
	ribosomalProteins = UniProtIDs[UniProtIDs['Protein names'].str.contains('Ribosom')]
	hsp70Proteins = UniProtIDs[UniProtIDs['Protein names'].str.contains('HSP(70|40|24)|Trigger factor')]
	groeProteins = UniProtIDs[UniProtIDs['Protein names'].str.contains('GroE(S|L)|Trigger factor')]
	lonProteins = UniProtIDs[UniProtIDs['Protein names'].str.contains('Lon protease')]
	clpProteins = UniProtIDs[UniProtIDs['Protein names'].str.contains('Clp(A|B)')]
	
	proteinComplexes = [ribosomalProteins, hsp70Proteins, groeProteins, lonProteins, clpProteins]
	cplxName = ['Ribosome', 'Hsp70', 'GroE', 'Lon', 'Clp']
	cplxStoichiometry = [1, 1, 7, 6, 6]
	
	for cplx in range(0,len(proteinComplexes)):
		m = cplxFormation(m, proteinComplexes[cplx], cplxName[cplx], cplxStoichiometry[cplx])
	
	"""
	Further we need to define the whole system of amino acid-tRNA ligases. First 
	of all, model is appended with respective molecular species, then the reactions are 
	defined.
	"""
	aminoTable = getAminoAcidTable(m)	
	aaLigases = UniProtIDs[UniProtIDs['Protein names'].str.contains('tRNA ligase')]
	aaLigases = aaLigases[~aaLigases['Protein names'].str.contains('inducible')]
	m = addAALigaseReactions(m, aminoTable, aaLigases, kcatLocation, kcatDefault)
	
	"""
	With amino acid-tRNA ligases and protein complexes in model, we can start describing 
	protein turnover in the following aspects: protein translation, protein folding, 
	protein degradation and dilution.
	kcats here are determined in units of 1/s.
	"""
	kcatComplexes = [20.0, 1.0, 1.0, 0.25, 2.0]
	
	for protein in UniProtIDs.index:
		m = addProteinSynthesisRxn(m, UniProtIDs['Entry'][protein],  UniProtIDs['Sequence'][protein], UniProtIDs['Length'][protein], kcatComplexes[0], aminoTable)
		for i in [1,2]: m = addChaperoneRxns(m, UniProtIDs['Entry'][protein], cplxName[i], UniProtIDs['Mass'][protein], kcatComplexes[i])
		for i in [3,4]: m = addDegradationRxns(m, UniProtIDs['Entry'][protein], cplxName[i], UniProtIDs['Sequence'][protein], kcatComplexes[i], aminoTable)
		m = addProteinUsageAndDilution(m, UniProtIDs['Entry'][protein])
	
	return m 
		
def splitGPRs(m):
	"""
	{
	input: CBMPy model instance; 
	output: CBMPy model instance, with split GPRs
	}
	A function for splitting GPRs in order to have gene names on rxn ID and to
	split OR relationships (isozymes) into different reactions.
	"""
	rxnList = m.getReactionIds()
	rxnList = [x for x in rxnList if m.getGPRforReaction(x) is not None]
	
	for rxn in rxnList:
		gprTree = m.getGPRforReaction(rxn).getTree()
		for key in gprTree.keys():
			if key[:3] == '_OR':
				ORKeys = gprTree[key].values()
				for ORrxn in ORKeys:
					try:
                                                gprString = []
                                                gprList = ORrxn.values()
                                                for gprValue in gprList:
                                                        try: gprString = gprString + gprValue.values()
                                                        except AttributeError: gprString.append(gprValue)
                                                gprString = '_'.join(gprString)
                                                gprString = gprString.replace('G_', '')
					except AttributeError: gprString = ORrxn.replace('G_', '')
					m = cloneRxn(m, rxn, gprString)
			elif key[:4] == '_AND':
				ANDValues = gprTree[key].values()
				try: gprString = '_'.join(ANDValues).replace('G_', '')
				except TypeError:
					gprString = []
					for value in gprTree[key].keys():
						if value[:2] == 'G_': gprString.append(value)
						else: gprString = gprString + gprTree[key][value].values() 
					gprString = '_'.join(gprString).replace('G_', '')
				m = cloneRxn(m, rxn, gprString)
                        elif key[:2] == 'G_':
				gprString = key.replace('G_', '')
				m = cloneRxn(m, rxn, gprString)
		
		m.deleteReactionAndBounds(rxn)
							
	return m
	
def splitReversibleRxns(m):
	"""
	{
	input: CBMPy model instance; 
	output: processed CBMPy model instance
	}
	Splits the reversible reactions of the model, except for the exchange ones.
	"""
	rxnsToSplit = [x for x in m.getReactionIds() if x.find('_EX_') < 0]
	m = cbmpy.CBTools.splitReversibleReactions(m, selected_reactions = rxnsToSplit)
	return m

