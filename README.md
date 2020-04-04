# Pipeline for extending conventional genome-scale models 

[Pranas Grigaitis](mailto:p.grigaitis@vu.nl), 2019-2020

Dept. Modelling of Biological Processes, BioQuant/COS Heidelberg, Heidelberg University, Heidelberg, Germany &
Systems Biology Lab, AIMMS, VU Amsterdam, Amsterdam, the Netherlands

---
## Dependencies
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for search of homologous protein sequences using the SwissProt protein sequence DB
* [CBMPy](https://github.com/SystemsBioinformatics/cbmpy) and its dependencies. 

For a fresh user, the following sequence of apt/dnf and Python packages could be installed:

`# apt install ncbi-blast+ python-pip`

`$ update_blastdb --decompress swissprot`

`$ pip install --user numpy pandas scipy sympy cbmpy python-libsbml matplotlib xlrd xlwt`

Also, a LP solver is required. Either [GLPK](https://www.gnu.org/software/glpk/) or [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) are supported by CBMPy. GLPK is open-source, and CPLEX is free for academic customers. Note: for using GLPK, please follow the instructions on installation, provided in the CBMPy [repository](https://github.com/SystemsBioinformatics/cbmpy).

---
## Setup
The `main.py` pipeline also contains all the information, needed to get the proteome information and kinetic data. The following organism-dependent parameters must be specified:

`organism =`: the name of biological species;
`proteomeID =`: the proteome ID from the [UniProt proteome database](https://www.uniprot.org/proteomes/).

Also, the pipeline uses [BRENDA](https://brenda-enzymes.org/) database to fetch the kinetic data. For use, one has to register (free of charge) and submit the user e-mail and password to the following variables:

`brendaEmail =`: the registration email of the user;
`brendaPassword =`: the password.

---
## Usage
The main pipeline could be run using:
`$ python main.py <modelLocation>`, where `<modelLocation>` is the location of the stoichiometric model in SBML (`.xml`) format. The output will be provided in the same folder.

Depending on the machine, the process might take up to 12 hours, but has to be run once.
