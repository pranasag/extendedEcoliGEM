# Pipeline for extending conventional genome-scale models 

[Pranas Grigaitis](mailto:p.grigaitis@vu.nl), 2019-2020

[Dept. Modelling of Biological Processes](https://www.cos.uni-heidelberg.de/index.php/u.kummer?l=), BioQuant/COS Heidelberg, Heidelberg University, Heidelberg, Germany &
[Systems Biology Lab](http://www.teusinkbruggemanlab.nl/), AIMMS, Vrije Universiteit Amsterdam, Amsterdam, the Netherlands

---
## SBML3 model
Fully working SBML3 FBCv2 model generated with this pipeline is available via [this link](https://surfdrive.surf.nl/files/index.php/s/4PuvHz482w7OuP7) or per request via email (the file is too large to be uploaded here).

---
## Dependencies
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for search of homologous protein sequences using the SwissProt protein sequence DB
* [CBMPy](https://github.com/SystemsBioinformatics/cbmpy) and its dependencies. 

For a fresh user, the following commands should be applied to install the required packages through `apt` (Debian-based GNU/Linux), as well as `pip` (Python packages):

`# apt install ncbi-blast+ python-pip`

`$ update_blastdb --decompress swissprot`

`$ pip install --user numpy pandas scipy sympy cbmpy python-libsbml matplotlib xlrd xlwt`

Also, a LP solver is required. Either [GLPK](https://www.gnu.org/software/glpk/) or [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) are supported by CBMPy. GLPK is open-source, and academic customers can acquire CPLEX for free. 
Note: for using GLPK, please follow the instructions on installation, provided in the CBMPy [repository](https://github.com/SystemsBioinformatics/cbmpy).

---
## Setup
The `main.py` pipeline also contains all the information, needed to get the proteome information and kinetic data. The following organism-dependent parameters must be specified:

`organism =`: the name of biological species;
`proteomeID =`: the proteome ID from the [UniProt proteome database](https://www.uniprot.org/proteomes/). 
Note: the `kcats.py` file uses the proteome ID to form a string `proteome`. This might be needed to be modified when working with different species than *Escherichia coli*.

Also, the pipeline uses [BRENDA](https://brenda-enzymes.org/) database to fetch the kinetic data. For use, one has to register (free of charge) and provide the user e-mail and password as the following variables:

`brendaEmail =`: the registration email of the user;
`brendaPassword =`: the password.
Note: due the fact that one has to submit the password as plain text, using these scripts should be performed in a strictly-local directory. 

---
## Usage
The main pipeline is started using:
`$ python main.py <modelLocation>`, where `<modelLocation>` is the location of the stoichiometric model in SBML (`.xml`) format. The output will be provided in the same folder.

Depending on the machine, the process might take some hours, but has to be run once to generate a working model.
