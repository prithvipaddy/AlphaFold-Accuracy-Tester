# AlphaFold-Accuracy-Tester
Python program to test the accuracy of DeepMind's 3D protein structure predictor, AlphaFold (https://alphafold.ebi.ac.uk).

The input for this program is a file with a list of UniProt IDs. An example has been added to the repo (ProteinList_HumanKinase.txt). 

The required pre-requisites for a system to run the program are as follows:
  1. Having NCBI-Blast+ installed with the 'pdbaa' database. Follow the instructions on "https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata" to download and install.
  2. Having PyMOL installed (https://pymol.org/dokuwiki/doku.php?id=installation).

The program runs a BLAST search to find all the PDB files that have at least a 70% match with each UniProt ID. Then it downloads AplhaFold's predicted model of each UniProt ID. It compares the predicted model with the existing PDB files and returns the RMSD of the two structures. 

The final output is a file (example attached: rmsdAll.txt) with the following elements in each line:
  1. UniProt ID
  2. PDB entry
  3. Match percentage of the PDB entry with respect to the UniProt ID after the BLAST search
  4. RMSD between the PDB file and the predicted structure of the UniProt ID
