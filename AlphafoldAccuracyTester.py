from pymol import cmd
import os
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import subprocess
import tempfile
import smtplib

#For this program to run, the following must be installed:
#   BLAST+ with the 'pdbaa' database
#   Pymol

testInput = "ProteinList_HumanKinases.txt"

#Downloads the fasta file for a UniProt ID
def downloadUniprotFasta(uniprot):
    link_pattern = 'http://www.uniprot.org/uniprot/{}.fasta?include=yes'
    link = link_pattern.format(uniprot)

    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)

    response = requests.get(link)
    # use session.get instead of requests.get as it allows for retries and 
    # prevents max retries exceeded error
    # https://stackoverflow.com/questions/23013220/max-retries-exceeded-with-url-in-requests

    fastafile = uniprot + ".fasta"
    open("<DESTINATION PATH>"+fastafile, "wb").write(response.content)

#Downloads the alphafold predicted model of a UniProt ID
def download_alphafold_prediction(uniprot_id):
    link_pattern = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'
    link = link_pattern.format(uniprot_id)      
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    response = session.get(link) 
    predictedFile = uniprot_id + "_predicted.pdb"
    open("<DESTINATION PATH>"+predictedFile, "wb").write(response.content) 

#Uses alphafold to predict the structure of a FASTA sequence
def predict(fastaf):
    x = "python3 docker/run_docker.py   --fasta_paths=" + fastaf 
    + " --max_template_date=2020-05-14  "
    + "--data_dir=/scratch1/prithvi/geneticdb/  --db_preset=reduced_dbs"
    os.system(x)

#Downloads PDB file from the RCSB database
def download_rcsb_pdb(protein):
    pdbfile = protein + ".pdb"
    pdb_url = "https://files.rcsb.org/download/" + pdbfile

    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter) 
    response = session.get(pdb_url)

    open("<DESTINATION PATH>"+pdbfile, "wb").write(response.content)

#Runs a BLAST search on 2 proteins and gives output in the form of 
# the protein and its percent identity
# the input 'path' is the path directory to a tempfile that stores the output
def blast_search(path,protein):
    command_line = ['blastp','-query',
                "<DESTINATION PATH>"+protein + '.fasta' 
                , '-out', path, '-outfmt',
                "10 sacc pident",'-db','pdbaa']
        #outfmt gives the format of the file (https://www.biostars.org/p/88944/)
    subprocess.call(command_line)

#Filters the results of the BLAST search based on percent identity. 
# It then creates a text file with the filtered matches.
def filt(path,uniprot):
      file = open(path,'r')
      f = file.readlines()
      out_filename = uniprot + "_matches.txt"
      output = open("<DESTINATION PATH>"+out_filename,"w")

      for i in f:
            splitline = i.split(',')
            pdb = splitline[0]
            matchPercent = splitline[1]
            matchPercent = float(matchPercent)
            if (matchPercent >= 70):
                  matchPercent = str(matchPercent)
                  output.write(pdb + "," + matchPercent + "\n")
      output.close()

#Uses pymol to find RMSD between 2 structures
def compare(alphafold,pdbaa):
    af = alphafold + "_predicted.pdb"
    pdb = pdbaa + ".pdb"
    #cmd.load(af) #!!!Gives error?
    cmd.load("<DESTINATION PATH>"+pdb)
    data  = cmd.align(alphafold, pdbaa)
    data = str(data)
    data = data[1:len(data)-1]
    data = data.split(", ")
    rmsd = data[0]
    return rmsd

#calculates avg RMSD
def meanRMSD():
    file = open("rmsd.txt",'r')
    f = file.readlines()
    c = 0
    tot = 0

    for i in f:
        spl = i.split(',')
        rmsd = spl[2]
        tot = tot + rmsd
        c = c + 1
    
    mean = tot / c
    return mean

def sendEmail():
    with smtplib.SMTP_SSL('smtp.gmail.com',465) as connection:
        sender = "example@gmail.com"
        password = "PASSWORD"
        connection.login(sender,password)
        receiver = "example@email.com"
        message = """Subject: Program complete

        The program has been completed on your device.
        """
        connection.sendmail(from_addr=sender,to_addrs=receiver,msg=message)

def main(filepath):
    file = open(filepath,'r')
    f = file.readlines()

    rmsdfile = open("rmsd.txt",'w') #This is the final file with RMSD values


    for i in f:
        uniprot_id = i.strip()
        
        downloadUniprotFasta(uniprot_id)
        download_alphafold_prediction(uniprot_id)
        
        af = uniprot_id + "_predicted.pdb"
        cmd.load("<DESTINATION PATH>"+af) 
        #Loads the alphafold predicted structure into pymol
        
        tempblast = tempfile.NamedTemporaryFile(mode = 'w+') 
        #Temporary file to store the output of the BLAST search

        tempblastname = tempblast.name #Filepath of the temporary file in order to filter results 
        #this tempfile name is fed as the path in blast_search() 
        # in order to be accessed and written in
        blast_search(tempblastname,uniprot_id)
        filt(tempblastname,uniprot_id) 
        tempblast.close()

        pdblist = open("<DESTINATION PATH>"+uniprot_id+"_matches.txt",'r')
        #List of filtered matches from the BLAST search that needs to be downloaded and compared
        
        p = pdblist.readlines()

        for j in p:
              pdbfile = j[0:4]
              spl = j.split(',')
              matchPercent = spl[1].strip()
              download_rcsb_pdb(pdbfile)
              rmsdfile.write(uniprot_id + ',' + pdbfile + ',' 
                             + matchPercent + ',' + compare(uniprot_id,pdbfile)+'\n')
    
    rmsdfile.close()

    print("Program complete. Opening rmsd.txt...")
    os.system("open rmsd.txt")
    sendEmail()

main(testInput)