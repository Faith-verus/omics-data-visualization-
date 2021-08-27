# Microbial-Genomic-Data-Analysis-
import Bio
import pylab
import urllib
import pandas as pd
import nglview as nv
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from collections import Counter
from Bio.Data import CodonTable
from Bio import SeqIO, SearchIO
from Bio.PDB import PDBParser,MMCIFParser
from Bio.SeqUtils import GC,molecular_weight

## Download genome sequence & protein sequence https://www.ncbi.nlm.nih.gov/assembly/GCF_015602625.1/

#Import sequence file 
from Bio import SeqIO

for record in SeqIO.parse("GCF_015602625.1_ASM1560262v1_genomic.fasta", "fasta"):
    print(record.id)

#Load sequence as a variable 

BSseq = SeqIO.read("GCF_015602625.1_ASM1560262v1_genomic.fasta","fasta")
BSseq

#ensure the file is biological 

type(BSseq)

#Retrive sequence ID
BSseq.id


for record in SeqIO.parse("GCF_015602625.1_ASM1560262v1_genomic.fasta", "fasta"):
    print(record)

#Store sequence for analysis 
BSseqfile = record.seq
BSseqfile

# length of sequence
len(BSseqfile)

#molecular weight
molecular_weight(BSseqfile)

#Calculate GC content as a percentage 
GC(BSseqfile)

BSprotein_seq = BSseqfile.translate()

#Calculate length of protein sequence 
len(BSprotein_seq)

# Listing the 10  most common amino acids
C_amino = Counter(BSprotein_seq)
C_amino.most_common(10)



# visualize all 20 amino acid occurrences in the form of a histogram

del C_amino["*"]
pylab.bar(C_amino.keys(), C_amino.values())
pylab.title("Protein Sequence Frequency")
pylab.xlabel("Amino Acid")
pylab.ylabel("Frequency")

#Split the protein sequence everytime there is a stop codon 

BSprotein_list = [str(i)for i in BSprotein_seq.split("*")]
BSprotein_list[:20]

# convert sequences to dataframe
large_BSproteins = [x for x in BSprotein_list if len(x)>10]
df = pd.DataFrame({"protein_seq": large_BSproteins})
df

# Add a column with sequence lengths
df["length"] = df["protein_seq"].apply(len)
df.head()

# sort sequence data

df.sort_values(by = ["length"], ascending=False)[:10]

# let's take a single protein from the table

largest_protein = df.nlargest(1,"length")
single_protein = largest_protein.iloc[0,0]
single_protein

# write to a file

with open("single_protein.fasta", "w") as file:
    file.write(">largest protein \n"+single_protein)

#Basic Local Alignment using NCBI blast tool 
# Read single_seq.fasta

readF = SeqIO.read("single_protein.fasta","fasta")
readF.seq

%%time 

# based on the server load this query might take 2-3 minutes to run

result_handle = NCBIWWW.qblast("blastp","pdb", readF.seq)
blast_qresult = SearchIO.read(result_handle,"blast-xml")

print(blast_qresult[0:5])

#fetch the id, description, evalue, bitscore & alignment of first hit

seqid = blast_qresult[0]

details = seqid[0]

print(f"\
Sequence ID:{seqid.id}\n\
description:{seqid.description}\n\
E value:    {details.evalue} \n\
Bit Score:  {details.bitscore}\n\
")


#id of protein 
seqid.id

# split seqid - remove all unecessary infor except for ID
seqid.id.split("|")[1]

#use ID to retrive file using URL 
# link format https://files.rcsb.org/download/6MFZ.pdb

urllib.request.urlretrieve("https://files.rcsb.org/download/6MFZ.pdb", "6MFZ.pdb")


parser = PDBParser()
structure = parser.get_structure("6MFZ","6MFZ.pdb")
structure

for chain in structure[0]:
    print(f"chain ID: {chain.id}")

view = nv.show_biopython(structure)
view

#creates an image of the protein and you can save it on your desktop 
view.render_image()

#GUI
nv.show_biopython(structure,gui=True)


