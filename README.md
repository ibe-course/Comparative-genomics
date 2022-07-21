# Comparative-genomics

# Begin by downloading your software
#Set up bioconda if you haven't already - https://bioconda.github.io/<br />
conda config --add channels defaults<br />
conda config --add channels bioconda<br />
conda config --add channels conda-forge<br />
conda config --set channel_priority strict<br />

#Install blast, mafft, trimal and hmmer using conda
conda install -c bioconda blast<br />
conda install -c bioconda mafft<br />
conda install -c bioconda trimal<br />
conda install -c bioconda hmmer<br />

# Exercise 1

#Make a blastable database from an infile containg amino acid sequences (for nucleotide sequences, use -dbtype nucl instead)
makeblastdb -in EukProt_v3mini.fasta -dbtype prot -out EukProt_v3miniDB

#search the EukProt_v3mini database using Homo sapiens Vinculin as a query, with tabular output and a cutoff value of 0.001. See blast parameters here: https://www.ncbi.nlm.nih.gov/books/NBK279684/
blastp -query Vinculinquery.fasta -db EukProt_v3miniDB -outfmt 6 -evalue 1e-3 -out VinculininEukaryotes.fasta

#retrieve the sequences corresponding to the hits
cut -f 2 VinculininEukaryotes.fasta > IDs.txt ; seqkit grep -n -f IDs.txt EukProt_v3mini.fasta > Eukaryotes_Vinculin.fasta

#remove duplicate entries based on sequence (i.e., identical Homo sapiens sequences will be trimmed, even though they have different header). Leaving off the 's' flag will remove duplicate sequences only if they also have identical fasta headers, and remove sequences shorter than 150 amino acids
seqkit seq -m 150 Eukaryotes_Vinculin.fasta | seqkit rmdup -s > Eukaryotes_Vinculin_reduced.fasta

#align the sequences in Eukaryotes_Vinculin_reduced.fasta using MAFFT with the accurate L-INS-I option
linsi Eukaryotes_Vinculin_reduced.fasta > Eukaryotes_Vinculin_reduced.linsi.fasta
#alternatively, use the mafft web server: https://mafft.cbrc.jp/alignment/server/ 

#Visualise the alignment using a viewer of your choice(Jalview, AliView, SeaView, Geneious...), or Alignment Viewer in a browser window
#https://alignmentviewer.org/

#trim the alignment using trimal, with the gappyout option, and phylip outformat readable by paml (and some other programs such as FastTree)
trimal -in Eukaryotes_Vinculin_reduced.linsi.fasta -out Eukaryotes_Vinculin_reduced.linsi.gappyout.phy -gappyout -phylip_paml
