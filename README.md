# Comparative-genomics

# Begin by downloading your software
[Set up bioconda if you haven't already](https://bioconda.github.io/)<br />
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
Install blast, mafft, trimal and hmmer using conda
```
conda install -c bioconda blast
conda install -c bioconda mafft
conda install -c bioconda trimal
conda install -c bioconda hmmer
```
 
Download files for Exercise 1 [here](https://www.dropbox.com/sh/3944t4l8xhq7rp8/AABF2IrgAmThs4nIvaIDETUta?dl=0)

Download Orthofinder output [here](https://www.dropbox.com/s/jpvp84gvh4g8fq4/Exercise2_Heterocephalus_Dataset1.zip?dl=0) and [here](https://www.dropbox.com/s/0bsu89iadyu3zig/Exercise2_Heterocephalus_Dataset2.zip?dl=0) They are about 1GB each, so don't worry if downloading them isn't feasible

# Exercise 1 - Proteins crucial for animal multicellularity

(Example freely adapted from Parfrey, L. W., & Lahr, D. J. 2013. Multicellularity arose several times in the evolution of eukaryotes. (Response to DOI 10.1002/bies. 201100187). Bioessays, 35: 339-347 https://doi.org/10.1002/bies.201200143)

Animals are descended from a last common ancestor that possessed complex multicellularity: a lifestyle characterized by spatial division of labour between multiple differented cell types, arising from a genetically regulated program of embryonic development. These features were made possible by sets of adhesion proteins (such as cadherins and - and -catenins), secretion of an anchoring extracellular membrane (including collagen), cell-cell signaling pathways (such as the TGF-b, Hippo or Wnt pathways), and transcription factors (such as NF-kappaB or homeodomain proteins) playing a regulatory role in all of these processes. 

In humans, alpha- and beta-catenin are crucial to the function of adhaerens junctions, which connect adjacent cells; -catenin is believed to anchor actin to catenin. The human alpha-catenin is a paralog of a protein called vinculin, which has a similar actin-anchoring function in focal adhesion junctions. The amoebozoan Dictyostelium discoideum forms multicellular fruiting bodies following cell aggregation – a form of multicellularity distinct from the clonal, embryonic multicellularity found in animals. However, Dictyostelium possesses two homologs of vinculin/alpha-catenin, one of which similarly appears to play a part in the formation of its multicellular structures. This led to a suggestion (Dickinson et al., 2012, BioEssays) that the last common ancestor of animals and amoebozoans possessed a multicellular lifestage. This hypothesis would receive support if the two paralogs in Dictyostelium, and those in animals, were descended from two paralogs already present in that shared ancestor. You will reconstruct the phylogeny of vinculin-like proteins to test this hypothesis.

•	You will begin by searching a paneukaryotic database with a human vinculin query. To save time in this exercise, only a subset of EukProt taxa have been included in this dataset (EukProt_v3mini), and we will be running a blast search against it rather than performing a much more computationally intensive Orthofinder search (see below). I’ve provided a list of information about the taxonomic affiliation of each species included, and the origins of its data, for subsequent ease of interpretation. Note that this database includes eukaryotes only – in most cases, you will also want to find out whether your protein of interest is present in prokaryotes and/or viruses as well.

Following the blast search, you will retrieve the subject sequences, align them using mafft (with the accurate l-ins-i option), trim them using trimal with the gappyout option, and construct a phylogeny using IQ-TREE.

Make a blastable database from an infile containg amino acid sequences (for nucleotide sequences, use -dbtype nucl instead)
```
makeblastdb -in EukProt_v3mini.fasta -dbtype prot -out EukProt_v3miniDB
```
search the EukProt_v3mini database using Homo sapiens Vinculin as a query, with tabular output and a cutoff value of 0.001. See blast parameters [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
```
blastp -query Vinculinquery.fasta -db EukProt_v3miniDB -outfmt 6 -evalue 1e-3 -out VinculininEukaryotes.fasta
```
retrieve the sequences corresponding to the hits
```
cut -f 2 VinculininEukaryotes.fasta > IDs.txt ; seqkit grep -n -f IDs.txt EukProt_v3mini.fasta > Eukaryotes_Vinculin.fasta
```
remove duplicate entries based on sequence (i.e., identical Homo sapiens sequences will be trimmed, even though they have different header). Leaving off the 's' flag will remove duplicate sequences only if they also have identical fasta headers, and remove sequences shorter than 150 amino acids
```
seqkit seq -m 150 Eukaryotes_Vinculin.fasta | seqkit rmdup -s > Eukaryotes_Vinculin_reduced.fasta
```
align the sequences in Eukaryotes_Vinculin_reduced.fasta using MAFFT with the accurate L-INS-I option
```
linsi Eukaryotes_Vinculin_reduced.fasta > Eukaryotes_Vinculin_reduced.linsi.fasta
```
alternatively, use [the mafft web server](https://mafft.cbrc.jp/alignment/server/)

Visualise the alignment using a viewer of your choice(Jalview, AliView, SeaView, Geneious...), or [Alignment Viewer](https://alignmentviewer.org/) in a browser window

trim the alignment using trimal, with the gappyout option, and phylip outformat readable by paml (and some other programs such as FastTree)
```
trimal -in Eukaryotes_Vinculin_reduced.linsi.fasta -out Eukaryotes_Vinculin_reduced.linsi.gappyout.phy -gappyout -phylip_paml
```
Reconstruct the phylogeny [IQ-TREE's webserver](http://iqtree.cibiv.univie.ac.at/)

Use the default options, and note that [ultrafast bootstrap support values are interpreted slightly differently to regular nonparametric bootstrap support values](https://academic.oup.com/mbe/article/30/5/1188/997508?login=false)

Check the log file output in a text viewer, and the maximum likelihood tree with support values in [iTOL](https://itol.embl.de/)

•	If you have extra time to spare/waste at some point, try the same exercise using the Drosophila extracellular matrix heparin sulfate proteoglycan Terribly Reduced Optic Lobes (Trol) as a query, as an illustration of an unpleasant situation that would normally require much more curation than the p53 example above (and is unlikely to give usable results). At minimum you would need to start by making the search more stringent and/or changing the maximum number of hits.

# Exercise 2 - insight into poikilothermy and blindness in the naked mole rat 
(example freely adapted from Kim, E., Fang, X., Fushan, A. et al. Genome sequencing reveals insights into physiology and longevity of the naked mole rat. Nature 479, 223–227 (2011). https://doi.org/10.1038/nature10533)

<img width="287" alt="image" src="https://user-images.githubusercontent.com/94364625/180106445-a9d50d7c-004d-47de-99e3-3e5a16887aed.png"><figcaption>Angry female naked mole rat. Credit: Buffenstein/Barshop Institute/UTHSCSA, Creative Commons CCBY 2.0</figcaption>


The naked mole rat Heterocephalus glaber (Rodentia, Castorimorpha) is an unusual mammal in several respects: it burrows deep underground, and has poor eyesight as a result; it is eusocial, with reproductive and non-reproductive castes; it is poikilothermic; it is hairless; its skin is insensitive to pain; and it is remarkably resistant to oxidative stress-related damage and cancer.

We will be looking at data from the naked mole rat's predicted proteome to understand the genetic basis for some of the traits described above. We will use Orthofinder to infer orthogroups in the proteome of the naked mole rat and some related mammals: the kangaroo rat Dipodomys spectabilis (Rodentia, Castorimorpha), the brown rat Rattus norvegicus (Rodentia, Muridae), the house mouse Mus musculus (Rodentia, Muridae), the human Homo sapiens (primates, Hominidae), and the platypus Ornithorhynchus anatinus (Monotremata, Ornithorhynchidae). Orthofinder2 is an orthogroup inference program that begins with sequence similarity searches between each pair of sequences, clusters the results, and provides phylogenetic reconstructions and information on gene duplication events for each orthogroup.

You will not be running Orthofinder yourselves; I have provided its output, but here is how to install and run it on all of the fasta sequences contained in the directory Heterocephalus_Dataset1, using diamond as the default search tool, or alternatively using BLAST
```
conda install -c bioconda orthofinder
orthofinder -f Heterocephalus_Dataset1
orthofinder -f Heterocephalus_Dataset1 -S blast
```

•	Have a look at the Orthofinder output. The folder names are fairly self-explanatory, and the individual files provide information on the numbers of genes assigned to orthogroups overall and in each species, gene duplication events,  orthologs determined based on the single gene phylogenies. Notice that Orthofinder distinguishes between ‘Gene trees’ and ‘Resolved gene trees’. The files may be opened in a text editing program, or (in the case of tab-delimited .tsv files) in a spreadsheet viewer such as Excel. The trees may be viewed by uploading them to https://itol.embl.de/upload.cgi 
Note that I generated this output using diamond as the tool for the initial all-against-all sequence similarity searches. This is much faster, but less accurate.

•	The mitochondrial proton carrier UCP1 is a major regulator of non-shivering thermoregulation in mammals. It is known to be regulated by fatty acids and purines, which interact with a regulatory site containing conserved G263 and P264 residues. Its proton pumping activity is dependent on a conserved H146 residue; mutating this residue to glutamine sharply reduces proton pumping capacity. 
The mouse UCP1 sequence is NP_033489.1. Identify the orthogroup that includes this protein in Orthogroups/Orthogroups.tsv, and retrieve the corresponding sequences from Orthogroup_Sequences. Align the sequences, and visualize the alignment to identify any notable features of naked mole rat UCP1. (Does looking at this alignment raise any questions about other parts of Orthofinder’s output?)

•	Humans have four visual opsin genes: rhodopsin (RHO, NP_000530.1), long-wave sensitive opsin (OPN1LW, NP_064445.2), medium-wave sensitive opsin (OPN1MW, NP_000504.1), and short-wave sensitive opsin (OPN1SW, NP_001699.1). Look for orthologs of these proteins in the naked mole rat proteome. (Optional: look at this table from the authors of the 2011 naked mole rat genome publication, and see if you agree with [their conclusions with regard to visual opsin loss in the naked mole rat](https://www.nature.com/articles/nature10533/tables/2)

You happen to have some output from an Orthofinder run that compared the mammals listed above with some additional animals (the platyhelminth Schistosoma mansoni, the sponge Amphimedon queenslandica, the ctenophore Mnemiopsis leidyi, the filasterean Capsaspora (a protist closely related to animals), the apicomplexan Plasmodium falciparum, and the plant Arabidopsis thaliana). This output should illustrate some potential sources of error that you may come across when using Orthofinder.
•	What is going on in OG0025093?
•	What is going on in OG0011927?
•	Consider how Orthofinder handles proteins involved in the folate biosynthesis pathway shown in Figure A of [this poster](https://finlaymagui.re/assets/posters/sgm_poster.pdf)

# Exercise 3 - Membrane fusion proteins in Archaea
(example freely adapted from Moi, D., Nishio, S., Li, X. et al. Discovery of archaeal fusexins homologous to eukaryotic HAP2/GCS1 gamete fusion proteins. Nat Commun 13, 3880 (2022). https://doi.org/10.1038/s41467-022-31564-1 )

<img width="202" alt="image" src="https://user-images.githubusercontent.com/94364625/180107201-3431d2f2-0e9e-44d9-bfb6-7b0e872ba2ed.png"><figcaption>Fédry et al. 2017 Cell graphical abstract</figcaption>


Based on the widespread presence of meiotic genes among eukaryotes, the last eukaryotic common ancestor is believed to have been capable of sexual reproduction. Its genome is also believed to have encoded a protein, HAP2 (Hapless 2, also known as Generative cell specific protein 1), that is involved in gamete fusion in a wide variety of extant eukaryotes and bears structural similarity to viral class II fusion proteins.

You are interested in finding out whether any archaea possess HAP2 homologs. Large quantities of sequence data from isolated archaea and metagenomic data are now available in repositories such as GenBank® and Uniclust, but in the interests of time and your hard drive space, we are going to limit this search to a very small subset of archaeal sequences derived from metagenome data.

•	Begin by searching the archaeal protein sequences with a set of curated eukaryotic HAP2 sequences, HAP2.fasta, as a query file.

Make a blastable database from the archaeal MAGs infile containg amino acid sequences
```
makeblastdb -in ArchaealMAGs_protein.fasta -dbtype prot -out ArchaealMAGsDB
```
Perform a blast search against the archaeal sequences using the eukaryotic HAP2 sequences as queries
```
blastp -query HAP2.fasta -db ArchaealMAGsDB -out HAP2inArchaealMAGs.fasta
```
If you don’t get any significant hits, it is still possible that the archaeal protein set actually does encode proteins homologous to HAP2, but so divergent in sequence that sequence similarity methods are not sensitive enough to detect them.

To perform a more sensitive search, you will construct a profile Hidden Markov Model (HMM). Briefly, this involves using a set of known HAP2 proteins to generate a statistical model of amino acid frequencies (and indel probabilities) at each site. This profile HMM can then be used to search a database for sequences that fit the profile.

•	Construct a multiple sequence alignment from HAP2.fasta, and use HMMER to build a profile HMM from it. Search the set of archaeal proteins with the profile HMM

Make an alignment using MAFFT with the accurate L-INS-I option from your HAP2 sequences
```
linsi HAP2.fasta > HAP2.linsi.fasta
```
Build a profile Hidden Markov Model from your alignment
```
hmmbuild HAP2.hmm HAP2.linsi.fasta
```
Search a fasta file of sequences for sequences that match an hmm profile
```
hmmsearch HAP2.hmm ArchaealMAGs_protein.fasta > HAP2inArchaealMAGs.out
```

•	We’ll use a global e-value of 10-3 and an alignment length of 100 amino acids as our cutoffs. If you get any significant hits, use InterProScan (https://www.ebi.ac.uk/interpro/result/InterProScan ) to predict their sequence features and protein family memberships.
•	Bonus: download ChimeraX, and predict the structures of your hits using Alphafold2 via CollabFold (ColabFold: Making protein folding accessible to all. Nature Methods (2022)):
Tools/Structure Prediction/Alphafold ; paste your sequences, separated by commas, and click ‘Predict’. This requires signing in with a Google account

This type of search is useful to try to identify very poorly conserved sequences, including transmembrane proteins such as SNAREs or mitochondrial membrane proteins, especially in genomes of distantly related organisms. It is also useful as a more sensitive way to try to rule out the presence of a protein in a given genome. 
More sophisticated versions use: 
-	structurally informed multiple sequence alignments; 
-	HMM searches against databases that are also composed of HMMs representing clusters of sequences; 
-	iterative HMM searches that take sequences identified in initial searches and incorporate them into new HMM profiles for subsequent searches.
Here we have searched against amino acid data, but HMMer can search against nucleotide data as well. As you can imagine, the output of HMM searches is dependent on taxon sampling. Domain prediction software such as PfamScan works by comparing query sequences against sets of profile HMMs for different domains.
