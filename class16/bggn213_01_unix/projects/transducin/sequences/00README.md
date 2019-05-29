# GalphaI (Transducin-like) mRNA sequence data

File 'data/seqdump.fasta' is from NCBI collaborator Stephen Altschul.

## Date Started:
2015-08-22

## Notes on supplied sequences:
These are mRNA sequence matches to their gene of interest obtained
 using TBLASTN 2.2.32 and the 'nr/nt nucleotide collection' database
 (limit to 100 seqs).

### Search Parameters used:
Program tblastn
Word size       3
Expect value    10
Hitlist size    100
Gapcosts        11,1
Matrix  BLOSUM62
Low Complexity Filter   Yes
Filter string   L;
Genetic Code    1
Window Size     40
Threshold       13
Composition-based stats 2

Database nr/nt
Posted date     Aug 13, 2015 2:35 AM
Number of letters       101,788,032,251
Number of sequences     31,670,722



## Description:
Our collaborator wants us to find 'functional' conserved regions in GalphaI for
targeting mutagenesis experiments to help determine the functional role of his gene.

Here we first:
   - Find PFAM database entry that best matches his supplied mRNA
   - Download corresponding PFAM alignment
   - Analyze sequence conservation per position
   - Produce a summary figure with conserved positions highlighted.



Extract first entry in 'data/seqdump.fasta'
## Q. What line does 1st and 2nd sequence start on?
grep -n "^>" ../data/seqdump.fasta | head -n2

# Answer: 1 and 45

## Extract first sequence to a new file called 'data/single_seq.fa'
head -n44 ../data/seqdump.fasta > ../data/single_seq.fa



## Search PDB to find if there are any structures of this protein available?
Used sequence to blastx against PDB data base with default parameters
See: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx

## Found excellent hits to G-Alpha-I1 (e.g. 1GIT etc.)

# Used 1GIT to search PFAM and found entry PF00503, see:
http://pfam.xfam.org/structure/1GIT#tabview=tab3
http://pfam.xfam.org//family/PF00503#tabview=tab3


## Conservation/Alignment analysis
Download seed alignment (to 'data/PF00503_seed.txt') from PFAM (version 28.0)
and analyzed with 'scripts/01_conservation.r'

See conservation profile per position in 'results/conservation.pdf'

