# Epichloe WGS scripts and supplementary files

This is a repository for scripts and files associated with the submitted MS describing assembly and annotation of whole genome sequence for Epichloe typhina and Epichloe clarkii.


# Readme for annotation files
For all these annotation filenames - "Ec" refers to Epichloe clarkii and "Et" to Epichloe typhina

Raw data (tsv) from Interproscan 
Ec_mb_itps.tsv
Et_mb_itps.tsv

Using R we created a new file for downstream analysis.

We retained only information from the following databases 
"CDD", "Gene3D", "PANTHER" , "Pfam", "SUPERFAMILY"
These seemed to be the most useful and complete.

We then reorganised the files so that each gene/protien is on a single line and added gene position (start/end). 
These files are tab deliminted.
genAnn_etmb.txt
genAnn_ecmb.txt


SignalP and effectorP produced the following fasta files.
Ec_signalP.fa
EffectorCandidates_Ec.fasta
EffectorCandidates_Et.fasta
Et_signalP.fa

Effector genes, effector probability and gene position files
EffectorCandidateGenesPos_Ec.txt
EffectorCandidateGenesPos_Et.txt
 
