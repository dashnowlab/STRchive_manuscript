#### location of TRGT vcfs
# /scratch/ucgd/lustre-work/quinlan/data-shared/datasets/HPRC/trgt-v0.8.0-9bd9f00
#/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/HPRC/trgt-v0.8.0-9bd9f00/HG00099.vcf.gz
### where I am working
# /uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/Git/STRchive_manuscript


#/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/HPRC/trgt-v0.8.0-9bd9f00/HG00099.vcf.gz

grep -E ',' 1099output.txt

awk '$2 ~ /,/' 1099output.txt


awk '$2 ~ /,/' 1099output.txt

awk '$6 ~ /n\(/' data/HPRC.allele_lengths.tsv

awk '$6 ~ /n\(/ {print $5}' data/HPRC.allele_lengths.tsv | sort -u


# Compound locus structure in HPRC.allele_lengths.tsv
CJD_PRNP
DM2_CNBP
FAME1_SAMD12
FAME2_STARD7
FAME3_MARCHF6
FAME4_YEATS2
FAME6_TNRC6A
FAME7_RAPGEF2
FRDA_FXN
HD_HTT
RCPS_EIF4A3
SCA31_BEAN1
SCA36_NOP56
SCA37_DAB1
SCA7_ATXN7
SCA8_ATXN8OS


grep -E ',' data/combined_output.txt


awk '$3 ~ /,/' data/combined_output.txt > comma_combined_output.tsv

awk '{print $6}' comma_combined_output.tsv