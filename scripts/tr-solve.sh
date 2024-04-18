#!/bin/bash

# Define the directory containing the VCF files
vcf_dir="/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/HPRC/trgt-v0.8.0-9bd9f00/"

# Define the output file
output_file="/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/Git/STRchive_manuscript/data/combined_tr-solve_output.tsv"

# Remove the output file if it exists
rm -f "$output_file"

# Loop through each VCF file in the directory
for vcf_file in "$vcf_dir"/*.vcf.gz; do
    # Get the filename without the directory path
    filename=$(basename "$vcf_file")

    # Run the Python script for the current VCF file and append the output to the output file
    python scripts/tr-solve_HPRC_100.py "$vcf_file" >> "$output_file" 2>&1

done


#awk '$6 ~ /n\(/ {print $5}' data/HPRC.allele_lengths.tsv | sort -u > data/locus_id_with_compound_structure.txt

### in python, done in ipynb:
### f = pd.read_csv('/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/Git/STRchive_manuscript/data/multiple_motif_combined_tr-solve_output.tsv',
### sep='\t', header = None)
#f.iloc[:, 5] = f.iloc[:, 5].str.replace("ID=", "")

# Manually review loci with multiple motifs to differentiate loci compound locus from where there are multiple motifs,
# non specific locus structures (CNG), and interruptions that do not represent variable motifs (info taken from STRchive)

### First, get easy output to manually review
# Step 1: Get unique IDs from column 5
#unique_ids = f.iloc[:, 5].unique()

# Create an empty list to store data for the new DataFrame
#data = []

# Step 2: Iterate over unique IDs and extract motifs
#for unique_id in unique_ids:
#    filtered_df = f[f.iloc[:, 5] == unique_id]
#    motifs_strings = filtered_df[2]

#    unique_motifs = set()
#    for motif_str in motifs_strings:
#        motifs = re.findall(r"'([^']*)'", motif_str)
#        unique_motifs.update(motifs)

    # Append data for the new DataFrame
#    data.append({'Gene': unique_id, 'UniqueMotifCount': len(unique_motifs), 'Motifs': ', '.join(unique_motifs)})

# Create the new DataFrame
# new_df = pd.DataFrame(data)

### Output with manual review notes
### Removed
#SCA7_ATXN7	2	GCC, CAG ***COMPOUND
#DM2_CNBP	3	AC, AGGCAGGCAGGCAGGCAGGCAGGCAGGCAGAC, AGGC ***COMPOUND
#FRDA_FXN	3	GAA, AAG, A *** COMPOUND
#SCA36_NOP56	2	GCCTGG, GCCTGC *** COMPOUND
#CCHS_PHOX2B	2	GCT, GCC ***GCN
#OPMD_PABPN1	2	GCA, GCG ***GCN
#SCA17_TBP	2	CAA, CAG *Previously noted interruption
#SCA8_ATXN8OS	2	CTG, CTA *Previously noted interruption
#NIID_NOTCH2NLC	2	GGA, GCG *Mentioned interruption
#HD_HTT	4	 non canonical: GCC, CCT,    "locus_structure":"(CAG)*CAACAG(CCG)*" known to have interruptions

### Considered 10
#FAME4_YEATS2	3	TTTAT, TTATG, TGTTA YUP
#FAME7_RAPGEF2	2	TTTAT, TTATG YUP
#FTDALS1_C9orf72	3	GCCCCG, ACCGCA, CCCCGGGCCCGC YUP
#DBQD2_XYLT1	7	GCC, CGCGG, GGGCGC, GC, AGG, AGGCGGG, AGGG YUP
#CANVAS_RFC1	17	GAAAG, GGGAA, AAAGG, GGAAAG, GGAA, AAAGGG, AGGAA, A, AAGAG, GAAAA, GAAGA, AAGGG, GGGAAGGAA, GGAAAA, GGAAA, AAAGA, AAGGA
#SCA27B_FGF14	2	AGC, GAA
#SCA31_BEAN1	3	ATAAC, ATAAA, A
#FAME2_STARD7	10	AAATAAATAAAATA, ATAAA, AAAAT, AAATAAAATAACATA, AAAAG, ACAAA, AACAT, ATAAC, AAATA, AACAC
#FAME1_SAMD12	2	ATAAC, AATAA
#SCA4_ZFHX3	3	CCGCCGCCACTGCCA, GCCACT, CCG

