#!/bin/bash
set -euo pipefail

PREFIX="GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/700/215/${PREFIX}"

# 1. Download files
wget -c "${BASE}/${PREFIX}_genomic.fna.gz"
wget -c "${BASE}/${PREFIX}_genomic.gtf.gz"
wget -c "${BASE}/${PREFIX}_assembly_report.txt"

# Optional: only needed if you want to use official NCBI protein sequences
wget -c "${BASE}/${PREFIX}_protein.faa.gz"

# 2. Uncompress
gunzip -f "${PREFIX}_genomic.fna.gz"
gunzip -f "${PREFIX}_genomic.gtf.gz"
gunzip -f "${PREFIX}_protein.faa.gz"

GENOME="${PREFIX}_genomic.fna"
GTF="${PREFIX}_genomic.gtf"
REPORT="${PREFIX}_assembly_report.txt"

# 3. Extract RefSeq accession IDs for Z and W chromosomes
awk -F '\t' '
    $0 !~ /^#/ && $2 == "assembled-molecule" && ($3 == "Z" || $3 == "W") {
        print $7
    }
' "$REPORT" > ZWname.txt

# 4. Extract autosomal chromosome IDs only
# This keeps assembled chromosomes with numeric names, e.g. 1, 2, 3...
awk -F '\t' '
    $0 !~ /^#/ && $2 == "assembled-molecule" && $3 ~ /^[0-9]+$/ {
        print $7
    }
' "$REPORT" > auto.txt

# 5. Extract autosomal fasta
seqtk subseq "$GENOME" auto.txt > auto.fna

# 6. Extract Z/W GTF annotations by exact matching of GTF column 1
awk 'NR==FNR{a[$1]; next} $1 in a' ZWname.txt "$GTF" > ZW.gtf

# 7. Translate Z/W coding sequences to protein sequences
gffread ZW.gtf -g "$GENOME" -y ZW.pro.fa

# 8. Build BLAST protein database for Z/W proteins
makeblastdb -in ZW.pro.fa -dbtype prot -out G.ZW

# 9. Build BWA index and fasta index for autosomes
bwa-mem2 index -p PHA auto.fna
samtools faidx auto.fna

# 10. Basic checks
echo "Number of Z/W sequences:"
wc -l ZWname.txt

echo "Number of autosomal sequences:"
wc -l auto.txt

echo "Number of autosomal fasta records:"
grep -c "^>" auto.fna

echo "Number of Z/W protein records:"
grep -c "^>" ZW.pro.fa
