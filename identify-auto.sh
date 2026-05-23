#!/bin/bash
set -euo pipefail

PREFIX="GCA_013401275.1_ASM1340127v1"
GENOME="${PREFIX}_genomic.fna"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/401/275/${PREFIX}"

REF_PRO="/path/PHASIANIDAE/ZW.pro.fa"

# 1. Download genome
wget -c "${BASE}/${PREFIX}_genomic.fna.gz"

# 2. Decompress, keeping original gz file
if [ ! -s "${GENOME}" ]; then
    gzip -dk "${GENOME}.gz"
fi

# 3. Build nucleotide BLAST database for target genome
makeblastdb -in "${GENOME}" -dbtype nucl -out target_genome

# 4. Search reference Z/W proteins against target genome
tblastn \
    -query "${REF_PRO}" \
    -db target_genome \
    -out out_ZW_vs_target.blast \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -evalue 1e-5 \
    -num_threads 5

# 5. Filter candidate Z/W scaffolds
# Criteria:
# pident >= 95
# alignment length >= 100 aa
# evalue <= 1e-10
# bitscore >= 50
awk '$3 >= 95 && $4 >= 100 && $11 <= 1e-10 && $12 >= 50 {print $2}' \
    out_ZW_vs_target.blast | sort -u > ZW_candidate.name.txt

# 6. Get all scaffold names
grep '^>' "${GENOME}" | awk '{print $1}' | sed 's/>//g' | sort -u > scaffold.name

# 7. Autosomal candidate scaffolds = all scaffolds minus Z/W candidates
comm -23 scaffold.name <(sort -u ZW_candidate.name.txt) > Auto.list.txt

# 8. Extract sequences
seqtk subseq "${GENOME}" Auto.list.txt > auto.fna
seqtk subseq "${GENOME}" ZW_candidate.name.txt > ZW1.fna

# Build nucleotide BLAST database from first-round candidate Z/W scaffolds
makeblastdb -in ZW1.fna -dbtype nucl -out ZW1

# Second-round BLAST: reference Z/W proteins against candidate target scaffolds
tblastn \
    -query /path/PHASIANIDAE/ZW.pro.fa \
    -db ZW1 \
    -out out2.blast \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 1e-5 \
    -num_threads 5

# Extract confirmed candidate Z/W scaffold names
awk '$3 > 95.00 {print $2}' out2.blast | sort -u > ZW2.name.txt

# Extract all scaffold names from target genome
grep '^>' GCA_013401275.1_ASM1340127v1_genomic.fna \
    | awk '{print $1}' \
    | sed 's/>//g' \
    | sort -u > scaffold.name

# Get autosomal candidate scaffold names: all scaffolds minus Z/W candidates
comm -23 scaffold.name <(sort -u ZW2.name.txt) > Auto.list.txt

# Extract autosomal and Z/W candidate fasta sequences
seqtk subseq GCA_013401275.1_ASM1340127v1_genomic.fna Auto.list.txt > auto.fna
seqtk subseq GCA_013401275.1_ASM1340127v1_genomic.fna ZW2.name.txt > ZW.fna

# Build BWA index for autosomal candidate genome
bwa-mem2 index -p PAN auto.fna
samtools faidx auto.fna

# Check sequence statistics
seqkit stats auto.fna
seqkit stats ZW.fna
