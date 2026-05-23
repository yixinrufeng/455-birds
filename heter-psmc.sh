#!/bin/bash
prefetch SRR952717
fastq-dump ./SRR952717/SRR952717.sra --split-3 --gzip --defline-qual '+' -O ./
repair.sh in=SRR952717_1.fastq.gz in2=SRR952717_2.fastq.gz out=SRR952717_1.re.fastq.gz out2=SRR952717_2.re.fastq.gz overwrite=f
trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 --clip_R1 5 --clip_R2 5 --fastqc --gzip --output_dir ./ --stringency 1 -e 0.1 ./SRR952717_1.re.fastq.gz ./SRR952717_2.re.fastq.gz
bwa-mem2 mem -t 6 -M -R '@RG\tID:SAMN02318032\tSM:SAMN02318032\tLB:WGS\tPL:Illumina' /path/ACA ./SRR952717_1.re_val_1.fq.gz ./SRR952717_2.re_val_2.fq.gz > ./SRR952717_bwa.sam
gatk --java-options '-Xmx70g' SortSam -I ./SRR952717_bwa.sam -O ./SRR952717.sorted.bam -SO coordinate --CREATE_INDEX true
rm ./SRR952717_bwa.sam
gatk --java-options '-Xmx70g' MarkDuplicates -I ./SRR952717.sorted.bam -O ./SRR952717.dedup.bam -M ./SRR952717.dedup.metrics.txt -REMOVE_DUPLICATES false --TAGGING_POLICY All -MAX_FILE_HANDLES 1000
rm SRR952717.sorted.*
mv SRR952717.dedup.bam SAMN02318032.bam
samtools flagstat -@ 6 ./SAMN02318032.bam > ./SAMN02318032.stat.txt
samtools depth -a ./SAMN02318032.bam > ./SAMN02318032.c.txt
cat ./SAMN02318032.c.txt | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ./SAMN02318032.txt
angsd -i SAMN02318032.bam -anc /path/auto.fna -ref /path/auto.fna -out SAMN02318032 -dosaf 1 -gl 2 -C 50 -minQ 20 -minmapq 30 -P 6
realSFS SAMN02318032.saf.idx > SAMN02318032.ml
wc -l SAMN02318032.c.txt >hangshu.txt
awk '$3 == 0' SAMN02318032.c.txt | wc -l >>hangshu.txt
awk 'NR==1 {a=$1} NR==2 {b=$1; print b/a}' hangshu.txt >per.txt
samtools mpileup -C50 -uf /path/auto.fna ./SAMN02318032.bam | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 32 | gzip > Achl.fq.gz
fq2psmcfa -q20 Achl.fq.gz > Achl.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o Achl.psmc Achl.psmcfa
psmc_plot.pl -u 0.0000000016988 -g 1.890638142 -R Achl Achl.psmc
