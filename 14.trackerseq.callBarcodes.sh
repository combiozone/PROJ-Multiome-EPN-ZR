#!/bin/bash

# raw data
rawFQ_dir='/path/folder'
n_cell=5000  # set expected cells
outdir='out_trackerbarcodes'


mkdir -p $outdir


fastqR1=$outdir/TrackerSeq.merged.R1.fastq.gz
fastqR2=$outdir/TrackerSeq.merged.R2.fastq.gz

cat $rawFQ_dir/*_R1_001.fastq.gz > $fastqR1
cat $rawFQ_dir/*_R2_001.fastq.gz > $fastqR2



## 1.Trimming and selecting (maq=20)
$Bbduk in1=$fastqR1 in2=$fastqR2 k=17 literal=GACTCTGGCTCACAAAT ktrim=r out=stdout.fq int=f skipr1 maq=20 |\
$Bbduk in=stdin.fq literal=CTGA k=4 restrictleft=4 ktrim=l out1=$outdir/trim_R1.fastq out2=$outdir/trim_R2.fastq outm1=$outdir/dis_R1.fastq outm2=$outdir/dis_R2.fastq int=t skipr1

## 2. Remove short reads and re-pair
# Remove short reads and re-pair R1 and R2 (for libraries using piggybac)
$Bbduk in=$outdir/trim_R2.fastq out=stdout.fq minlength=37 maxlength=37| \
$Repair in1=stdin.fq in2=$outdir/trim_R1.fastq out1=$outdir/pair_R1.fastq out2=$outdir/pair_R2.fastq repair

## 3. identify correct cell barcodes
$UMI_Tools whitelist --stdin $fastqR1 \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--set-cell-number=${n_cell} \
--log2stderr > $outdir/whitelist.txt

## 4. Add the cell barcodes from step 3 to R2 reads 
$UMI_Tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--stdin $outdir/pair_R1.fastq \
--stdout $outdir/pair_R1_extracted.fastq.gz \
--read2-in $outdir/pair_R2.fastq \
--read2-out=$outdir/dataset_barcode_extracted.fastq \
--filter-cell-barcode \
--error-correct-cell \
--whitelist=$outdir/whitelist.txt




