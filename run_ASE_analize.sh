#!/bin/bash

# === Set path variables ===
CONTAINER_PATH="/Container/path"
SOFTWARE_PATH="/Software/path"

# === Define SIF container commands ===
SIF_GATK="singularity exec -B /data:/data ${CONTAINER_PATH}/gatk_latest.sif"
SIF_ASCAN="singularity exec -B /data:/data ${CONTAINER_PATH}/ascan.sif"
PICARD="${SOFTWARE_PATH}/picard.jar"
GENEIASE="${SOFTWARE_PATH}/geneiase-1.0.1/bin/geneiase"
PHASER_PY="${SOFTWARE_PATH}/phaser-master/phaser/phaser.py"
PHASER_AE_PY="${SOFTWARE_PATH}/phaser-master/phaser/phaser_gene_ae.py"


## Pipeline 01 GATK ASEReadCounter + geneiase
# === Function 1: Add Read Group information ===

function AddReadGroups() {
    echo ">>> Adding Read Groups"
    cd 01data
    while read sample; do
        java -jar $PICARD AddOrReplaceReadGroups \
            I=${sample}.sort.bam \
            O=${sample}.sort.rg.bam \
            RGID=$sample RGPU=unit1 RGLB=lib1 RGPL=illumina sampleM=$sample
    done < sampleid.txt
    cd ..
}

# === Function 2: Generate reference index + Run GATK ASEReadCounter ===
function RunASEReadCounter() {
    echo ">>> Running GATK ASEReadCounter"
    mkdir -p 02count
    cd 00ref
    $SIF_GATK gatk CreateSequenceDictionary -R ref.fa
    java -jar $PICARD CreateSequenceDictionary R=ref.fa O=ref.dict
    samtools faidx ref.fa
    cd ..

    while read sample; do
        $SIF_GATK gatk ASEReadCounter \
            -R 00ref/ref.fa \
            -I 01data/${sample}.sort.bam \
            -V 01data/sample_all_filter.mask.clean.vcf.gz \
            -O 02count/${sample}.output.csv
    done < 01data/sampleid.txt
}

# === Function 3: Run geneiASE analysis ===
function RunGeneiASE() {
    echo ">>> Running geneiASE analysis"
    mkdir -p 03ASEresult
    cd 03ASEresult

    while read sample; do
        awk -F',' 'NR>1{print $1,$2,$2}' OFS='	' ../02count/${sample}.output.csv > ${sample}.output.bed

        bedtools intersect -a ${sample}.output.bed \
                           -b ../00ref/ref.gff3 -wb > ${sample}.gene.bed

        Rscript ../scripts/00ASE_Read_Counter.R \
                ../02count/${sample}.output.csv \
                ${sample}.gene.bed $sample

        $GENEIASE -t static -i ${sample}_GeneiASE_input.tab -b 100

        Rscript ../scripts/00result_table.R \
                ${sample}_GeneiASE_input.tab.static.gene.pval.tab \
                ${sample}_merga.csv ${sample}_result.csv
    done < ../01data/sampleid.txt

    cd ..
}


## Pipeline 02 aScan Pipeline (https://github.com/Federico77z/aScan)
#  Run aScan analysis

function RunAScan() {
    echo ">>> Running aScan"
    mkdir -p 04ascan
    cd 04ascan

    while read sample; do
        $SIF_ASCAN aScan \
            --rna ../01data/${sample}.sort.bam \
            --vcf ../01data/sample_SNP_filter.mask.vcf \
            --gtf ../00ref/ref.gtf \
            -p 36
    done < ../01data/sampleid.txt

    cd ..
}


## Pipeline 03 phASER +DESeq2 Pipeline (https://github.com/secastel/phaser)
# === Function 1: Run phASER analysis ===


function RunPhaser() {
    echo ">>> Running phASER"
    mkdir -p 06phaser
    cd 06phaser

    N=6  # Number of concurrent threads
    i=0

    while read sample; do
        (
            echo ">>> Processing $sample"
            python $PHASER_PY \
                --vcf ../01data/sample_SNP_filter.mask.vcf.gz \
                --bam ../01data/${sample}.sort.bam \
                --paired_end 1 \
                --mapq 20 \
                --baseq 10 \
                --sample sample \
                --threads 12 \
                --o ${sample}_ase > ../logs/${sample}.log 2>&1

            if [ $? -ne 0 ]; then
                echo "[ERROR] Sample ${sample} failed." >> ../logs/error_samples.log
            else
                echo "[OK] Sample ${sample} completed."
            fi
        ) &

        ((i++))
        if (( i % N == 0 )); then wait; fi
    done < ../01data/sampleid.txt
    wait
    cd ..
}

# === Function 2: Parse phASER AE results ===
function RunPhaserAE() {
    echo ">>> Generating gene-level AE tables"
    cd 06phaser

    while read sample; do
        python $PHASER_AE_PY \
            --haplotypic_counts ${sample}_ase.haplotypic_counts.txt \
            --features ../00ref/ref.gene.bed \
            --o ${sample}_gene_ae.txt
    done < ../01data/sampleid.txt

    cd ..
}
function Main() {
    AddReadGroups
    RunASEReadCounter
    RunGeneiASE
    RunAScan
    RunPhaser
    RunPhaserAE
}

Main
