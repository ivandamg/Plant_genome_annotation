# Plant_genome_annotation
Pipeline for plant genome annotation

## Steps

### 1. purge haplotigs assembly.
https://bitbucket.org/mroachawri/purge_haplotigs/src/master/
## 0. Activate

    conda create -n purge_haplotigs -c conda-forge -c bioconda samtools bedtools r-base r-ggplot2 minimap2
    conda activate purge_haplotigs

## 1. Map pacbio subreads

    minimap2 -t 4 -ax map-pb hap1.fasta ../01_Hifi/Combined_clean.fq --secondary=no | samtools sort -m 1G -o hap1_aligned.bam -T hap1_tmp.ali

## 2. Purge haplotigs

    purge_haplotigs  hist  -b aligned.bam  -g genome.fasta  [ -t threads ]

    sbatch --partition=pibu_el8 --job-name=purge2hap1 --time=0-10:00:00 --mem-per-cpu=60G --ntasks=12 --cpus-per-task=1 --output=purge2.out --error=purge2.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/20_GenomeAssembly/04_LFassembly; conda activate purge_haplotigs; purge_haplotigs  hist  -b hap1_aligned.bam  -g hap1.fasta"


## 3. Analyse coverage

    purge_haplotigs  cov  -i aligned.bam.genecov  -l <integer>  -m <integer>  -h <integer> [-o coverage_stats.csv -j 80  -s 80 ]
    sbatch --partition=pibu_el8 --job-name=purge2hap1 --time=0-10:00:00 --mem-per-cpu=60G --ntasks=12 --cpus-per-task=1 --output=purge2.out --error=purge2.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/20_GenomeAssembly/04_LFassembly; conda activate purge_haplotigs; purge_haplotigs  hist  -b hap1_aligned.bam  -g hap1.fasta"

## 4. Run purging pipeline

    purge_haplotigs  purge  -g genome.fasta  -c coverage_stats.csv



### 2. Repeat mask the assembly

### 3. Use RNAseq evidence for genemark prediction

### 4. Braker pipeline
