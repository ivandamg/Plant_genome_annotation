# Plant_genome_annotation
Pipeline for plant genome annotation



## Soft Masking with TE tools
https://github.com/Dfam-consortium/TETools
Use Docker container in windows 


1. copy files into container
2. BuildDatabase -name genome1 genome1.fa
3. RepeatModeler -database genome1 [-LTRStruct] [-threads 16]
3b. continue interrupted modeler:
4. RepeatClassifier
5. RepeatMasker genome1.fa [-lib library.fa] [-pa 8]

runcoseg.pl -d -m 50 -c ALU.cons -s ALU.seqs -i ALU.ins









### 1. purge haplotigs assembly.
https://bitbucket.org/mroachawri/purge_haplotigs/src/master/
## 0. Activate

    conda create -n purge_haplotigs -c conda-forge -c bioconda samtools bedtools r-base r-ggplot2 minimap2
    conda activate purge_haplotigs

## 1. Map pacbio subreads

    minimap2 -t 4 -ax map-pb hap1.fasta ../01_Hifi/Combined_clean.fq --secondary=no | samtools sort -m 1G -o hap1_aligned.bam -T hap1_tmp.ali
    sbatch --partition=pibu_el8 --job-name=purgeVGPhap1 --time=0-10:00:00 --mem-per-cpu=60G --ntasks=12 --cpus-per-task=1 --output=purge2.out --error=purge2.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/20_GenomeAssembly/03_VGP/01_hap1;conda activate purge_haplotigs;minimap2 -t 12 -ax map-pb VGPhap1_only20.fasta ../../01_Hifi/Combined_clean.fq --secondary=no | samtools sort -m 1G -o VGPhap1_aligned.bam -T VGPhap1_tmp.ali"

## 2. Purge haplotigs

    purge_haplotigs  hist  -b aligned.bam  -g genome.fasta  [ -t threads ]

    sbatch --partition=pibu_el8 --job-name=purge2hap1 --time=0-10:00:00 --mem-per-cpu=60G --ntasks=12 --cpus-per-task=1 --output=purge2.out --error=purge2.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/20_GenomeAssembly/03_VGP/01_hap1; conda activate purge_haplotigs; purge_haplotigs  hist  -b VGPhap1_aligned.bam  -g VGPhap1_only20.fasta"


## 3. Analyse coverage

    purge_haplotigs  cov  -i aligned.bam.genecov  -l <integer>  -m <integer>  -h <integer> [-o coverage_stats.csv -j 80  -s 80 ]
    sbatch --partition=pibu_el8 --job-name=purge2hap1 --time=0-10:00:00 --mem-per-cpu=60G --ntasks=12 --cpus-per-task=1 --output=purge2.out --error=purge2.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/20_GenomeAssembly/04_LFassembly; conda activate purge_haplotigs; purge_haplotigs  hist  -b hap1_aligned.bam  -g hap1.fasta"

## 4. Run purging pipeline

    purge_haplotigs  purge  -g genome.fasta  -c coverage_stats.csv



### 2. Repeat mask the assembly

### 3. Use RNAseq evidence for genemark prediction

### 4. Braker pipeline
