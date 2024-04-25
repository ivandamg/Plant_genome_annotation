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

## 1. Trim data with fastp 

Clean reads with fastp

                        for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-02:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/00_rawreads ; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz"; sleep  1; done



## 2. align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html

           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H1Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2-build hap1.fasta.masked hap1_hisat_index -p 12"

   2. Map reads to genome, sort and compress


                sbatch --partition=pibu_el8 --job-name=H1Hisatmap3 --time=3-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2 --phred33 --dta -x hap1_hisat_index -1 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_1_trimmed.fastq.gz -2 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_2_trimmed.fastq.gz -S threeControlSamples_Hap1.sam -p 12"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1__SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load SAMtools; samtools view --threads 12 -b -o threeControlSamples_Hap1.bam threeControlSamples_Hap1.sam; samtools sort -m 7G -o threeControlSamples_Hap1_sorted.bam -T threeControlSamples_Hap1_temp --threads 12 threeControlSamples_Hap1.bam"


## BRAKER3 annotation with GALAXY

Run BRAKER3pipeline in galaxy it needs as input a genome assembly, proteins from related species and mapped RNAseq reads to the assembly.

it gives as output a .gff3 file

convert gff to transcripts and proteins with gffread on local

        #convert to proteins
        /Users/mateusgo/ZZ_Software/gffread/gffread -y proteins.fa -g hap1.fasta.masked hap1_Braker.gff3

        # convert to transcripts
        /Users/mateusgo/ZZ_Software/gffread/gffread -w transcripts.fa -g hap1.fasta.masked hap1_Braker.gff3

Then busco the annotation on cluster

            sbatch --partition=pibu_el8 --job-name=LFh1protBusco --time=0-04:00:00 --mem-per-cpu=50G --ntasks=36 --cpus-per-task=1 --output=BuscoFunHap1.out --error=BuscoFunHap1.error --mail-type=END,FAIL --wrap "module load BUSCO; cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/07_Busco/v2/; busco -o Funannotate_hap1 -i Argania_spinosa_DR.proteins.fa -l embryophyta_odb10 --cpu 36 -m proteins -f"


sbatch --partition=pibu_el8 --job-name=BRAKER3h1protBusco --time=0-04:00:00 --mem-per-cpu=50G --ntasks=36 --cpus-per-task=1 --output=BuscoBRAKERHap1.out --error=BuscoBRAKERHap1.error --mail-type=END,FAIL --wrap "module load BUSCO; cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/07_Busco/v2/; busco -o BRAKER3_hap1 -i hap1_Braker_proteins.fa -l embryophyta_odb10 --cpu 36 -m proteins -f"










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
