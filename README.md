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

# A. Published data.

## 0. Download public data with SRAtools

          # Download raw data
          for i in $(echo SRX18026584 SRX18026583 SRX18026582 SRX18026581 SRX18026580 SRX18026579 SRX18026578 SRX18026577 SRX18026576 SRX18026575 SRX18026574 SRX18026573 SRX18026572 SRX18026571 SRX18026570 SRX18026569 SRX18026568 SRX18026567 SRX18026566 SRX18026565 SRX18026564 SRX18026563 SRX18026562 SRX18026561 SRX18026560 SRX18026559 SRX18026558 SRX18026557 SRX18026556 SRX18026555); do echo $i; sbatch --partition=pshort_el8 --job-name=SRAtools1--time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools1.out --error=SRAtools1.error --mail-type=END,FAIL --wrap "module load SRA-Toolkit;cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; prefetch $i"; done

          # mv to parent folder
          mv SRR*/*.sra . 

          # split en paired fastq files
          for i in $(ls SRR*.sra); do echo $i;sbatch --partition=pshort_el8 --job-name=SRAtools2--time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools2.out --error=SRAtools2.error --mail-type=END,FAIL --wrap "module load SRA-Toolkit;cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; fastq-dump --split-files $(echo $i | cut -d'.' -f1)"; done

          # gzip fastq files
          for i in $(ls SRR*.fastq); do echo $i; sbatch --partition=pshort_el8 --job-name=gzip --time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools1.out --error=SRAtools1.error --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; gzip $i"; done
          
## 1. Trim data with fastp 

Clean reads with fastp

                        cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads;
                        for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-04:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads ; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz"; sleep  1; done

## 2. align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html

           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H2Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap2_HiSat2index.log --error=Hap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/40_S_spinosum_FinalFinal/03_BRAKER/Ref_RnaSeq/ref_Genomes; module load HISAT2; hisat2-build S_spinosum_hap2.fa.masked hap2_hisat_index -p 12"

   2. Map reads to genome, sort and compress



              # run from /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData
            sbatch --partition=pibu_el8 --job-name=H1Hisatmap3 --time=3-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/40_S_spinosum_FinalFinal/03_BRAKER/Ref_RnaSeq/ref_Genomes/hap2_hisat_index -1 SRR22045430_1_trimmed.fastq.gz,SRR22045431_1_trimmed.fastq.gz,SRR22045432_1_trimmed.fastq.gz,SRR22045433_1_trimmed.fastq.gz,SRR22045434_1_trimmed.fastq.gz,SRR22045435_1_trimmed.fastq.gz,SRR22045436_1_trimmed.fastq.gz,SRR22045437_1_trimmed.fastq.gz,SRR22045438_1_trimmed.fastq.gz,SRR22045439_1_trimmed.fastq.gz,SRR22045440_1_trimmed.fastq.gz,SRR22045441_1_trimmed.fastq.gz,SRR22045442_1_trimmed.fastq.gz,SRR22045443_1_trimmed.fastq.gz,SRR22045444_1_trimmed.fastq.gz,SRR22045445_1_trimmed.fastq.gz,SRR22045446_1_trimmed.fastq.gz,SRR22045447_1_trimmed.fastq.gz,SRR22045448_1_trimmed.fastq.gz,SRR22045449_1_trimmed.fastq.gz,SRR22045450_1_trimmed.fastq.gz,SRR22045451_1_trimmed.fastq.gz,SRR22045452_1_trimmed.fastq.gz,SRR22045453_1_trimmed.fastq.gz,SRR22045454_1_trimmed.fastq.gz,SRR22045455_1_trimmed.fastq.gz,SRR22045456_1_trimmed.fastq.gz,SRR22045457_1_trimmed.fastq.gz,SRR22045458_1_trimmed.fastq.gz -2 SRR22045430_2_trimmed.fastq.gz,SRR22045431_2_trimmed.fastq.gz,SRR22045432_2_trimmed.fastq.gz,SRR22045433_2_trimmed.fastq.gz,SRR22045434_2_trimmed.fastq.gz,SRR22045435_2_trimmed.fastq.gz,SRR22045436_2_trimmed.fastq.gz,SRR22045437_2_trimmed.fastq.gz,SRR22045438_2_trimmed.fastq.gz,SRR22045439_2_trimmed.fastq.gz,SRR22045440_2_trimmed.fastq.gz,SRR22045441_2_trimmed.fastq.gz,SRR22045442_2_trimmed.fastq.gz,SRR22045443_2_trimmed.fastq.gz,SRR22045444_2_trimmed.fastq.gz,SRR22045445_2_trimmed.fastq.gz,SRR22045446_2_trimmed.fastq.gz,SRR22045447_2_trimmed.fastq.gz,SRR22045448_2_trimmed.fastq.gz,SRR22045449_2_trimmed.fastq.gz,SRR22045450_2_trimmed.fastq.gz,SRR22045451_2_trimmed.fastq.gz,SRR22045452_2_trimmed.fastq.gz,SRR22045453_2_trimmed.fastq.gz,SRR22045454_2_trimmed.fastq.gz,SRR22045455_2_trimmed.fastq.gz,SRR22045456_2_trimmed.fastq.gz,SRR22045457_2_trimmed.fastq.gz,SRR22045458_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/40_S_spinosum_FinalFinal/03_BRAKER/Ref_RnaSeq/02_PublishedData/PRJNA863910_30samples_Hap2.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/01_hap1/; module load SAMtools; samtools view --threads 48 -b -o PRJNA863910_30samples_Hap1.bam PRJNA863910_30samples_Hap1.sam; samtools sort -o PRJNA863910_30samples_Hap1_sorted.bam -T PRJNA863910_30samples_Hap1_temp --threads 48 PRJNA863910_30samples_Hap1.bam"






# B. Own data

## 1. Trim data with fastp 

Clean reads with fastp

                        for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-02:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/00_rawreads ; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz"; sleep  1; done



## 2. align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html

           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H1Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2-build hap1.fasta.masked hap1_hisat_index -p 12"

   2. Map reads to genome, sort and compress


                sbatch --partition=pibu_el8 --job-name=FrH1Hisatmap3 --time=3-21:00:00 --mem-per-cpu=16G --ntasks=48 --cpus-per-task=1 --output=FrHap1_HiSat2index.log --error=FrHap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/01_hap1/hap1_hisat_index -1 1A_1_trimmed.fastq.gz,2A_1_trimmed.fastq.gz,3A_1_trimmed.fastq.gz,4A_1_trimmed.fastq.gz,5A_1_trimmed.fastq.gz,6A_1_trimmed.fastq.gz,7A_1_trimmed.fastq.gz,8A_1_trimmed.fastq.gz,9A_1_trimmed.fastq.gz,10A_1_trimmed.fastq.gz,11A_1_trimmed.fastq.gz,12A_1_trimmed.fastq.gz -2 1A_2_trimmed.fastq.gz,2A_2_trimmed.fastq.gz,3A_2_trimmed.fastq.gz,4A_2_trimmed.fastq.gz,5A_2_trimmed.fastq.gz,6A_2_trimmed.fastq.gz,7A_2_trimmed.fastq.gz,8A_2_trimmed.fastq.gz,9A_2_trimmed.fastq.gz,10A_2_trimmed.fastq.gz,11A_2_trimmed.fastq.gz,12A_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/01_hap1/Ufribourg_12samples_Hap1.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=FrH2SAMTOOLS --time=0-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=FrHap2_SAMTOOLS.log --error=FrHap2_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/02_hap2/; module load SAMtools; samtools view --threads 48 -b -o Ufribourg_12samples_Hap2.bam Ufribourg_12samples_Hap2.sam; samtools sort -m 64G -o Ufribourg_12samples_Hap2_sorted.bam -T Ufribourg_12samples_Hap2_temp --threads 48 Ufribourg_12samples_Hap2.bam"


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
