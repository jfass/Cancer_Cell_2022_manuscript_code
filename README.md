# Cancer_Cell_2022_manuscript_code

Code associated with the Cancer Cell 2022 article "Transcriptomic profiles of neoantigen-reactive T cells in human gastrointestinal cancers" by Chunhong Zheng, et al.

Copyright Providence Health Services

This code is distributed under a dual license. For academic, non-commerical use of the sofware the GNU Lesser General Public License (LGPLv3) open source license may be used. For commercial use of the software, please contact Dr. Brady Bernard <brady.bernard@providence.org>.

Code examples have been edited for clarity and brevity.

# Somatic mutation calling

GATK best practices were followed for germline (blood) and tumor tissue Whole Exome Sequencing (WES) data (paired-end Illumina reads, 150 x 150bp), aligning with BWA MEM, followed by mutation calling using Mutect2, VarScan, Strelka, and SomaticSniper, and annotation with ANNOVAR and snpEff:

```bash
# BWA MEM (separately for NORMAL_ID and TUMOR_ID)
$(DOCKER_RUN) ppmp-bwa mem \
  -M \
  -R '@RG\tID:'$(SAMPLE)'\tSM:'$(SAMPLE)'\tPL:ILLUMINA\tLB:lib1\tPU:flowcell-barcode.lane' \
  $(REFSEQ) \
  > $(NORMAL_ID.unsorted.sam
# ... GATK Best Practices ... 
# Mutect2
$(DOCKER_RUN) ppmp-gatk4 gatk --java-options "-Xmx12g" Mutect2 \
  -R $(REFSEQ) \
  -I $(TUMOR_ID)_reads.bam \
  -I $(NORMAL_ID)_reads.bam \
  -tumor $(TUMOR_ID) -normal $(NORMAL_ID) \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --germline-resource $(REFDIR)/hg19/af-only-gnomad.hg19.sorted.vcf.gz \
  --panel-of-normals $(REFDIR)/hg19/mutect2_pon.vcf.gz \
  -O mutect_raw.vcf
# VarScan
$(DOCKER_RUN) ppmp-varscan \
  somatic \
  $(NORMAL_ID).pileup \
  $(TUMOR_ID).pileup \
  varscan/varscan \
  --tumor-purity .5 \
  --output-vcf 1 \
  --min-reads2 2 \
  --min-coverage 4 \
  --min-var-freq .05 \
  --strand-filter 0 \
# Strelka - commands from Docker entrypoint script
/opt/strelka/bin/configureStrelkaSomaticWorkflow.py \
  --tumorBam /tmp/tumor.bam --normalBam /tmp/normal.bam \
  --ref /tmp/reference/ref.fasta \
  --runDir /tmp/out \
  --exome
/tmp/out/runWorkflow.py -m local -j $(nproc)
# SomaticSniper
$(DOCKER_RUN) ppmp-somaticsniper \
  bam-somaticsniper \
  -L -G -F \
  vcf -f $(REFSEQ) \
  $(TUMOR_ID)_reads.bam $(NORMAL_ID)_reads.bam \
  somaticsniper_raw.vcf
# combine calls
$(DOCKER_RUN) ppmp-gatk \
  -T CombineVariants \
  -R $(REFSEQ) \
  -V:varscan varscan_raw.vcf \
  -V:mutect mutect_raw.vcf  \
  -V:strelka strelka_raw.vcf \
  -V:somaticsniper somaticsniper_raw.vcf \
  -o combined.vcf
# snpEff
$(DOCKER_RUN) ppmp-snpeff \
  java -Xmx4G -jar /opt/snpEff/snpEff.jar -v -t \
  -datadir $(SNPEFF_DATA) -nodownload \
  $(SNPEFF_REF)\
  combined.vcf \
  > combined.snpeff.vcf 
# ANNOVAR
$(DOCKER_RUN) ppmp-tnannovar convert2annovar.pl \
  -format vcf4old \
  combined.vcf \
  --withzyg \
  -outfile combined.av
ppmp-tnannovar table_annovar.pl \
  combined.av \
  $(ANNOVAR_DATA) -buildver hg19 -v \
  -out combined.sum \
  -remove -protocol refGene,knownGene,ensGene,snp138NonFlagged,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,cosmic70 -operation g,g,g,f,f,f,f,f,f,f -nastring NA
```

# HLA-calling from WES or RNA-Seq

```bash
# PHLAT - command from Docker entrypoint script
/usr/bin/python2.7 -O /opt/phlat-1.0/dist/PHLAT.py \
  -1 /tmp/sample_1.fastq.gz -2 /tmp/sample_2.fastq.gz \
  -index /tmp/reference/index4phlat \
  -b2url /opt/bowtie2-2.2.3/bowtie2 \
  -tag sample -p 16 -e /opt/phlat-1.0 -o /tmp/output
# xHLA - command from Docker entrypoint script
samtools view -h -F 3840 alignments.bam \
  | perl -ane '@B=split("",sprintf("%08b",$F[1]));
               print if $F[0]=~m/^@/ or $B[5] or $B[4] or ($F[2] eq "chr6" and $F[3]>28477796 and $F[3]<33448355) or $F[2]=~m/^chr6.*hap/;' \
  | samtools view -bS - \
  | samtools sort -n -@ 12 -m 16G - \
  | bedtools bamtofastq -i - -fq /tmp/output/R1.fq -fq2 /tmp/output/R2.fq 2> /tmp/output/b2fq.log
bwa mem -t 12 /tmp/ref/GRCh38.primary_assembly.genome.fa /tmp/output/R1.fq /tmp/output/R2.fq 2> /tmp/output/mem.log \
  | samtools view -b - \
  | samtools sort -@ 12 -m 16G -o /tmp/output/remap.bam -
samtools index /tmp/output/remap.bam
samtools view -b /tmp/output/remap.bam chr6:29844528-33100696 > /tmp/output/remap.hla.bam
samtools index /tmp/output/remap.hla.bam
python3 /opt/HLA/bin/run.py --sample_id sample --input_bam_path /tmp/output/remap.hla.bam --output_path /tmp/output
# HLA*LA
$(DOCKER_RUN) rna-hla-la \
  --BAM /tmp/alignments.bam \
  --graph ../graphs/PRG_MHC_GRCh38_withIMGT \
  --sampleID $(NEWSAMPLE) \
  --maxThreads 12 \
  --workingDir /tmp/output \
  --picard_sam2fastq_bin /opt/picard.jar 
```

# RNA-Seq

## Alignment

```bash
$(DOCKER_RUN) rna-star \
  --genomeDir $(REFDIR)/$(REF)/STAR_genome \
  --twopassMode Basic --readFilesIn \
  $(SAMPLE)_R1.fastq.gz $(SAMPLE)_R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix twopass/$(SAMPLE)_paired \
  --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --runThreadN 3 \
  --outSAMunmapped Within
```

## Transcript-level expression

```bash
# Cufflinks - command from Docker entrypoint script
/opt/cufflinks/cufflinks \
  -p 32 -G /Ref/ann.gtf --library-type fr-firststrand \
  --output-dir /tmp/ /work/$(SAMPLE).deduped.bam
```

# Single-cell (Smart-seq2)

## TraCeR

```bash
docker run --rm -v /home/ubuntu/:/scratch -w /scratch teichlab/tracer assemble \
    -p 64 -s Hsap ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    ${sample} tracer.results \
    > ${sample}.log
```

## MiXCR

```bash
mixcr analyze amplicon --species hsa --starting-material rna \
    --5-end v-primers \
    --3-end c-primers \
    --adapters adapters-present \
    --contig-assembly \
    --impute-germline-on-export \
    --report ${outdir}/${sample}.report \
    ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    ${outdir}/${sample} \
    > ${sample}.log
```

## Gene Expression (in Seurat)

```R
fn <- list.files( path="tracer.results", pattern="abundance.tsv$", recursive=T, full.names=T )
samples <- sub( "\\/.*$", "", sub( ".*?\\/", "", fn ) )  # derive cell names from folders created by TraCeR
Ns <- length( samples )
# get gene names, and exclude assembled TCR names, from first sample:
temp <- read.table( fn[1], header=T, sep="\t", colClasses=c("character",rep("numeric",4)), row.names=1 )
tcrs <- grep( "^TRACER", rownames( temp ) )
txids <- rownames( temp )[ -tcrs ]
Nt <- length( txids )
# pull in counts, omitting final rows (TCR assemblies that vary per sample), takes a few minutes
est_counts <- as.data.frame(
  do.call(
    cbind,
    lapply( fn,
            function(x) {
              read.table( x,
                          header=T,
                          sep="\t",
                          colClasses=c("character",rep("numeric",4))
                          )[1:Nt,4]
              }
            )
    )
  )
rownames( est_counts ) <- txids
colnames( est_counts ) <- samples
# pull in gene names
ENST_to_ENSG <- read.table( "ENST_to_ENSG.biomart.tsv", header=T, as.is=T, sep="\t" )
# next line takes a few minutes
counts.genes <- merge( x=est_counts, y=ENST_to_ENSG, by.x=0, by.y="Transcript.stable.ID", all.x=T, all.y=F )
colnames( counts.genes )[1] <- "Transcript.stable.ID"
# summarize to gene level; also takes a few minutes
counts <- do.call( rbind,
                   lapply( unique( counts.genes$Gene.stable.ID ),
                           function(x) {
                             colSums( counts.genes[ which(counts.genes$Gene.stable.ID == x), 2:(Ns+1) ] )
                             }
                        )
                 )
# add gene rownames
rownames( counts ) <- unique( counts.genes$Gene.stable.ID )
# find unidentified genes (no name, or NA):
toDel <- c( which( rownames( counts ) == '' ),
            which( is.na( rownames( counts ) ) )
            )
# ... and remove them:
counts <- counts[ -toDel, ]
# prepend Gene Symbols to Ensembl IDs where possible
ENSG_to_HGNC <- read.table( "ENSG_to_HGNC.biomart.tsv", header=T, as.is=T, sep="\t" )
temp <- ENSG_to_HGNC$HGNC.symbol[ match( rownames( counts ), ENSG_to_HGNC$Gene.stable.ID ) ]
rownames( counts ) <- paste( temp, rownames( counts ), sep=":" )
# finally, create Seurat object:
library( Seurat )
gex.raw <- CreateSeuratObject( counts )
```

# Single-cell (10X Genomics)

## RNA-Seq / CITE-Seq

```bash
./cellranger-3.1.0/cellranger count \
  --id=${sample}-count \
  --libraries=libraries.${sample}.csv \
  --feature-ref=feature-reference.csv \
  --transcriptome=refdata-cellranger-GRCh38-3.0.0 \
  --localcores=24 --localmem=128
```

```R
library( Seurat )
temp <- Read10X( "${sample}-count/outs/count/filtered_feature_bc_matrix" )[[1]]
gex <- CreateSeuratObject( temp )
```

## TCR-Seq

```bash
./cellranger-3.1.0/cellranger vdj \
  --id=${sample}-vdj \
  --fastqs=MKFQ/outs/fastq_path/H* \
  --reference=refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 \
  --sample=${sample} \
  --localcores=24 --localmem=128
```

