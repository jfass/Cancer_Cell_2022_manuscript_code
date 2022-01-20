# Cancer_Cell_2022_manuscript_code

Code associated with the Cancer Cell 2022 article "Transcriptomic profiles of neoantigen-reactive T cells in human gastrointestinal cancers" by Chunhong Zheng, et al.

Copyright Providence Health Services

This code is distributed under a dual license. For academic, non-commerical use of the sofware the GNU Lesser General Public License (LGPLv3) open source license may be used. For commercial use of the software, please contact Dr. Brady Bernard <brady.bernard@providence.org>.

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

