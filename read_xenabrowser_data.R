# read TCGA-TARGET-GTEx data----

# read metadata
metadata <- read.csv("data/xenabrowser/TcgaTargetGTEX_phenotype.txt", sep = "\t")
metadata$sample <- gsub("-", ".", metadata$sample)
rownames(metadata) <- metadata$sample

dim(metadata)

# read countdata
counts <- read.csv("data/xenabrowser/TcgaTargetGtex_gene_expected_count.txt", sep = "\t")
dim(counts)
counts[1:10,1:10]
row.names(counts) <- counts$sample
counts <- counts[,-1]

samples <- intersect(colnames(counts), metadata$sample)
counts <- counts[,samples]
dim(counts)
metadata <- metadata[samples,]

all(rownames(metadata) %in% colnames(counts))
all(rownames(metadata) == colnames(counts))

# create Summarized Experiment object----
SE_xenabrowser <- SummarizedExperiment(assays = list(counts=counts), colData = metadata)
saveRDS(SE_xenabrowser, "data/xenabrowser/TcgaTargetGTEX_gene_expected_count.Rds")

# read SE object data----
data_file <- "data/xenabrowser/TcgaTargetGTEX_gene_expected_count.Rds"
SE_xenabrowser <- readRDS(data_file)

# subsetting data
SE_NB <- SE_xenabrowser[,colData(SE_xenabrowser)$detailed_category=="Neuroblastoma"]
SE_NB
saveRDS(SE_NB, "data/xenabrowser/NB_gene_expected_count.Rds")

SE_GTEx <- SE_xenabrowser[,colData(SE_xenabrowser)$X_study=="GTEX"]
SE_GTEx
saveRDS(SE_GTEx, "data/xenabrowser/GTEx_gene_expected_count.Rds")

SE_TARGET <- SE_xenabrowser[,colData(SE_xenabrowser)$X_study=="TARGET"]
SE_TARGET
saveRDS(SE_TARGET, "data/xenabrowser/TARGET_gene_expected_count.Rds")

SE_TCGA <- SE_xenabrowser[,colData(SE_xenabrowser)$X_study=="TCGA"]
SE_TCGA
saveRDS(SE_TCGA, "data/xenabrowser/TCGA_gene_expected_count.Rds")

