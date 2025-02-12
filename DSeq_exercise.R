if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("DESeq2")
install.packages("tidyverse")

library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1
# read in cts data
counts_data <- read.csv("counts_data.csv")
head(counts_data)

# read in sample data
col_Data <- read.csv("sample_info.csv")

# row names in colData match col names in counts?
all(colnames(counts_data) %in% rownames(col_Data))

# also in same order?
all(colnames(counts_data) == rownames(col_Data))

# Step 2: construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_Data,
                              design = ~ dexamethasone)

dds
# pre-filter, keep rows that hav at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
dds$dexamethasone

# run DESeq and save results
dds <- DESeq(dds)
res <- results(dds)

# explore results
res
summary(res)

# can change p-val used, say, to 0.01
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)
