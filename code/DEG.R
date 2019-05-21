library(DESeq2)
library(openxlsx)

Filename <- c('Run29_D8','Run32_Day20','Run32_Day8_TIM3p','Run32_Day8_TIM3n')

#make sure you have installed above packages.:D
Deseq.Func <- function(ReadCountTable, DesignTable, Filename, Treat, Ctr, 
                       BaseMean= 50, logFC = 1, Padj = 0.05){
  
  #gene pre-filtering
  ReadCountTable <- ReadCountTable[rowSums(ReadCountTable) > 10,]
  
  #perform DEGA
  dds <- DESeqDataSetFromMatrix(countData = ReadCountTable,
                                colData = DesignTable,
                                design = ~ Condition)
  dds <- DESeq(dds)
  
  #get raw DEG list
  res <- as.data.frame(results(dds, contrast = c("Condition", Treat, Ctr)))
  res$geneName <- rownames(res)
  
  #export norm/rlog, rawDEG data
  write.table(x = rownames_to_column(as.data.frame(counts(dds, normalized = T)),'Gene'), 
            file = paste0('output/DEG/','Norm_', Filename, '.csv'), sep = '\t', row.names = F, quote = F)

  write.table(x = rownames_to_column(as.data.frame(assay(rlog(dds, blind = F))),'Gene'),
            file = paste0('output/DEG/','Rlog_', Filename,'.csv'), sep = '\t', row.names = F, quote = F)
  
  write.xlsx(x = na.omit(res), file = paste0('output/DEG/','Raw_DEG_', Filename,'_', Treat, '_', Ctr, '.xlsx'))
  write.xlsx(x = na.omit(res[which(res$baseMean >= BaseMean & abs(res$log2FoldChange) >= logFC &
                                     res$padj < Padj),]), file = paste0('output/DEG/','Filter_DEG_', 
                                                                        Filename, '_', Treat, '_', Ctr,'.xlsx'))
  
}


for(file in Filename){
  ReadCountTable <- read.csv(list.files(path = 'data/',pattern = paste0(file,'.*Count*.'), 
                               recursive=TRUE, full.names = T), sep = '\t', row.names = 1, header = T)
  DesignTable <- read.csv(list.files(path = 'data/',pattern = paste0(file,'.*Design*.'), 
                                        recursive=TRUE, full.names = T), sep = '\t', row.names = 1, header = T)
  
  Deseq.Func(ReadCountTable, DesignTable, file, 'WT', 'KO', 
             BaseMean= 50, logFC = 1, Padj = 0.05)
}

