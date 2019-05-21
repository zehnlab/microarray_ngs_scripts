library(pheatmap)
library(openxlsx)
library(ggplot2)

#import data
Filename <- c('Run29_D8','Run32_Day20','Run32_Day8_TIM3p','Run32_Day8_TIM3n')


Heatmap.Func <- function(Heatmap_Matrix, Annotation_Table, SelectGene, Filename, Gaps, Height, Width){
  
  # re-order Heatmap 
  Heatmap_df <- NULL
  Heatmap_df <- as.data.frame(lapply(levels(Annotation_Table$Group), function (x) {
    rbind(Heatmap_df, Heatmap_Matrix[rownames(Heatmap_Matrix) %in% SelectGene,
                                     rownames(subset(Annotation_Table, Group == x))])
  }))
  
  pdf(paste0(Filename,'.pdf'), onefile = F, height = Height, width = Width)
  pheatmap::pheatmap(Heatmap_df,
                     show_colnames = T, cluster_rows = T, cluster_cols = F,show_rownames = T, 
                     cellheight = 10, cutree_cols = 2, gaps_col = Gaps,fontsize = 6,
                     annotation = Annotation_Table['Group'])
  
  dev.off()
}


Data.Prepare <- function(Heatmap_Matrix_1, Annotation_Table_1, Heatmap_Matrix_2 = NULL, 
                         Annotation_Table_2 = NULL, Scale = FALSE,
                         SelectGene, Filename, Gaps, Height, Width){
  
  # New Annotation table
  Annotation_Table <- rbind(Annotation_Table_1['Group'], Annotation_Table_2['Group'])
  # Scale or Merge dataset
  if(Scale){Heatmap_Matrix <- t(scale(t(Heatmap_Matrix_1)))}
  Heatmap.Func(Heatmap_Matrix, Annotation_Table, SelectGene, Filename, Gaps, Height, Width)
}


for(file in Filename){
  
  RlogTable <- read.csv(list.files(path = 'output/DEG/',pattern = paste0('Rlog_',file), 
                                     recursive=TRUE, full.names = T), sep = '\t', row.names = 1, header = T)
  Annotation_Table <- read.csv(list.files(path = 'data/',pattern = paste0(file,'.*Design*.'), 
                                          recursive=TRUE, full.names = T), sep = '\t', row.names = 1, header = T)
  
  DEG_List <- read.xlsx(list.files(path = 'output/DEG/',pattern = paste0('Filter_DEG_',file), 
                                  recursive=TRUE, full.names = T))
  #We are not interested in Ig genes
  SelectGene <- DEG_List[!grepl(paste0(c("(Tr[a-b]+)","(Ig[h-k]+)"), collapse = '|'),DEG_List$geneName, perl=TRUE),'geneName']
    
  Data.Prepare(RlogTable, Annotation_Table, Heatmap_Matrix_2 = NULL, 
               Annotation_Table_2 = NULL, Scale = T,
               SelectGene, paste0('output/Heatmaps/Rlog_Scale_',file), Gaps = c(5), Height = ifelse(length(SelectGene) > 150, 55, 40), Width = 4)
}


