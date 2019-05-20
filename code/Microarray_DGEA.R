library(limma)


#read files
inputFiles_array1 <- list.files(path = '/home/angela/small_ex/RS-313_Array1', 
                                pattern='RS-313_*', recursive=TRUE, full.names = T)
inputFiles_array2 <- list.files(path = '/home/angela/small_ex/RS-350_Array2', 
                                pattern='RS-350_*', recursive=TRUE, full.names = T)

target_file_array1 <- read.table('~/small_ex/RS-313_Array1/Array1_Design.txt', header = T, sep = '\t')
target_file_array2 <- read.table('~/small_ex/RS-350_Array2/Array2_Design.txt', header = T, sep = '\t')


MicroExp_Func <- function(target_file, filepath, source, correct_method, normalize_method, offset){
  
  inputData <- read.maimages(sapply(target_file[1],function(x) paste0(x, '.txt')), path=filepath, source=source, green.only = T,
                             annotation=c("Row", "Col", "ProbeName", "SystematicName","GeneName"))
  tmp <- inputData$genes$GeneName %in% 
    inputData$genes[which(inputData$genes$GeneName == inputData$genes$ProbeName | inputData$genes$GeneName == inputData$genes$SystematicName),'GeneName']
  
  Filter_input <- new('EListRaw', list(targets = inputData$targets, genes = subset(inputData$genes, !tmp),
                                     source = inputData$source, E = subset(inputData$E, !tmp),
                                     Eb = subset(inputData$Eb, !tmp)))
  
  corretedData <- limma::backgroundCorrect(Filter_input, method=correct_method, offset = offset)
  normData <- normalizeBetweenArrays(corretedData, method=normalize_method)
  normData.ave <- avereps(normData, ID=normData$genes$ProbeName)
  
  # design <- model.matrix(~ 0 + Infection + GeneArraySet, data=target_file, 
  #                        contrasts.arg=list(Infection=diag(nlevels(target_file$Infection)), 
  #                                           GeneArraySet=diag(nlevels(target_file$GeneArraySet))))
  design <- model.matrix( ~ 0 + target_file$Infection)
  colnames(design) <- levels(target_file$Infection)
  fit <- lmFit(normData.ave, design)
  
  contrast.matrix <- makeContrasts('Chronic-Acute',levels=c('Acute','Chronic'))
  fit.eBayes <- eBayes(contrasts.fit(fit, contrasts = contrast.matrix))
  plotSA(fit.eBayes)
  print(fit.eBayes$df.prior)
  output <- topTable(fit.eBayes, adjust="BH", coef="Chronic-Acute", genelist=normData.ave$genes, number=Inf,
                     sort.by = 'P')
  
  output <- aggregate(. ~ GeneName, mean, data = output[,c('GeneName','logFC','AveExpr','t',
                                           'P.Value','adj.P.Val','B')])
  
  return(output)
  
}

Array2_Chronic_Acute <- MicroExp_Func(target_file_array2[target_file_array2$Day == 'd28',], "/home/angela/small_ex/RS-350_Array2", 
                                      "agilent.median","normexp","quantile",16)
write.xlsx(Array2_Chronic_Acute, '~/MING_V9T/Microarray/D28_Array2_DEGList_Chronic_Acute_raw.xlsx')

Array2_DEG <- as.data.frame(Array2_Chronic_Acute[which(abs(Array2_Chronic_Acute$logFC) > 1 & Array2_Chronic_Acute$adj.P.Val < 0.05),])

write.table(Array2_DEG[order(Array2_DEG$logFC, decreasing = T),], "RS-350_Array2/D28_Array2_Chronic_Acute_DEGList.csv",
            sep = '\t', row.names=F, quote = F)


Array1_Chronic_Acute <- MicroExp_Func(target_file_array1, "/home/angela/small_ex/RS-313_Array1",
                          "agilent.median","normexp","quantile",16)

write.xlsx(Array1_Chronic_Acute, '~/MING_V9T/Microarray/D28+D8_Array1_DEGList_Chronic_Acute_raw.xlsx')
Array1_DEG <- as.data.frame(Array1_Chronic_Acute[which(abs(Array1_Chronic_Acute$logFC) > 1 & Array1_Chronic_Acute$adj.P.Val < 0.05),])

write.table(Array1_DEG[order(Array1_DEG$logFC, decreasing = T),], "RS-313_Array1/D28+D8_Array1_Chronic_Acute_DEGList.csv",
            sep = '\t', row.names=F, quote = F)

