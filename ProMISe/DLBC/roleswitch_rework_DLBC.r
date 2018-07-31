#roleswitch rework

setwd('/home/node05/바탕화면/miRNA target performance/ProMISe/DLBC')

library(PBE)
library(ROCR)
library(preprocessCore)
library(Roleswitch)


devide_val = function(table,val_list)
{
  val_collection = list()
  for(i in 1:11)
  {
    val = val_list[[i]]
    #ROC and AUC calculated
    gpval = val[intersect(getGoldPositive(novel_miRs[i]),rownames(table))]
    gnval = val[intersect(getGoldNegative(novel_miRs[i]),rownames(table))]
    val = c(gpval, gnval)
    
    val_collection[[novel_miRs[i]]] = val
  }
  return(val_collection)
}

devide_label = function(table,val_list)
{
  label_collection = list()
  for(i in 1:11)
  {
    val = val_list[[i]]
    #ROC and AUC calculated
    gpval = val[intersect(getGoldPositive(novel_miRs[i]),rownames(table))]
    gnval = val[intersect(getGoldNegative(novel_miRs[i]),rownames(table))]
    val = 1-c(gpval, gnval)
    label = array(0, dim = length(val))
    label[1:length(gpval)]=1
    label = factor(label)
    
    label_collection[[novel_miRs[i]]] = label
  }
  return(label_collection)
}

#function for calculate raw values of ROC
draw_AUC = function(table,val,method_type,miR)
{
  #ROC and AUC calculated
  gpval = val[intersect(getGoldPositive(miR),rownames(sq_mean))]
  gnval = val[intersect(getGoldNegative(miR),rownames(sq_mean))]
  val = c(gpval, gnval)
  label = array(0, dim = length(val))
  label[1:length(gpval)]=1
  label = factor(label)
  if(method_type=='negative')
  {
    pred = prediction(predictions = 1-val, labels = label)
  }
  else if(method_type=='positive')
  {
    pred = prediction(predictions = val, labels = label)
  }
  perf1 = performance(pred, 'tpr', 'fpr')
  X_FPR = as.numeric(perf1@x.values[[1]])
  Y_TPR = as.numeric(perf1@y.values[[1]])
  perf2 = performance(pred, 'auc')
  AUC = perf2@y.values[[1]]
  roleswitch_test = list(perf1=perf1,X_FPR = X_FPR, Y_TPR = Y_TPR,perf2 = perf2,AUC=AUC)
  
}

getGoldPositive = function (miRNA) 
{
  data(mtb)
  seqtar = getBindingTarget(miRNA)
  GP1 = mtb$Target.Gene[which(mtb$miRNA == miRNA & mtb$Support.Type == "Functional MTI")]
  GP1 = intersect(seqtar, GP1)
  GP2 =mtb$Target.Gene[which(mtb$miRNA == miRNA & mtb$Support.Type == "Non-Functional MTI")]
  GP2 = intersect(seqtar, GP2)
  GP = unique(c(GP1, GP2))
  return(GP)
}

getNewGoldPositive = function(miRNA)
{
  GP = as.character(read.table(paste('/home/node05/바탕화면/miRNA target performance/ProMISe/new_Golds/',miRNA,'_GoldPositives.txt',sep = ''))[,1])
  return(GP)
}
getNewGoldNegative = function(miRNA)
{
  GN = as.character(read.table(paste('/home/node05/바탕화면/miRNA target performance/ProMISe/new_Golds/',miRNA,'_GoldNegatives.txt',sep = ''))[,1])
  return(GN)
}

mRNA_data = read.delim(file = '/home/node05/바탕화면/miRNA target performance/ProMISe/DLBC/RSEM_genes_normalized__data/RSEM_genes_normalized_data.txt',sep = '\t')
edit_table = mRNA_data[2:length(mRNA_data[,1]),]
Gene.Symbol = as.character(edit_table$Hybridization.REF)

A=strsplit(Gene.Symbol[1], split = "|", fixed = T)
A[[1]][1]
GeneSymbol=sapply(Gene.Symbol, FUN = function(x){A=strsplit(x, split = "|", fixed = T); A[[1]][1]}, simplify = T, USE.NAMES = F)
edit_table=cbind(GeneSymbol,edit_table[,2:length(edit_table[1,])])
edit_table = edit_table[!(edit_table$GeneSymbol=="?"),]
GeneSymbol = as.character(edit_table$GeneSymbol)
edit_table_1 = subset(x = edit_table, select = -1)
edit_table_2 = apply(edit_table_1, 2, as.numeric)
rownames(edit_table_2) = GeneSymbol
colnames(edit_table_2) = gsub(".","-",colnames(edit_table_2),fixed = T) #find and change
edit_table_2 = edit_table_2[, order(colnames(edit_table_2))]

miRNA_body = read.delim(file = "mature_miRNA_exp.txt",header = F)
miRNA_col = readLines("sampleID.txt")
miRNA_row = readLines("miRNAID.txt")
miRNA_body = data.matrix(miRNA_body)

rownames(miRNA_body) = miRNA_row
colnames(miRNA_body) = miRNA_col
miRNA_body = miRNA_body[, order(colnames(miRNA_body))]



a = unlist(sapply(colnames(edit_table_2), function(x) unlist(strsplit(x,split="-"))[3],simplify = T,USE.NAMES = F))
b = unlist(sapply(colnames(miRNA_body), function(x) unlist(strsplit(x,split="-"))[3],simplify = T,USE.NAMES = F))
edit_table_2 = edit_table_2[,-29]

id_convert = read.delim("miRNAid_convert.txt",header = F)
rownames(miRNA_body)

MIMAT = as.character(id_convert$V1)
miR = as.character(id_convert$V2)

MIMAT2MIR = function(x)
{
  idx = which(MIMAT == x)
  if(length(idx)==0){return(x)}
  return(miR[idx])
}

newRow = unlist(sapply(rownames(miRNA_body), function(x){MIMAT2MIR(x)}, simplify = T, USE.NAMES = F))
rownames(miRNA_body) = newRow

novel_miRs = c("155-5p", "29a-3p", "34a-5p", "125a-5p","145-5p","29b-3p","1-3p", "21-5p","29c-3p","221-3p","204-5p")
novel_miRs = paste("hsa-miR-",novel_miRs,sep = "")
novel_miR_table = miRNA_body[novel_miRs,]

Gold_Positives = c()
Gold_Negatives = c()
NG_Positives = c()
NG_Negatives = c()
for(i in novel_miRs)
{
  Gold_Positives = c(Gold_Positives,getGoldPositive(i))
  Gold_Negatives = c(Gold_Negatives,getGoldNegative(i))
  NG_Positives = c(NG_Positives,getNewGoldPositive(i))
  NG_Negatives = c(NG_Negatives,getNewGoldNegative(i))
}
Gold_Positives = unique(Gold_Positives)
Gold_Negatives = unique(Gold_Negatives)
NG_Positives = unique(NG_Positives)
NG_Negatives = unique(NG_Negatives)

PARTmiR = t(novel_miR_table)
PARTGP = t(edit_table_2[intersect(Gold_Positives,rownames(edit_table_2)),])
PARTGN = t(edit_table_2[intersect(Gold_Negatives,rownames(edit_table_2)),])
PARTNP = t(edit_table_2[intersect(NG_Positives,rownames(edit_table_2)),])
PARTNN = t(edit_table_2[intersect(NG_Negatives,rownames(edit_table_2)),])
PARTmRNA = cbind(PARTGP,PARTGN)
PARTmRNA = PARTmRNA[,unique(colnames(PARTmRNA))]
PARTmRNA_News = cbind(PARTNP,PARTNN)
PARTmRNA_News = PARTmRNA_News[,unique(colnames(PARTmRNA_News))]
                  
                  
#quantile normalize data
QN_miR = normalize.quantiles(t(PARTmiR)) # Sora: Must be transposed
QN_miR = t(QN_miR) # Sora: Added by Sora
colnames(QN_miR) = colnames(PARTmiR)
rownames(QN_miR) = rownames(PARTmiR)
QN_mRNA = normalize.quantiles(t(PARTmRNA)) # Sora: Must be transposed
QN_mRNA = t(QN_mRNA) # Sora: Added by Sora
colnames(QN_mRNA) = colnames(PARTmRNA)
rownames(QN_mRNA) = rownames(PARTmRNA)
QN_mRNA_News = normalize.quantiles(t(PARTmRNA_News)) # Sora: Must be transposed
QN_mRNA_News = t(QN_mRNA_News) # Sora: Added by Sora
colnames(QN_mRNA_News) = colnames(PARTmRNA_News)
rownames(QN_mRNA_News) = rownames(PARTmRNA_News)

c = matrix(0,nrow = ncol(QN_mRNA),ncol = ncol(QN_miR))
colnames(c) = colnames(QN_miR)
rownames(c) = colnames(QN_mRNA)

c_New = matrix(0,nrow = ncol(QN_mRNA_News),ncol = ncol(QN_miR))
colnames(c_New) = colnames(QN_miR)
rownames(c_New) = colnames(QN_mRNA_News)

for(mir in colnames(c))
{
  c[intersect(c(getGoldPositive(mir),getGoldNegative(mir)),rownames(c)),mir] = 1
}

for(mir in colnames(c_New))
{
  c_New[intersect(c(getNewGoldPositive(mir),getNewGoldNegative(mir)),rownames(c_New)),mir] = 1
}


rownames(c) = c(1:length(colnames(QN_mRNA)))
colnames(c) = c(1:length(colnames(QN_miR)))
rownames(c_New) = c(1:length(colnames(QN_mRNA_News)))
colnames(c_New) = c(1:length(colnames(QN_miR)))

mother_list = list()
log_mother_list = list()
new_mother_list = list()
new_log_mother_list = list()

for(i in 1:length(rownames(QN_miR)))
{
  x.o = matrix(QN_mRNA[i,],dimnames = list(c(1:length(colnames(QN_mRNA))),'mRNA'))
  z.o = matrix(QN_miR[i,],dimnames = list(c(1:length(colnames(QN_miR))),'miRNA'))
  mother_list[[substr(rownames(QN_miR)[i],1,25)]] = roleswitch(x.o,z.o,c)$p.xz
  x.o = log2(x.o+1) # Sora: Add pseudocount 1
  z.o = log2(z.o+1)
  log_mother_list[[substr(rownames(QN_miR)[i],1,25)]] = roleswitch(x.o,z.o,c)$p.xz
  
  nx.o = matrix(QN_mRNA_News[i,],dimnames = list(c(1:length(colnames(QN_mRNA_News))),'mRNA'))
  nz.o = matrix(QN_miR[i,],dimnames = list(c(1:length(colnames(QN_miR))),'miRNA'))
  new_mother_list[[substr(rownames(QN_miR)[i],1,25)]] = roleswitch(nx.o,nz.o,c_New)$p.xz
  nx.o = log2(nx.o+1) # Sora: Add pseudocount 1
  nz.o = log2(nz.o+1)
  new_log_mother_list[[substr(rownames(QN_miR)[i],1,25)]] = roleswitch(nx.o,nz.o,c_New)$p.xz
}
                  
################
x = matrix(colMeans(QN_mRNA),dimnames = list(c(1:length(colnames(QN_mRNA))),'mRNA'))
z = matrix(colMeans(QN_miR),dimnames = list(c(1:length(colnames(QN_miR))),'miRNA'))
miRLAB_style = roleswitch(x,z,c)$p.xz
rownames(miRLAB_style) = colnames(QN_mRNA)
colnames(miRLAB_style) = colnames(QN_miR)
                  
mother_names = names(mother_list)
log_mother_names = names(log_mother_list)
new_mother_names = names(new_mother_list)
new_log_mother_names = names(new_log_mother_list)

sq_mean = matrix(0,length(rownames(c)),length(colnames(c)))
log_sq_mean = matrix(0,length(rownames(c)),length(colnames(c)))
new_sq_mean = matrix(0,length(rownames(c_New)),length(colnames(c_New)))
new_log_sq_mean = matrix(0,length(rownames(c_New)),length(colnames(c_New)))

for(i in mother_names)
{
  sq_mean = sq_mean + log2(mother_list[[i]]+1.0e-08)
  log_sq_mean = log_sq_mean + log2(log_mother_list[[i]]+1.0e-08)
  new_sq_mean = new_sq_mean + log2(new_mother_list[[i]]+1.0e-08)
  new_log_sq_mean = new_log_sq_mean + log2(new_log_mother_list[[i]]+1.0e-08)
}
sq_mean = 2^(sq_mean/length(mother_names))
rownames(sq_mean) = colnames(QN_mRNA)
colnames(sq_mean) = colnames(QN_miR)

log_sq_mean = 2^(log_sq_mean/length(log_mother_names))
rownames(log_sq_mean) = colnames(QN_mRNA)
colnames(log_sq_mean) = colnames(QN_miR)

new_sq_mean = 2^(new_sq_mean/length(new_mother_names))
rownames(new_sq_mean) = colnames(QN_mRNA_News)
colnames(new_sq_mean) = colnames(QN_miR)

new_log_sq_mean = 2^(new_log_sq_mean/length(new_log_mother_names))
rownames(new_log_sq_mean) = colnames(QN_mRNA_News)
colnames(new_log_sq_mean) = colnames(QN_miR)


val_list = list()
log_val_list = list()
new_val_list = list()
new_log_val_list = list()

for(i in 1:11)
{
  val_miRLAB[[i]] = miRLAB_style[intersect(c(getGoldPositive(novel_miRs[i]),getGoldNegative(novel_miRs[i])),rownames(miRLAB_style)),i]
  val_list[[i]] = sq_mean[intersect(c(getGoldPositive(novel_miRs[i]),getGoldNegative(novel_miRs[i])),rownames(sq_mean)),i]
  log_val_list[[i]] = log_sq_mean[intersect(c(getGoldPositive(novel_miRs[i]),getGoldNegative(novel_miRs[i])),rownames(sq_mean)),i]
  new_val_list[[i]] = new_sq_mean[intersect(c(getNewGoldPositive(novel_miRs[i]),getNewGoldNegative(novel_miRs[i])),rownames(new_sq_mean)),i]
  new_log_val_list[[i]] = new_log_sq_mean[intersect(c(getNewGoldPositive(novel_miRs[i]),getNewGoldNegative(novel_miRs[i])),rownames(new_sq_mean)),i]
}


AUCs = c()
log_AUCs = c()
new_AUCs = c()
new_log_AUCs = c()

for(i in 1:11)
{
  A = draw_AUC(sq_mean,val_list[[i]],'positive',novel_miRs[i])
  B = draw_AUC(log_sq_mean,log_val_list[[i]],'positive',novel_miRs[i])
  C = draw_AUC(new_sq_mean,new_val_list[[i]],'positive',novel_miRs[i])
  D = draw_AUC(new_log_sq_mean,new_log_val_list[[i]],'positive',novel_miRs[i])
  
  png(paste(novel_miRs[i],'_ROC.png',sep = ''))
  plot(A$perf1)
  abline(0,1)
  AUCs = c(AUCs,A$AUC)
  title(main = paste(novel_miRs[i],' ROC',sep = ''))
  text(0.2,1,paste('AUC: ',round(A$AUC,digits = 6),sep = ''))
  dev.off()
  
  png(paste(novel_miRs[i],'_ROC_log2.png',sep = ''))
  plot(B$perf1)
  abline(0,1)
  log_AUCs = c(log_AUCs,B$AUC)
  title(main = paste('log2 ',novel_miRs[i],' ROC',sep = ''))
  text(0.2,1,paste('AUC: ',round(B$AUC,digits = 6),sep = ''))
  dev.off()
  
  png(paste(novel_miRs[i],'_ROC_New_Golds.png',sep = ''))
  plot(C$perf1)
  abline(0,1)
  new_AUCs = c(new_AUCs,C$AUC)
  title(main = paste(novel_miRs[i],' ROC (New Golds)',sep = ''))
  text(0.2,1,paste('AUC: ',round(C$AUC,digits = 6),sep = ''))
  dev.off()
  
  png(paste(novel_miRs[i],'_ROC_New_Golds_log2.png',sep = ''))
  plot(D$perf1)
  abline(0,1)
  new_log_AUCs = c(new_log_AUCs,D$AUC)
  title(main = paste('log2 ',novel_miRs[i],' ROC (New Golds)',sep = ''))
  text(0.2,1,paste('AUC: ',round(D$AUC,digits = 6),sep = ''))
  dev.off()
  
}
max(AUCs)
min(AUCs)

max(log_AUCs)
min(log_AUCs)

max(new_AUCs)
min(new_AUCs)

max(new_log_AUCs)
min(new_log_AUCs)

AUC_Grade = function(X)
{
  Y = c()
  for(i in 1:length(X))
  {
    if(X[i]>0.9){Y[i] = 'A'}
    else if(X[i]>0.8&&X[i]<=0.9){Y[i] = 'B'}
    else if(X[i]>0.7&&X[i]<=0.8){Y[i] = 'C'}
    else if(X[i]>0.6&&X[i]<=0.7){Y[i] = 'D'}
    else{Y[i] = 'F'}
  }
  return(Y)
}
AUC_Grade(AUCs)
AUC_Grade(log_AUCs)

#drawing sum-up ROCs
mir_val = devide_val(miRLAB_style,val_miRLAB)
old_raw_val = devide_val(sq_mean,val_list)
old_log_val = devide_val(log_sq_mean,log_val_list)
new_raw_val = devide_val(new_sq_mean,new_val_list)
new_log_val = devide_val(new_log_sq_mean,new_log_val_list)

mir_label = devide_label(miRLAB_style,val_miRLAB)
old_raw_label = devide_label(sq_mean,val_list)
old_log_label = devide_label(log_sq_mean,log_val_list)
new_raw_label = devide_label(new_sq_mean,new_val_list)
new_log_label = devide_label(new_log_sq_mean,new_log_val_list)

mir_pred = prediction(mir_val, mir_label)
old_raw = prediction(old_raw_val,old_raw_label)
old_log = prediction(old_log_val,old_log_label)
new_raw = prediction(new_raw_val,new_raw_label)
new_log = prediction(new_log_val,new_log_label)

tma = 0
tor = 0
tol = 0
tnr = 0
tnl = 0
for(i in 1:11)
{
  tma = tma + performance(mir_pred,'auc')@y.values[[i]][1]
  tor = tor + performance(old_raw,'auc')@y.values[[i]][1]
  tol = tol + performance(old_log,'auc')@y.values[[i]][1]
  tnr = tnr + performance(new_raw,'auc')@y.values[[i]][1]
  tnl = tnl + performance(new_log,'auc')@y.values[[i]][1]
}
mir_AUC = tma/11
old_raw_AUC = tor/11
old_log_AUC = tol/11
new_raw_AUC = tnr/11
new_log_AUC = tnl/11

png('ROC_sum_up_miRLAB_style.png')
plot(performance(mir_pred,'tpr','fpr'),col='gray')
par(new=T)
plot(performance(mir_pred,'tpr','fpr'),avg='horizontal',xlab='',ylab='')
abline(0,1,lty = "dotted")
text(0.2,1,paste('AUC: ',round(mir_AUC,digits = 6),sep = ''))
title(main = 'ROC sum up 3DB-Golds raw miRLAB style (11 miRs)')
dev.off()

png('ROC_sum_up_3DB_raw.png')
plot(performance(old_raw,'tpr','fpr'),col='gray')
par(new=T)
plot(performance(old_raw,'tpr','fpr'),avg='horizontal',xlab='',ylab='')
abline(0,1,lty = "dotted")
text(0.2,1,paste('AUC: ',round(old_raw_AUC,digits = 6),sep = ''))
title(main = 'ROC sum up 3DB-Golds raw (11 miRs)')
dev.off()

png('ROC_sum_up_3DB_log.png')
plot(performance(old_log,'tpr','fpr'),col='gray')
par(new=T)
plot(performance(old_log,'tpr','fpr'),avg='horizontal',xlab='',ylab='')
abline(0,1,lty = "dotted")
text(0.2,1,paste('AUC: ',round(old_log_AUC,digits = 6),sep = ''))
title(main = 'ROC sum up 3DB-Golds log (11 miRs)')
dev.off()

png('ROC_sum_up_2DB_raw.png')
plot(performance(new_raw,'tpr','fpr'),col='gray')
par(new=T)
plot(performance(new_raw,'tpr','fpr'),avg='horizontal',xlab='',ylab='')
abline(0,1,lty = "dotted")
text(0.2,1,paste('AUC: ',round(new_raw_AUC,digits = 6),sep = ''))
title(main = 'ROC sum up 2DB-Golds raw (11 miRs)')
dev.off()

png('ROC_sum_up_2DB_log.png')
plot(performance(new_log,'tpr','fpr'),col='gray')
par(new=T)
plot(performance(new_log,'tpr','fpr'),avg='horizontal',xlab='',ylab='')
abline(0,1,lty = "dotted")
text(0.2,1,paste('AUC: ',round(new_log_AUC,digits = 6),sep = ''))
title(main = 'ROC sum up 2DB-Golds log (11 miRs)')
dev.off() 
