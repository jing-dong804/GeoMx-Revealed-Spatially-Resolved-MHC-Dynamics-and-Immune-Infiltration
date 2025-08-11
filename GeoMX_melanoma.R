
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[1] test the immune infiltration
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.1] BASE inference
## please note the renormalization by using DeduplicatedReads

rm(list=ls())
myDir1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/"
myoutf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

myinf = dir(myDir1)
myinf = myinf[grep("xlsx", myinf)]
myinf1 = myinf[grep("_ITX1_", myinf)]
myinf1 = paste(myDir1, myinf1, sep="")
myinf2 = myinf[grep("_ITX2_", myinf)]
myinf2 = paste(myDir1, myinf2, sep="")
library(gdata)
library(readxl)


subset1 = c("ROI_ID_2","ROICoordinateX", "ROICoordinateY",  "RawReads", "AlignedReads", "DeduplicatedReads", "TrimmedReads", 
"StitchedReads", "disease_control_", "TREATMENT", "LTB24", "SPEC_CATEGORY", 
"DCB6", "LTB12", "PROG_STATUS_BY_SCAN", "VITAL", "Core_Type", 
"response", "CPID", "RACE", "MUTATION", "AGE_AT_DX",
"OS_FROM_START_OF_ITX", "Prior_Checkpoint_Blockade", "SEX")

stage.col1 = c("IIA", "IIIC", "IVB (H&N)", "IIIB", "IV", "IB", "IIIA", "UNKNOWN", "III", "IIB", 
"I", "IIC", "IVA")
       
data.List = list(NULL)
clin.List = list(NULL)
for(k in 1:length(myinf1))
{
	info = read_excel(myinf1[k], sheet=1)
	info = as.data.frame(info)
	row.names(info) = gsub(" ", "", info$SegmentDisplayName)
	se = which(colnames(info)%in%stage.col1)
	xx = info[, se]
	for(i in 1:ncol(xx))
	{
		xx[,i] = ifelse(xx[,i]=="True", 1, 0)
	}
	tmp = apply(xx, 1, which.max)
	tmp = colnames(xx)[tmp]
	stage= tmp
	info = info[, subset1]
	info$stage = stage
	clin.List[[k]] = info 
	
	data = read_excel(myinf1[k], sheet=4)
	data = as.data.frame(data)
	tmp = data[,1]
	data = data[, -1]
	row.names(data) = tmp
	colnames(data) = gsub(" ", "", colnames(data) )
	data.List[[k]] = data
}

for(k in 1:length(data.List))
{
	data = data.List[[k]]
	nn = clin.List[[k]]$DeduplicatedReads
	for(i in 1:ncol(data))
	{
		data[,i] = data[,i]*1e6/nn[i]
	}
	data.List[[k]] = data
}
clin.List1 = clin.List
data.List1 = data.List





subset2 = c("ROI_ID_2", "ROICoordinateX", "ROICoordinateY", "RawReads", "AlignedReads", "DeduplicatedReads", "TrimmedReads", 
"StitchedReads", "LDH_high_gt_240_at_beginig_of_treatment", "Performance_Status_at_treatment", "OS_days",
"LTB12", "LTB24", "DCB6", "primary_lymph_MET", "PFS_days", "Age_at_the_time_of_biopsy", "Sex",  "Prior_immunetherapy",
"Best_Overall_response", "Death_due_to_other_cause", "Melanoma_Type_Uveal_Cutaneous_mucosal", "Mutation_Status_BRAF_CKIK_NRAS", "PFS_index__1_Progression_",
"NormalizationFactor")
stage.col2 = c("NA",  "A", "B", "C", "D")


data.List = list(NULL)
clin.List = list(NULL)
for(k in 1:length(myinf2))
{
	info = read_excel(myinf2[k], sheet=1)
	info = as.data.frame(info)
	row.names(info) = gsub(" ", "", info$SegmentDisplayName)
	if(k!=3)
	{
		se = which(colnames(info)%in%stage.col2)
		xx = info[, se]
		for(i in 1:ncol(xx))
		{
			xx[,i] = ifelse(xx[,i]=="True", 1, 0)
		}
		tmp = apply(xx, 1, which.max)
		tmp = colnames(xx)[tmp]
		stage= tmp
	}else
	{
		stage = rep(NA, nrow(info))
	}
	xx = subset2[subset2%in%colnames(info)]
	info = info[, xx]
	info$stage = stage
	clin.List[[k]] = info 
	
	data = read_excel(myinf2[k], sheet=4)
	data = as.data.frame(data)
	tmp = data[,1]
	data = data[, -1]
	row.names(data) = tmp
	colnames(data) = gsub(" ", "", colnames(data) )
	data.List[[k]] = data
}
for(k in 1:length(data.List))
{
	data = data.List[[k]]
	if(k!=3)
	{
		nn = clin.List[[k]]$DeduplicatedReads
	}else
	{
		nn = apply(data, 2, sum)
	}
	for(i in 1:ncol(data))
	{
		data[,i] = data[,i]*1e6/nn[i]
	}
	data.List[[k]] = data
}
clin.List2 = clin.List
data.List2 = data.List

save(data.List1, data.List2, clin.List1, clin.List2, file= myoutf1)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[2] test the predictive power of genes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.1] calculate the average gene expression across all ROIs to obtain patient expression profiles  --> associate with response 
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

load(myinf1)

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
	data = data.List[[k]]	
	xx = colnames(data)
	for(i in 1:length(xx))
	{
		xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
	}
	colnames(data) = xx
	mypp = unique(colnames(data))
	mymat = matrix(0, nrow(data), length(mypp))
	row.names(mymat) = row.names(data)
	colnames(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(colnames(data)==mypp[i])
		mymat[,i] = apply(data[,se, drop=F], 1, mean)		
	}
	data.List[[k]] = mymat
}

#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]

res.List = list(NULL)
for(k in 1:3)
{
	data = data.List[[k]]
	data = log2(data+1)
	se = which(colnames(data)%in%pos.sam)
	dat1 = data[,se]
	se = which(colnames(data)%in%neg.sam)
	dat2 = data[,se]
	myavg1 = apply(dat1, 1, mean, na.rm=T)
	myavg2 = apply(dat2, 1, mean, na.rm=T)
	myvar1 = apply(dat1, 1, var, na.rm=T)
	myvar2 = apply(dat2, 1, var, na.rm=T)
	n1 = apply(!is.na(dat1), 1, sum)
	n2 = apply(!is.na(dat2), 1, sum)
	tscore = (myavg1-myavg2)/sqrt(myvar1/n1+myvar2/n2)
	df= (myvar1/n1+myvar2/n2)^2/((myvar1/n1)^2/(n1-1) + (myvar2/n2)^2/(n2-1))
	tscore[is.na(tscore)] = 0
	pval = pt(-abs(tscore), df)*2
	diff = myavg1-myavg2
	qval = p.adjust(pval, method="BH")
	res = data.frame(myavg1, myavg2, diff, tscore, pval, qval)
	row.names(res) = row.names(data)
	res = res[order(res$pval), ]
	res.List[[k]] = res
}
names(res.List) = c("CD45", "CD68", "S100")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## save feature for models
myoutf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P1_MHC_genes.txt"


mygen = NULL
for(k in 1:3)
{
  tmp = res.List[[k]]
  se = which(tmp$pval<0.01)
  tmp = tmp[se,]
  mygen = unique(c(mygen, row.names(tmp)))
}
mygen = mygen[grep("HLA-|B2M|CD74", mygen)]


for(k in 1:3)
{
  mygen = mygen[mygen%in%row.names(data.List[[k]])]
}


for(k in 1:3)
{
  tmp = t(data.List[[k]][mygen,])
  row.names(tmp) = paste(names(res.List)[k], row.names(tmp), sep="__")
  if(k==1)
  {
    mod.dat = tmp
  }else
  {
    mod.dat = rbind(mod.dat, tmp)
  }
  
  
}
mod.dat = log10(mod.dat)


write.table(mod.dat, myoutf1, sep="\t", quote=F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[3] immune cell marker genes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.1] Antigen-related pathways
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GSEA/c2.all.v2023.2.Hs.symbols.csv"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

tmpf = "~/TEMP/tmp.txt"

gset = read.table(myinf1, sep=",", header=T, row.names=1, quote="", check.names=F)
se = grep("REACTOME_|KEGG_", colnames(gset))
gset = gset[,se]
se = grep("ANTIGEN", colnames(gset))
reg = gset[,se]
xx = apply(reg, 2, sum)
reg = reg[, which(xx>=20)]


## check overlapping among pathways
xx = apply(reg,1,sum)
tmp <- reg[which(xx>0),]
tmp_MHCI <- tmp[, 1:7]
xx = apply(tmp_MHCI,1,sum)
tmp_MHCI <- tmp_MHCI[which(xx>0),] # 631 -> 535

tmp_MHCI <- tmp_MHCI[which(xx>1),] # 631 -> 381
## 30 out of 126 genes in REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION appear in other pathways

mark.List = list(NULL)
for(k in 1:ncol(reg))
{
  mark.List[[k]] = row.names(reg)[reg[,k]>0]
}
names(mark.List) = colnames(reg)
imm.gen = unique(unlist(mark.List))

#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1
sapply(data.List,  nrow)
sum(imm.gen%in% row.names(data.List[[1]]))
sum(imm.gen%in% row.names(data.List[[2]]))
sum(imm.gen%in% row.names(data.List[[3]]))

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
  data = data.List[[k]]	
  xx = colnames(data)
  for(i in 1:length(xx))
  {
    xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
  }
  colnames(data) = xx
  mypp = unique(colnames(data))
  mymat = matrix(0, nrow(data), length(mypp))
  row.names(mymat) = row.names(data)
  colnames(mymat) = mypp
  for(i in 1:length(mypp))
  {
    se = which(colnames(data)==mypp[i])
    mymat[,i] = apply(data[,se, drop=F], 1, mean)		
  }
  data.List[[k]] = mymat
}

#-----------------------
for(k in 1:length(clin.List))
{
  tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
  if(k==1)
  {
    info = tmp
  }else
  {
    info = rbind(info, tmp)
  }
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
  xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#--------------
source("~/WorSpa/system/myRprogram/base5.R")
##
for(k in 1:3)
{
  data = data.List[[k]]
  data = log2(data+1)
  reg = matrix(0, nrow(data), length(mark.List))
  row.names(reg) = row.names(data)
  colnames(reg) = names(mark.List)
  for(i in 1:length(mark.List))
  {
    se = which(row.names(reg)%in%mark.List[[i]])
    reg[se, i] = 1
  }
  xx = apply(reg, 2, sum)
  se = which(xx>=5)
  reg = reg[, se]
  
  xx = base5(data, reg, perm=1000, tmpf, median.norm=T)
  score = xx[[1]]
  colnames(score) = colnames(reg)
  row.names(score) = colnames(data)
  if(k==1) {dat.CD45 = score}
  if(k==2) {dat.CD68 = score}
  if(k==3) {dat.S100 = score}
}
##---------------

#---------------------
res.List = list(NULL)
for(Opt in 1:3)
{
  if(Opt==1) {data = dat.CD45}
  if(Opt==2) {data = dat.CD68}
  if(Opt==3) {data = dat.S100}
  
  se = which(row.names(data)%in%pos.sam)
  dat1 = data[se, ]
  se = which(row.names(data)%in%neg.sam)
  dat2 = data[se, ]
  
  res = matrix(0, ncol(dat1), 5)
  row.names(res) = colnames(dat1)
  colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
  res = as.data.frame(res)
  res[,1] = apply(dat1, 2, mean)
  res[,2] = apply(dat2, 2, mean)
  for(k in 1:ncol(dat1))
  {
    x1 = dat1[,k]
    x2 = dat2[,k]
    tmp = t.test(x1, x2)
    res[k,3] = tmp$statistic
    res[k,4] = tmp$p.value
    tmp = wilcox.test(x1, x2)
    res[k,5] = tmp$p.value
  }
  res.List[[Opt]] = res[order(res$p.w), 3:5]
}
names(res.List) = c("CD45", "CD68", "S100")
res.List


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.2] MHC -- GO
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GSEA/c5.all.v2023.2.Hs.symbols.csv"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

tmpf = "/mount/ictr1/chenglab/jdong/2lab/TEMP/tmp.txt"

gset = read.table(myinf1, sep=",", header=T, row.names=1, quote="", check.names=F)
se = grep("MHC", colnames(gset))
reg = gset[,se]

xx = apply(reg, 2, sum)
reg = reg[, which(xx>=10)]

mark.List = list(NULL)
for(k in 1:ncol(reg))
{
  mark.List[[k]] = row.names(reg)[reg[,k]>0]
}
names(mark.List) = colnames(reg)
imm.gen = unique(unlist(mark.List))

#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1
sapply(data.List,  nrow)
sum(imm.gen%in% row.names(data.List[[1]]))
sum(imm.gen%in% row.names(data.List[[2]]))
sum(imm.gen%in% row.names(data.List[[3]]))

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
  data = data.List[[k]]	
  xx = colnames(data)
  for(i in 1:length(xx))
  {
    xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
  }
  colnames(data) = xx
  mypp = unique(colnames(data))
  mymat = matrix(0, nrow(data), length(mypp))
  row.names(mymat) = row.names(data)
  colnames(mymat) = mypp
  for(i in 1:length(mypp))
  {
    se = which(colnames(data)==mypp[i])
    mymat[,i] = apply(data[,se, drop=F], 1, mean)		
  }
  data.List[[k]] = mymat
}

#-----------------------
for(k in 1:length(clin.List))
{
  tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
  if(k==1)
  {
    info = tmp
  }else
  {
    info = rbind(info, tmp)
  }
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
  xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#--------------
source("/mount/ictr1/chenglab/jdong/Tools/base5.R")
##
for(k in 1:3)
{
  data = data.List[[k]]
  data = log2(data+1)
  reg = matrix(0, nrow(data), length(mark.List))
  row.names(reg) = row.names(data)
  colnames(reg) = names(mark.List)
  for(i in 1:length(mark.List))
  {
    se = which(row.names(reg)%in%mark.List[[i]])
    reg[se, i] = 1
  }
  xx = apply(reg, 2, sum)
  se = which(xx>=5)
  reg = reg[, se]
  
  xx = base5(data, reg, perm=1000, tmpf, median.norm=T)
  score = xx[[1]]
  colnames(score) = colnames(reg)
  row.names(score) = colnames(data)
  if(k==1) {dat.CD45 = score}
  if(k==2) {dat.CD68 = score}
  if(k==3) {dat.S100 = score}
}
##---------------

#---------------------
res.List = list(NULL)
for(Opt in 1:3)
{
  if(Opt==1) {data = dat.CD45}
  if(Opt==2) {data = dat.CD68}
  if(Opt==3) {data = dat.S100}
  
  se = which(row.names(data)%in%pos.sam)
  dat1 = data[se, ]
  se = which(row.names(data)%in%neg.sam)
  dat2 = data[se, ]
  
  res = matrix(0, ncol(dat1), 5)
  row.names(res) = colnames(dat1)
  colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
  res = as.data.frame(res)
  res[,1] = apply(dat1, 2, mean)
  res[,2] = apply(dat2, 2, mean)
  for(k in 1:ncol(dat1))
  {
    x1 = dat1[,k]
    x2 = dat2[,k]
    tmp = t.test(x1, x2)
    res[k,3] = tmp$statistic
    res[k,4] = tmp$p.value
    tmp = wilcox.test(x1, x2)
    res[k,5] = tmp$p.value
  }
  res.List[[Opt]] = res[order(res$p.w), 3:5]
}
names(res.List) = c("CD45", "CD68", "S100")
res.List

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Save feature for models
myoutf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P2_MHCI_II.txt"


for(Opt in 1:3)
{
  if(Opt==1) {data = dat.CD45}
  if(Opt==2) {data = dat.CD68}
  if(Opt==3) {data = dat.S100}
  
  subset = c("GOCC_MHC_CLASS_I_PROTEIN_COMPLEX", "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX")
  tmp = data[,subset]
  colnames(tmp) = c("GOCC.MHC1", "GOCC.MHC2")
  if(Opt==1)
  {
    row.names(tmp) = paste("CD45__", row.names(tmp), sep="")
    mod.dat = tmp
  }
  if(Opt==2)
  {
    row.names(tmp) = paste("CD68__", row.names(tmp), sep="")
    mod.dat = rbind(mod.dat, tmp)
  }
  if(Opt==3)
  {
    row.names(tmp) = paste("S100__", row.names(tmp), sep="")
    mod.dat = rbind(mod.dat, tmp)
  }
}


write.table(mod.dat, myoutf1, sep="\t", quote=F)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[3.3] Panglaodb
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Database/Panglaodb/PanglaoDB_markers_27_Mar_2020.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

tmpf = "~/TEMP/tmp.txt"


data = read.table(myinf1,  sep="\t", header=T, quote="")

se = which(data$organ=="Immune system")
data = data[se,]
se = which(data$species %in% c("Mm Hs", "Hs"))
data = data[se,]
xx = table(data$cell.type)
xx = xx[xx>50]
mycel = names(xx)
se = which(data$cell.type %in% mycel)
data = data[se,]
mark.List = list(NULL)
for(k in 1:length(mycel))
{
	mark.List[[k]] = unique(data$official.gene.symbol[data$cell.type == mycel[k]])
}
names(mark.List) = mycel
imm.gen = unique(unlist(mark.List))

#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1
sapply(data.List,  nrow)
sum(imm.gen%in% row.names(data.List[[1]]))
sum(imm.gen%in% row.names(data.List[[2]]))
sum(imm.gen%in% row.names(data.List[[3]]))

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
	data = data.List[[k]]	
	xx = colnames(data)
	for(i in 1:length(xx))
	{
		xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
	}
	colnames(data) = xx
	mypp = unique(colnames(data))
	mymat = matrix(0, nrow(data), length(mypp))
	row.names(mymat) = row.names(data)
	colnames(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(colnames(data)==mypp[i])
		mymat[,i] = apply(data[,se, drop=F], 1, mean)		
	}
	data.List[[k]] = mymat
}

#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#--------------
source("~/WorSpa/system/myRprogram/base5.R")
##
for(k in 1:3)
{
	data = data.List[[k]]
	data = log2(data+1)
	reg = matrix(0, nrow(data), length(mark.List))
	row.names(reg) = row.names(data)
	colnames(reg) = names(mark.List)
	for(i in 1:length(mark.List))
	{
		se = which(row.names(reg)%in%mark.List[[i]])
		reg[se, i] = 1
	}
	xx = apply(reg, 2, sum)
	se = which(xx>=10)
	reg = reg[, se]
	
	xx = base5(data, reg, perm=1000, tmpf, median.norm=T)
	score = xx[[1]]
	colnames(score) = colnames(reg)
	row.names(score) = colnames(data)
	if(k==1) {dat.CD45 = score}
	if(k==2) {dat.CD68 = score}
	if(k==3) {dat.S100 = score}
}
##---------------

#---------------------
res.List = list(NULL)
for(Opt in 1:3)
{
	if(Opt==1) {data = dat.CD45}
	if(Opt==2) {data = dat.CD68}
	if(Opt==3) {data = dat.S100}

	se = which(row.names(data)%in%pos.sam)
	dat1 = data[se, ]
	se = which(row.names(data)%in%neg.sam)
	dat2 = data[se, ]

	res = matrix(0, ncol(dat1), 5)
	row.names(res) = colnames(dat1)
	colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
	res = as.data.frame(res)
	res[,1] = apply(dat1, 2, mean)
	res[,2] = apply(dat2, 2, mean)
	for(k in 1:ncol(dat1))
	{
		x1 = dat1[,k]
		x2 = dat2[,k]
		tmp = t.test(x1, x2)
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(x1, x2)
		res[k,5] = tmp$p.value
	}
	res.List[[Opt]] = res[order(res$p.w), 3:5]
}
names(res.List) = c("CD45", "CD68", "S100")
res.List


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.4] Charoentong2017
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/ImmGen/Charoentong2017_CellReport_TableS6_ImmCelSpeGenes.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

tmpf = "~/TEMP/tmp.txt"


data = read.table(myinf1,  sep="\t", header=T, quote="")
sort(table(data$Cell.type))
mycel = unique(data$Cell.type)

mark.List = list(NULL)
for(k in 1:length(mycel))
{
	mark.List[[k]] = unique(data$Metagene[data$Cell.type == mycel[k]])
}
names(mark.List) = mycel
imm.gen = unique(unlist(mark.List))

#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1
sapply(data.List,  nrow)
sum(imm.gen%in% row.names(data.List[[1]]))
sum(imm.gen%in% row.names(data.List[[2]]))
sum(imm.gen%in% row.names(data.List[[3]]))

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
	data = data.List[[k]]	
	xx = colnames(data)
	for(i in 1:length(xx))
	{
		xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
	}
	colnames(data) = xx
	mypp = unique(colnames(data))
	mymat = matrix(0, nrow(data), length(mypp))
	row.names(mymat) = row.names(data)
	colnames(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(colnames(data)==mypp[i])
		mymat[,i] = apply(data[,se, drop=F], 1, mean)		
	}
	data.List[[k]] = mymat
}

#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#--------------
source("~/WorSpa/system/myRprogram/base5.R")
##
for(k in 1:3)
{
	data = data.List[[k]]
	data = log2(data+1)
	reg = matrix(0, nrow(data), length(mark.List))
	row.names(reg) = row.names(data)
	colnames(reg) = names(mark.List)
	for(i in 1:length(mark.List))
	{
		se = which(row.names(reg)%in%mark.List[[i]])
		reg[se, i] = 1
	}
	xx = apply(reg, 2, sum)
	se = which(xx>=5)
	reg = reg[, se]
	
	xx = base5(data, reg, perm=1000, tmpf, median.norm=T)
	score = xx[[1]]
	colnames(score) = colnames(reg)
	row.names(score) = colnames(data)
	if(k==1) {dat.CD45 = score}
	if(k==2) {dat.CD68 = score}
	if(k==3) {dat.S100 = score}
}
##---------------

#---------------------
res.List = list(NULL)
for(Opt in 1:3)
{
	if(Opt==1) {data = dat.CD45}
	if(Opt==2) {data = dat.CD68}
	if(Opt==3) {data = dat.S100}

	se = which(row.names(data)%in%pos.sam)
	dat1 = data[se, ]
	se = which(row.names(data)%in%neg.sam)
	dat2 = data[se, ]

	res = matrix(0, ncol(dat1), 5)
	row.names(res) = colnames(dat1)
	colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
	res = as.data.frame(res)
	res[,1] = apply(dat1, 2, mean)
	res[,2] = apply(dat2, 2, mean)
	for(k in 1:ncol(dat1))
	{
		x1 = dat1[,k]
		x2 = dat2[,k]
		tmp = t.test(x1, x2)
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(x1, x2)
		res[k,5] = tmp$p.value
	}
	res.List[[Opt]] = res[order(res$p.w), 3:5]
}
names(res.List) = c("CD45", "CD68", "S100")
res.List

#-------------
xx = mark.List[[which(names(mark.List)=="Effector memeory CD8 T cell")]]
xx[xx%in%row.names(data.List[[1]])]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.5] Charoentong2017: calculate region-specific score --> convert into patient-level by avg/max/min
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/ImmGen/Charoentong2017_CellReport_TableS6_ImmCelSpeGenes.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

tmpf = "/mount/ictr1/chenglab/jdong/2lab/TEMP/tmp.txt"


data = read.table(myinf1,  sep="\t", header=T, quote="")
sort(table(data$Cell.type))
mycel = unique(data$Cell.type)

mark.List = list(NULL)
for(k in 1:length(mycel))
{
	mark.List[[k]] = unique(data$Metagene[data$Cell.type == mycel[k]])
}
names(mark.List) = mycel
imm.gen = unique(unlist(mark.List))

#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1
sapply(data.List,  nrow)
sum(imm.gen%in% row.names(data.List[[1]]))
sum(imm.gen%in% row.names(data.List[[2]]))
sum(imm.gen%in% row.names(data.List[[3]]))

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#--------------
source("/mount/ictr1/chenglab/jdong/Tools/base5.R")
##
for(k in 1:3)
{
	data = data.List[[k]]
	data = log2(data+1)
	reg = matrix(0, nrow(data), length(mark.List))
	row.names(reg) = row.names(data)
	colnames(reg) = names(mark.List)
	for(i in 1:length(mark.List))
	{
		se = which(row.names(reg)%in%mark.List[[i]])
		reg[se, i] = 1
	}
	xx = apply(reg, 2, sum)
	se = which(xx>=5)
	reg = reg[, se]
	
	xx = base5(data, reg, perm=1000, tmpf, median.norm=T)
	score = xx[[1]]
	colnames(score) = colnames(reg)
	row.names(score) = colnames(data)
	if(k==1) {dat.CD45 = score}
	if(k==2) {dat.CD68 = score}
	if(k==3) {dat.S100 = score}
}
##---------------

#---------------------
res.List = list(NULL)
for(Opt in 1:3)
{
	if(Opt==1) {data = dat.CD45}
	if(Opt==2) {data = dat.CD68}
	if(Opt==3) {data = dat.S100}
	
	nn = nrow(data)
	mytag = unlist(strsplit(row.names(data), "\\|"))[(1:nn)*3-1]
	mypp = unique(mytag)
	mymat = matrix(0, length(mypp), ncol(data))
	colnames(mymat) = colnames(data)
	row.names(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(mytag==mypp[i])
		mymat[i,] = apply(data[se, , drop=F], 2, min)	## mean/max/min	
	}

	se = which(row.names(mymat)%in%pos.sam)
	dat1 = mymat[se, ]
	se = which(row.names(mymat)%in%neg.sam)
	dat2 = mymat[se, ]

	res = matrix(0, ncol(dat1), 5)
	row.names(res) = colnames(dat1)
	colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
	res = as.data.frame(res)
	res[,1] = apply(dat1, 2, mean)
	res[,2] = apply(dat2, 2, mean)
	for(k in 1:ncol(dat1))
	{
		x1 = dat1[,k]
		x2 = dat2[,k]
		tmp = t.test(x1, x2)
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(x1, x2)
		res[k,5] = tmp$p.value
	}
	res.List[[Opt]] = res[order(res$p.w), 3:5]
}
names(res.List) = c("CD45", "CD68", "S100")
res.List

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Save feature for models
myoutf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P3_ImmSigScore.txt"


for(Opt in 1:3)
{
  if(Opt==1) {data = dat.CD45}
  if(Opt==2) {data = dat.CD68}
  if(Opt==3) {data = dat.S100}
  
  
  nn = nrow(data)
  mytag = unlist(strsplit(row.names(data), "\\|"))[(1:nn)*3-1]
  mypp = unique(mytag)
  mymat = matrix(0, length(mypp), ncol(data))
  colnames(mymat) = colnames(data)
  row.names(mymat) = mypp
  for(i in 1:length(mypp))
  {
    se = which(mytag==mypp[i])
    mymat[i,] = apply(data[se, , drop=F], 2, min)	## mean/max/min	
  }
  se = which(colnames(mymat)=="Activated dendritic cell")
  aDC.min = mymat[,se]
  se = which(colnames(mymat)=="Effector memeory CD8 T cell")
  emCD8T.min = mymat[,se]
  
  for(i in 1:length(mypp))
  {
    se = which(mytag==mypp[i])
    mymat[i,] = apply(data[se, , drop=F], 2, max)	## mean/max/min	
  }
  se = which(colnames(mymat)=="Activated dendritic cell")
  aDC.max = mymat[,se]
  se = which(colnames(mymat)=="Effector memeory CD8 T cell")
  emCD8T.max = mymat[,se]
  
  for(i in 1:length(mypp))
  {
    se = which(mytag==mypp[i])
    mymat[i,] = apply(data[se, , drop=F], 2, mean)	## mean/max/min	
  }
  se = which(colnames(mymat)=="Activated dendritic cell")
  aDC.avg = mymat[,se]
  se = which(colnames(mymat)=="Effector memeory CD8 T cell")
  emCD8T.avg = mymat[,se]
  
  if(Opt==1)
  {
    tmp = data.frame(aDC.avg, aDC.max, aDC.min, emCD8T.avg, emCD8T.max, emCD8T.min)
    row.names(tmp) = paste("CD45__", mypp, sep="")
    mod.dat = tmp
  }
  if(Opt==2)
  {
    tmp = data.frame(aDC.avg, aDC.max, aDC.min, emCD8T.avg, emCD8T.max, emCD8T.min)
    row.names(tmp) = paste("CD68__", mypp, sep="")
    mod.dat = rbind(mod.dat, tmp)
  }
  if(Opt==3)
  {
    tmp = data.frame(aDC.avg, aDC.max, aDC.min, emCD8T.avg, emCD8T.max, emCD8T.min)
    row.names(tmp) = paste("S100__", mypp, sep="")
    mod.dat = rbind(mod.dat, tmp)
  }
}


write.table(mod.dat, myoutf1, sep="\t", quote=F)



#-------------------------------
## combine different regions
myList = list(CD45=dat.CD45, CD68 = dat.CD68,S100=dat.S100)
sapply(myList, nrow)
for(k in 1:length(myList))
{
	data = myList[[k]]
	nn = nrow(data)
	mytag = unlist(strsplit(row.names(data), "\\|"))[(1:nn)*3-1]
	mypp = unique(mytag)
	mymat = matrix(0, length(mypp), ncol(data))
	colnames(mymat) = colnames(data)
	row.names(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(mytag==mypp[i])
		mymat[i,] = apply(data[se, , drop=F], 2, mean)	## mean/max/min	
	}
	myList[[k]] = mymat
}
sapply(myList, nrow)
		CD45 CD68 S100 
		  46   65   88 
		CD45 CD68 S100 
		  37   48   53 
#--------------	
datx = myList[[3]]
daty = myList[[1]]
comxx = intersect(row.names(datx), row.names(daty))
datx = datx[comxx,]
daty = daty[comxx,]
comxx = intersect(colnames(datx), colnames(daty))
datx = datx[,comxx]
daty = daty[,comxx]
mymat = datx - daty

se = which(row.names(mymat)%in%pos.sam)
dat1 = mymat[se, ]
se = which(row.names(mymat)%in%neg.sam)
dat2 = mymat[se, ]

res = matrix(0, ncol(dat1), 5)
row.names(res) = colnames(dat1)
colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
res = as.data.frame(res)
res[,1] = apply(dat1, 2, mean)
res[,2] = apply(dat2, 2, mean)
for(k in 1:ncol(dat1))
{
	x1 = dat1[,k]
	x2 = dat2[,k]
	tmp = t.test(x1, x2)
	res[k,3] = tmp$statistic
	res[k,4] = tmp$p.value
	tmp = wilcox.test(x1, x2)
	res[k,5] = tmp$p.value
}
res = res[order(res$p.w), ]
res




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4] deconvolution
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.1] patient-specific analysis for S100
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/annotation/UCSC_hg19_Gene_Length.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

myoutf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_PatAvg_deconvolution.rda"

library("immunedeconv")

data = read.table(myinf1,  sep="\t", header=T, quote="")
mylen = data$avg.tr.len
names(mylen) = row.names(data)
#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1
for(k in 1:length(data.List))
{
	data = data.List[[k]]	
	xx = colnames(data)
	for(i in 1:length(xx))
	{
		xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
	}
	colnames(data) = xx
	mypp = unique(colnames(data))
	mymat = matrix(0, nrow(data), length(mypp))
	row.names(mymat) = row.names(data)
	colnames(mymat) = mypp
	for(i in 1:length(mypp))
	{
		se = which(colnames(data)==mypp[i])
		mymat[,i] = apply(data[,se, drop=F], 1, mean)		
	}
	data.List[[k]] = mymat
}

for(k in 1:length(data.List))
{
	data = data.List[[k]]
	xx = mylen[row.names(data)]
	se = which(!is.na(xx))
	data = data[se,]
	xx = xx[se]
	for(i in 1:ncol(data))
	{
		tmp = data[,i]/xx
		data[,i] = tmp*1e6/sum(tmp)
	}
	data.List[[k]] = data
}
sapply(data.List, nrow)

myList = list(NULL)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Deconvolution
df = data.List[[3]]
dim(df)

# Method-1 TIMER2.0
cancertype <- "skcm"
immunedeconv::timer_available_cancers
myres = as.data.frame(deconvolute(df, "timer",indications=rep(tolower(cancertype),ncol(df))))
xx = as.vector(myres[,1])
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

myList[[1]] = res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method-2 Cibersort
set_cibersort_binary("~/WorSpa/system/myRprogram/CIBERSORT.R")
set_cibersort_mat("~/WorSpa/system/myRprogram/LM22.txt")
myres = as.data.frame(deconvolute(df, "cibersort"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

myList[[2]] = res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method-3 Run quanTIseq
myres = as.data.frame(deconvolute(df, "quantiseq"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

myList[[3]] = res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method-4 xCell
#Run xCell
myres = as.data.frame(deconvolute(df, "xcell"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

myList[[4]] = res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method-5 #Run EPIC
myres = as.data.frame(deconvolute(df, "epic"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

myList[[5]] = res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method-6 #Run MCP-counter
myres = as.data.frame(deconvolute(df, "mcp_counter"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx
myList[[6]] = res

names(myList) = c("timer", "cibersort", "quantiseq", "xcell", "epic", "mcp_counter")
S100.List = myList


save(S100.List,  file=myoutf1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.2] Evaluate the signficance
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_PatAvg_deconvolution.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"

load(myinf1)
## S100.List
#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1

#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]


#-----------------------
myList = S100.List

res.List = list(NULL)
for(p in 1:length(myList))
{
	data = myList[[p]]
	xx = apply(data!=0, 2, sum)
	data = data[, which(xx>=10)]
	se = which(row.names(data)%in%pos.sam)
	dat1 = data[se, ]
	se = which(row.names(data)%in%neg.sam)
	dat2 = data[se, ]

	res = matrix(0, ncol(dat1), 5)
	row.names(res) = colnames(dat1)
	colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
	res = as.data.frame(res)
	res[,1] = apply(dat1, 2, mean)
	res[,2] = apply(dat2, 2, mean)
	for(k in 1:ncol(dat1))
	{
		x1 = dat1[,k]
		x2 = dat2[,k]
		tmp = t.test(x1, x2)
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(x1, x2)
		res[k,5] = tmp$p.value
	}
	res.List[[p]] = res[order(res$p.w), 3:5]
}
names(res.List) = names(myList)
res.List



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.3] region-specific analysis for S100
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/annotation/UCSC_hg19_Gene_Length.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"


library("immunedeconv")

data = read.table(myinf1,  sep="\t", header=T, quote="")
mylen = data$avg.tr.len
names(mylen) = row.names(data)
#--------------------
load(myinf2)
data.List = data.List1
clin.List = clin.List1

#------------------
## transform ROI expression profiles into patient expression profiles
data.List = data.List1
clin.List = clin.List1

for(k in 1:length(data.List))
{
	data = data.List[[k]]
	xx = mylen[row.names(data)]
	se = which(!is.na(xx))
	data = data[se,]
	xx = xx[se]
	for(i in 1:ncol(data))
	{
		tmp = data[,i]/xx
		data[,i] = tmp*1e6/sum(tmp)
	}
	data.List[[k]] = data
}
sapply(data.List, nrow)

#-----------------------
for(k in 1:length(clin.List))
{
	tmp = unique(clin.List[[k]][, c("ROI_ID_2", "response", "AGE_AT_DX")])
	if(k==1)
	{
		info = tmp
	}else
	{
		info = rbind(info, tmp)
	}
}
info = unique(info)
xx = row.names(info)
for(i in 1:length(xx))
{
	xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
row.names(info) = xx

pos.sam = row.names(info)[info$response=="yes"]
neg.sam = row.names(info)[info$response=="no"]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Deconvolution
# Method-2 Cibersort
df = data.List[[3]]
set_cibersort_binary("~/WorSpa/system/myRprogram/CIBERSORT.R")
set_cibersort_mat("~/WorSpa/system/myRprogram/LM22.txt")
myres = as.data.frame(deconvolute(df, "cibersort"))
xx = myres[,1]
res = t(myres[, 2:ncol(myres)])
colnames(res) = xx

data = res
xx = row.names(data)
nn = length(xx)
mytag = unlist(strsplit(xx, "\\|"))[(1:nn)*3-1]
mypat = unique(mytag)

tmp = matrix(0, length(mypat), ncol(data))
row.names(tmp) = mypat
colnames(tmp) = colnames(data)
tmp = as.data.frame(tmp)
dat.max = dat.min = dat.avg = tmp
for(k in 1:length(mypat))
{
	se = which(mytag==mypat[k])
	tmp = data[se,,drop=F]
	dat.max[k,] = as.numeric(apply(tmp,2 , max))	
	dat.min[k,] = as.numeric(apply(tmp,2 , min))	
	dat.avg[k,] = as.numeric(apply(tmp,2 , mean))	
}
myList = list(max=dat.max, min = dat.min, avg=dat.avg)

#-----------------------
res.List = list(NULL)
for(p in 1:length(myList))
{
	data = myList[[p]]
	xx = apply(data!=0, 2, sum)
	data = data[, which(xx>=10)]
	se = which(row.names(data)%in%pos.sam)
	dat1 = data[se, ]
	se = which(row.names(data)%in%neg.sam)
	dat2 = data[se, ]

	res = matrix(0, ncol(dat1), 5)
	row.names(res) = colnames(dat1)
	colnames(res) = c("avgT", "avgN", "tscore", "p.t", "p.w")
	res = as.data.frame(res)
	res[,1] = apply(dat1, 2, mean)
	res[,2] = apply(dat2, 2, mean)
	for(k in 1:ncol(dat1))
	{
		x1 = dat1[,k]
		x2 = dat2[,k]
		tmp = t.test(x1, x2)
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(x1, x2)
		res[k,5] = tmp$p.value
	}
	res.List[[p]] = res[order(res$p.w), 3:5]
}
names(res.List) = names(myList)
res.List



## [5] Models -- RF
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.1] Collect features
## HLA genes
## GOCC MHC-I and MHC-II
## DC Signatures and CD8T signatures (avg/max/min)
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P1_MHC_genes.txt"
myinf2 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P2_MHCI_II.txt"
myinf3 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_P3_ImmSigScore.txt"


myoutf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_AllFeatures.txt"


dat1 = read.table(myinf1, sep="\t", header=T, row.names=1)
dat2 = read.table(myinf2, sep="\t", header=T, row.names=1)
dat3 = read.table(myinf3, sep="\t", header=T, row.names=1)
dim(dat1)
dim(dat2)
dim(dat3)


res = cbind(dat1, dat2, dat3)
write.table(res, myoutf1, sep="\t", quote=F)


#+++++++++++++++++++++++++++++++++++++++++
[5.2] play with the data
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_AllFeatures.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"


load(myinf2)
clin.List = clin.List1


data = read.table(myinf1, sep="\t", header=T, row.names=1)
raw.data = data


subset = c("PID", "response")
res.List = list(NULL)
for(k in 1:3)
{
  data = raw.data
  info = clin.List[[k]]
  xx = paste(names(clin.List)[k], "__", sep="") 
  se = grep(xx, row.names(data))
  data = data[se,]
  row.names(data) = gsub(xx, "", row.names(data))
  
  PID = row.names(info)
  for(i in 1:length(PID))
  {
    PID[i] = unlist(strsplit(PID[i], "\\|"))[2]
  }
  info = cbind(PID, info)
  info = unique(info[,subset])
  row.names(info) = info$PID
  comxx = intersect(row.names(data), row.names(info))
  data = data[comxx,]
  info = info[comxx,]
  
  se = which(info$response == "yes")
  dat1 = data[se,]
  se = which(info$response == "no")
  dat2 = data[se,]
  
  res = matrix(0, ncol(data), 3)
  row.names(res) = colnames(data)
  colnames(res) = c("tscore", "P.t", "P.w")
  res = as.data.frame(res)
  for(i in 1:ncol(data))
  {
    xx1 = as.numeric(dat1[,i])
    xx2 = as.numeric(dat2[,i])
    tmp = t.test(xx1, xx2)
    res[i,1] = tmp$statistic
    res[i,2] = tmp$p.value
    tmp = wilcox.test(xx1, xx2)
    res[i,3] = tmp$p.value		
  }
  res.List[[k]] = res
}
names(res.List) = names(clin.List)
res.List


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.3] RF model
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_AllFeatures.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"


library(randomForest)
library(ROCR)


load(myinf2)
clin.List = clin.List1


data = read.table(myinf1, sep="\t", header=T, row.names=1)
raw.data = data


subset = c("PID", "response", "AGE_AT_DX", "SEX", "RACE", "OS_FROM_START_OF_ITX", "VITAL", "PROG_STATUS_BY_SCAN")


Opt =0
if(Opt==0){	se = 1:ncol(raw.data)	}
if(Opt==1){	se = c(1:7, grep("avg", colnames(raw.data)))	}
if(Opt==2){	se = c(1:7, grep("max", colnames(raw.data)))	}
if(Opt==3){	se = c(1:7, grep("min", colnames(raw.data)))	}
sub.raw.data = raw.data[,se]


res.List = list(NULL)
myauc = rep(0, 3)
names(myauc) = names(clin.List)
for(k in 1:3)
{
  data = sub.raw.data
  xx = paste(names(clin.List)[k], "__", sep="") 
  se = grep(xx, row.names(data))
  data = data[se,]
  row.names(data) = gsub(xx, "", row.names(data))
  
  info = clin.List[[k]]
  PID = row.names(info)
  for(i in 1:length(PID))
  {
    PID[i] = unlist(strsplit(PID[i], "\\|"))[2]
  }
  info = cbind(PID, info)
  info = unique(info[,subset])
  colnames(info)[3] = "AGE"
  colnames(info)[6:8] = c("OS.time", "OS.event", "PROG")
  row.names(info) = info$PID
  
  
  comxx = intersect(row.names(data), row.names(info))
  data = data[comxx,]
  info = info[comxx,]
  
  myy = as.factor(info$response)
  mydat = cbind(myy, data)
  score = rep(0, nrow(mydat))
  cat = mydat$myy
  
  for(i in 1:nrow(mydat))
  {
    tr = mydat[-i,]	
    te = mydat[i,,drop=F]
    fit <- randomForest(myy~., data=tr)
    tmp = predict(fit, te[,-1, drop=F], type="prob")
    score[i] = as.numeric(tmp[1, "yes"])
  }
  
  
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[k] = as.numeric(auc.perf@y.values)
  
  
  fit <- randomForest(myy~., importance=T, data=mydat)
  xx = as.data.frame(fit$importance)
  xx= xx[order(xx$MeanDecreaseGini, decreasing=T), ]
  res.List[[k]] = xx
}
names(res.List) = names(clin.List)
myauc
res.List


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.4] RF model  -- combine all components
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/jdong/2lab/Spatial_Signatures/GeoMX_MS/Discovery_Cohort/Aung_GSE233305_Model_Discovery_AllFeatures.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"


library(randomForest)
library(ROCR)


load(myinf2)
clin.List = clin.List1
xx = row.names(clin.List[[1]])
for(i in 1:length(xx))
{
  xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
xx1 = unique(xx)
xx = row.names(clin.List[[2]])
for(i in 1:length(xx))
{
  xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
xx2 = unique(xx)
xx = row.names(clin.List[[3]])
for(i in 1:length(xx))
{
  xx[i] = unlist(strsplit(xx[i], "\\|"))[2]
}
xx3 = unique(xx)
length(xx1)  ## 37
length(xx2)	## 48
length(xx3)	## 53


mypat = intersect(intersect(xx1, xx2), xx3)
length(mypat) ## 33


#------------------------
subset = c("PID", "response", "AGE_AT_DX", "SEX", "RACE", "OS_FROM_START_OF_ITX", "VITAL", "PROG_STATUS_BY_SCAN")
info = clin.List[[1]]
PID = row.names(info)
for(i in 1:length(PID))
{
  PID[i] = unlist(strsplit(PID[i], "\\|"))[2]
}
info = cbind(PID, info)
info = unique(info[,subset])
colnames(info)[3] = "AGE"
colnames(info)[6:8] = c("OS.time", "OS.event", "PROG")
row.names(info) = info$PID
se = which(row.names(info)%in%mypat)
info = info[se,]




#------------------------
data = read.table(myinf1, sep="\t", header=T, row.names=1)
se = grep("CD45", row.names(data))
tmp = data[se,]
row.names(tmp) = gsub("CD45__", "", row.names(tmp))
dat1 = tmp[mypat,]
se = grep("CD68", row.names(data))
tmp = data[se,]
row.names(tmp) = gsub("CD68__", "", row.names(tmp))
dat2 = tmp[mypat,]
se = grep("S100", row.names(data))
tmp = data[se,]
row.names(tmp) = gsub("S100__", "", row.names(tmp))
dat3 = tmp[mypat,]
colnames(dat1) = paste("CD45__", colnames(dat1), sep="")
colnames(dat2) = paste("CD68__", colnames(dat2), sep="")
colnames(dat3) = paste("S100__", colnames(dat3), sep="")
data = cbind(dat1, dat2, dat3)


comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

# try only significant features --> worse performance ---------------------------------------------------------------------------------
#colnames(data) 

#data = data[, c("CD45__GOCC.MHC2","CD45__emCD8T.avg" ,"CD45__emCD8T.max", "CD45__emCD8T.min", "CD68__CD74", "CD68__B2M", "CD68__HLA.DRA",
#                "CD68__HLA.DQB1",  "CD68__GOCC.MHC2", "CD68__aDC.avg", "CD68__aDC.max", "CD68__aDC.min", "S100__HLA.F", "S100__GOCC.MHC1",
#                "S100__aDC.avg", "S100__aDC.max", "S100__aDC.min" )]

raw.data = data


#------------------
res.List = list(NULL)
myauc = rep(0, 4)
for(Opt in 1:4)
{
  if(Opt==1){	se = 1:ncol(raw.data)	}
  if(Opt==2){	se = c(c(1, 4:8, 12:13), grep("avg", colnames(raw.data)))	}
  if(Opt==3){	se = c(c(1, 4:8, 12:13), grep("max", colnames(raw.data)))	}
  if(Opt==4){	se = c(c(1, 4:8, 12:13), grep("min", colnames(raw.data)))	}
  sub.raw.data = raw.data[,se]
  
  
  data = sub.raw.data
  myy = as.factor(info$response)
  mydat = cbind(myy, data)
  score = rep(0, nrow(mydat))
  cat = mydat$myy
  
  for(i in 1:nrow(mydat))
  {
    tr = mydat[-i,]	
    te = mydat[i,,drop=F]
    fit <- randomForest(myy~., data=tr)
    tmp = predict(fit, te[,-1, drop=F], type="prob")
    score[i] = as.numeric(tmp[1, "yes"])
  }
  
  
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[Opt] = as.numeric(auc.perf@y.values)
  
  
  fit <- randomForest(myy~., importance=T, data=mydat)
  xx = as.data.frame(fit$importance)
  xx= xx[order(xx$MeanDecreaseGini, decreasing=T), 3:4]
  res.List[[Opt]] = xx
}
names(myauc) = names(res.List)  = c("All", "avg", "max", "min")




#@@@@&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Check the variation between different ROIs in S100B from the same patients
## → high intratumoral heterogeneity (ITH) of melanoma

rm(list = ls())
myinf = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Melanoma/Aung_GSE233305_GeoMx/data/Aung_GSE233305_GeoMx_processed.rda"
load(myinf)
data.List = data.List1
clin.List = clin.List1

## S100B
S100B <- t(as.data.frame(data.List[["S100"]]))
S100B <- as.data.frame(S100B)

S100B$PatientID <- sapply(rownames(S100B), function(x) strsplit(x, "\\|")[[1]][2])
S100B$ROI <- ifelse(grepl("-1-", rownames(S100B)), "ROI1", "ROI2")

S100B_long <- S100B %>%
  pivot_longer(
    cols = -c(PatientID, ROI),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  pivot_wider(
    names_from = ROI,
    values_from = expression
  ) %>%
  mutate(
    RowID = paste(gene, PatientID, sep = "_")
  )

## filter patients with both ROIs
data <- S100B_long %>% filter(!is.na(S100B_long$ROI1) & !is.na(S100B_long$ROI2))

# correlation
correlations <- data %>%
  group_by(PatientID) %>%
  summarize(correlation = cor(ROI1, ROI2, use = "pairwise.complete.obs"))

# Differences & overall variability
var <- data %>%
  mutate(abs_diff = abs(ROI1 - ROI2)) %>%
  group_by(PatientID) %>%
  summarize(mean = mean(abs_diff, na.rm = TRUE),
            sd = sd(abs_diff, na.rm = TRUE))

data$diff <- abs(data$ROI1 - data$ROI2)

global_var <- S100B_long %>%
  group_by(gene) %>%
  summarize(mean = mean(diff, na.rm = TRUE),
            sd = sd(diff, na.rm = TRUE))

# Shannon Entropy & paired comparison
cal_entropy <- function(P) {
  -sum(P * log(P))
}

norm_to_prob <- function(x) {
  x / sum(x)
}

entropy_scores <- data %>%
  mutate(
    ROI1_Prob = norm_to_prob(ROI1),
    ROI2_Prob = norm_to_prob(ROI2)
  ) %>%
  group_by(PatientID) %>%
  summarize(
    ROI1_Entropy = cal_entropy(ROI1_Prob),
    ROI2_Entropy = cal_entropy(ROI2_Prob)
  )

res <- wilcox.test(entropy_scores$ROI1_Entropy, entropy_scores$ROI2_Entropy, paired = TRUE)

## wilcox test group by patient 
res_gene <- data %>%
  group_by(PatientID) %>%
  summarise(
    p_value = list(wilcox.test(ROI1, ROI2, paired = TRUE)$p.value)) %>%
  unnest(cols = p_value) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj)

res_gene$p_adj <- format(res_gene$p_adj, scientific = TRUE,digits = 3)
#

## CD45
CD45 <- t(as.data.frame(data.List[["CD45"]]))
CD45 <- as.data.frame(CD45)

CD45$PatientID <- sapply(rownames(CD45), function(x) strsplit(x, "\\|")[[1]][2])
CD45$ROI <- ifelse(grepl("-1-", rownames(CD45)), "ROI1", "ROI2")

CD45_long <- CD45 %>%
  pivot_longer(
    cols = -c(PatientID, ROI),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  pivot_wider(
    names_from = ROI,
    values_from = expression
  ) %>%
  mutate(
    RowID = paste(gene, PatientID, sep = "_")
  )

## filter patients with both ROIs
data <- CD45_long %>% filter(!is.na(CD45_long$ROI1) & !is.na(CD45_long$ROI2)) ## 9 patients

# correlation
correlations <- data %>%
  group_by(PatientID) %>%
  summarize(correlation = cor(ROI1, ROI2, use = "complete.obs"))

## wilcox test group by patient 
res_gene <- data %>%
  group_by(PatientID) %>%
  summarise(
    p_value = list(wilcox.test(ROI1, ROI2, paired = TRUE)$p.value)) %>%
  unnest(cols = p_value) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj)

res_gene$p_adj <- format(res_gene$p_adj, scientific = TRUE,digits = 3)
  
## CD68
CD68 <- t(as.data.frame(data.List[["CD68"]]))
CD68 <- as.data.frame(CD68)

CD68$PatientID <- sapply(rownames(CD68), function(x) strsplit(x, "\\|")[[1]][2])
CD68$ROI <- ifelse(grepl("-1-", rownames(CD68)), "ROI1", "ROI2")

CD68_long <- CD68 %>%
  pivot_longer(
    cols = -c(PatientID, ROI),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  pivot_wider(
    names_from = ROI,
    values_from = expression
  ) %>%
  mutate(
    RowID = paste(gene, PatientID, sep = "_")
  )

## filter patients with both ROIs
data <- CD68_long %>% filter(!is.na(CD68_long$ROI1) & !is.na(CD68_long$ROI2)) ## 17 patients

# correlation
correlations <- data %>%
  group_by(PatientID) %>%
  summarize(correlation = cor(ROI1, ROI2, use = "complete.obs"))

## wilcox test group by patient 
res_gene <- data %>%
  group_by(PatientID) %>%
  summarise(
    p_value = list(wilcox.test(ROI1, ROI2, paired = TRUE)$p.value)) %>%
  unnest(cols = p_value) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj)

res_gene$p_adj <- format(res_gene$p_adj, scientific = TRUE,digits = 3)
res_gene$p_value <- format(res_gene$p_value, scientific = TRUE,digits = 3)