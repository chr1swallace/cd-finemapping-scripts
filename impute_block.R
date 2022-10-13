#!/usr/bin/env Rscript

library(randomFunctions)
library(devtools)
library(magrittr)
## devtools::load_all("/home/cew54/RP/coloc")
source("dirs.rb")
library(data.table)

args <- getArgs(default=list(block="chr19_block6",redo=FALSE))

files=file.path(DIR,"bytrait","CD_DeLange_28067908_1-hg38.tsv.gz",paste0(args$block,".rds"))

## already imputed
ofiles=sub(".rds","_imputed.rds",files)
of_exists=file.exists(ofiles)
files=files[!of_exists]
if(!length(files))
  stop("no files to impute, stopping ",args$block)

sizes=file.info(files)$size
files=files[sizes>0]
if(!length(files))
  stop("no files to impute, stopping ",args$block)

## load target snps
target_snps=fread(file.path(DIR,"reference","byblock",paste0(args$block,".csv")))
snps=with(target_snps,paste(V1,V2,sep=":"))

## get MAF, LD stats
r2.thr=0.8
ss=strsplit(snps,":") %>% do.call("rbind",.) %>%as.data.frame()
tmp=tempfile()
tmp.plink=tempfile()
write.table(ss, file=tmp, quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")
paste("./make_ref_subset.rb ",args$block,tmp,tmp.plink) %>% system()

## reset variant identifier in bim file
file.rename(paste0(tmp.plink,".bim"),paste0(tmp.plink,".bim.orig"))
paste0("cat ",tmp.plink,".bim.orig | awk 'BEGIN{ OFS=\"\t\" } {print $1,$1 \":\" $4,$3,$4,$5,$6}' > ",tmp.plink,".bim") %>%
  system()

## read
library(annotSnpStats)
ref=annot.read.plink(tmp.plink)

MAF=col.summary(ref)$MAF; names(MAF)=sub("chr","",colnames(ref))
drop=which(MAF==0)
if(length(drop)) {
  MAF=MAF[-drop]
  ref=ref[,-drop]
  target_snps=target_snps[-drop,]
  snps=names(MAF)
}
LD=cor(as(ref,"numeric")) ## look at adding lambda here
dimnames(LD) %<>% lapply(., function(x) sub("chr","",x))
lambda=2/sqrt(503)
LDl=(1-lambda) * LD + lambda * diag(nrow(LD))

################################################################################

for(ifile in files) {
  message("processing ",basename(dirname(ifile)))
  ofile=sub(".rds","_imputed.rds",ifile)
  if(file.exists(ofile) & !args$redo)
    next

  x=readRDS(ifile)
  drop=which(x$varbeta==0 | is.infinite(x$beta))
  if(length(drop)) {
    x$snp=x$snp[-drop]
    x$beta=x$beta[-drop]
    x$varbeta=x$varbeta[-drop]
  }
  str(x)

  ## find snps to impute
  snps_toimp=setdiff(colnames(LD),x$snp[x$ambig==FALSE])
  message("snps to impute: ",length(snps_toimp)," out of ",ncol(LD))
  if(length(snps_toimp)==0)
    saveRDS(x,file=ofile)

  ## find tags (thin a little if needed)
  tags=setdiff(colnames(LD),snps_toimp)
  LD2=LD[tags,tags]^2
  diag(LD2)=0

  drop=which(LD2>r2.thr, arr.ind=TRUE)
  tagdrop=rep(FALSE,length(tags))
  while(nrow(drop)) {
    idrop=sample(drop[,"row"],1)
    tagdrop[idrop]=TRUE
    drop=drop[ drop[,"row"]!=idrop & drop[,"col"]!=idrop , ]
  }
  tags=tags[!tagdrop]
  message("tags selected with max r2 <= ",r2.thr," : ",length(tags))

  cl=LDl[snps_toimp,tags]
  Z_tags=(x$beta/sqrt(x$varbeta))[ match(tags,x$snp) ]
  ## A=solve(LDl[tags,tags], Z_tags)
  ## Z_imp=cl %*% matrix(A,ncol=1)
  Cl=solve(LDl[tags,tags])
  Z_imp=cl %*% Cl %*% matrix(Z_tags,ncol=1)
  se_imp=1/sqrt(2*MAF[rownames(Z_imp)]*(1-MAF[rownames(Z_imp)])*x$N)
  beta_imp=Z_imp*se_imp

  ## prediction quality
  r2pred=diag(cl %*% Cl %*% t(cl))

  ## where alleles were ambiguous, use Z_imp to fix sign
  idx=match(snps_toimp,x$snp)
  x$beta[idx[!is.na(idx)]]=sign(beta_imp[!is.na(idx)]) * abs(x$beta[idx[!is.na(idx)]])
  Z_imp=Z_imp[is.na(idx),1]
  se_imp=se_imp[is.na(idx)]
  beta_imp=beta_imp[is.na(idx)]
  r2pred=r2pred[is.na(idx)]

  ## append to dataset
  x$r2pred=c(rep(NA,length(x$snp)),r2pred)
  x$snp=c(x$snp,names(Z_imp))
  x$varbeta=c(x$varbeta, se_imp^2)
  x$beta=c(x$beta, beta_imp)
  x$Z=x$beta/sqrt(x$varbeta)
  ## x$imputed=ifelse(x$snp %in% snps_toimp,ifelse(x$ambig)

  ## save
  saveRDS(x,file=ofile)
}
