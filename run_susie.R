#!/usr/bin/env Rscript

library(randomFunctions)
library(devtools)
library(magrittr)
devtools::load_all("/home/cew54/RP/coloc")
source("~/E/dirs.rb")
source("~/E/common.R")

args <- getArgs(default=list(block="chr21_block20",redo=TRUE)) # highlighted by Victoria

idir <- file.path(DIR,"bytrait")
stopfile=file.path(DIR,"byblock",paste0(args$block,".stop"))
tmpfile_mafld=file.path(DIR,"byblock",paste0("tmp_mafld_",args$block,".RData"))
outfile_single=file.path(DIR,"byblock",paste0("datasingle_",args$block,".RData"))
outfile_susie=file.path(DIR,"byblock",paste0("datasusie_",args$block,".RData"))

## files may be imputed or not. Use imputed if available
files=file.path(idir,"CD_DeLange_28067908_1-hg38.tsv.gz",paste0(args$block,".rds"))
## replace by imputed files, if they exist
  imp_files=sub(".rds","_imputed.rds",files)
  imp_ex=file.exists(imp_files)
  files=ifelse(imp_ex,imp_files,files)

## NB originally written to process multiple traits in parallel. In this use
## case, a single trait, DATA will be a list of length 1. Safer than rewriting
## all the code.

DATA <- lapply(files, readRDS)

################################################################################

DATA[]=lapply(DATA[], function(d) {
  wh <- which(is.na(d$beta) | is.na(d$varbeta) | d$varbeta==0 |
              is.infinite(d$beta))
  if("MAF" %in% names(d))
    wh=c(wh, which(is.na(d$MAF) | d$MAF==0 | d$MAF==1)) %>% unique()
  if("r2pred" %in% names(d))
    wh=c(wh,which(d$r2pred < 0.5)) %>% unique()
  if(length(wh)) {
    nm <- intersect(names(d),c("beta","varbeta","snp","MAF","Z","imputed","r2pred"))
    d[nm] <- lapply(d[nm], function(x) x[-wh])
  }
  ## remove dup snps
  if(length(wh <- which(duplicated(d$snp)))) {
    nm <- intersect(names(d),c("beta","varbeta","snp","MAF","Z","imputed","r2pred"))
    d[nm] <- lapply(d[nm], function(x) x[-wh])
  }
  d$z=d$beta/sqrt(d$varbeta)
  ## add position
  if(!("position" %in% names(d)))
    d$position=as.numeric(sub(".*:","",d$snp))
  d
})

## only run if min p < 5e-8 for one or more disease traits
## thin datasets to only min p < 1e-6
keep=sapply(DATA, function(d) any(abs(d$z)>4.89)) # minimum p=1e-6 or less
trait.keep=sapply(DATA[traits], function(d) any(abs(d$z)>qnorm(5e-8/2,lower=FALSE))) %>%
                                        any() # minimum p=5e-8 or less
if(!all(keep))
  DATA=DATA[keep]
if(!trait.keep || !length(intersect(names(DATA),traits)))  {
  ## paste("touch", outfile) %>% system()
  ## TODO REVERT WHEN ALIGN FINISHED
  paste("touch", stopfile) %>% system()
  stop("no traits with significance, stopping ",args$block)
}

lapply(DATA, function(d) summary(abs(d$z)))
lapply(DATA, function(d) summary(d$MAF))
message("after limiting to min p < 1e-6, datasets remaining: ",length(DATA))

## save data for basic analysis using all available trait data
save(DATA, file=outfile_single)
################################################################################

tnm=intersect(traits,names(DATA)) # disease traits
enm=setdiff(names(DATA), traits) # eqtls, biomarkers

snps=DATA[[1]]$snp

## get MAF, LD stats
message("reading MAF, LD")
runorload(tmpfile_mafld, c("MAF","LD"), depfile=NULL, {
  ss=strsplit(snps,":") %>% do.call("rbind",.) %>%as.data.frame()
  ss$V1 %<>% paste0("chr",.)
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

  message("calculating MAF, LD (can be slow)")
  MAF=col.summary(ref)$MAF; names(MAF)=sub("chr","",colnames(ref))
  drop=which(MAF==0)
  if(length(drop)) {
    MAF=MAF[-drop]
    ref=ref[,-drop]
    snps=names(MAF)
  }
  LD=cor(as(ref,"numeric"))
  ## LD=ld(ref,ref,stat="R",symmetric=TRUE)
  dimnames(LD) %<>% lapply(., function(x) sub("chr","",x))
})


################################################################################

## function to limit datasets to snp subsets
snps=colnames(LD)
message("limiting ",length(DATA)," datasets to common snps: ",length(snps))
limsnps <- function(d,snps) {
    wh <- match(snps, d$snp) %>% setdiff(.,NA)
    nm <- intersect(names(d),c("beta","varbeta","snp","MAF","z","r2pred","imputed","Z","position"))
    d[nm] <- lapply(d[nm], function(x) x[wh])
    d
}
DATA <- lapply(DATA, limsnps, snps=snps)

## only run if minp < 1e-6
mp=sapply(DATA,function(d) pnorm(-max(abs(d$z)))*2)
if(any(mp > 1e-6))
  DATA=DATA[mp < 1e-6]

## add/replace MAF where missing/incomplete
fixmaf=function(d) {
  if("MAF" %in% names(d) && length(d$MAF)==length(d$snp) && !any(is.na(d$MAF)))
    return(d)
  d$MAF=MAF[d$snp]
  d
}
DATA %<>% lapply(., fixmaf)

SUSIE.FM=PROPOS=vector("list",length(DATA))
names(SUSIE.FM)=names(PROPOS)=names(DATA)
odir=file.path(DIR,"byblock_bytrait",args$block)
tmpfiles=file.path(odir,paste0(names(DATA),".RData"))
tmp.exists=file.exists(tmpfiles)
table(tmp.exists)

## run susie once per dataset
d=c(DATA[[1]], list(LD=LD[d$snp,d$snp], MAF=MAF[d$snp]))

if(!file.exists(odir))
  dir.create(odir)

i=1
  if(file.exists(tmpfiles[i])) # perhaps created by another script in parallel
    return(NULL)
  message("trait ",i," / ",length(tmpfiles),"\tsnps: ",length(d$snp))
  message(names(DATA)[i])
  message("\t",date())
  int=intersect(c("beta","varbeta","snp","N","s","type","position"),names(DATA[[i]]))
  d=DATA[[i]][int]
  d <- c(d, list(LD=LD[d$snp,d$snp], MAF=MAF[d$snp]))
  ## prop_pos=check_alignment(d)
  prop_pos=check_alignment(d,do_plot=FALSE)
  ## conditional - can we skip susie? is there dodgy data?
  cond=finemap.signals(d, method="cond")
  if(length(cond) <= 1) { # can skip
    result=finemap.signals(d, method="single", return.pp=TRUE)
    meth="cond, single"
  }
  if(length(cond) > 1) {
    s=try(runsusie(d,maxit=MAXIT,trimpeaks=3,
                   trimz=1,repeat_until_convergence=FALSE,check_prior=FALSE),silent=TRUE)
    if("try-error" %in% class(s)) { # do it the old way
      message("trait ",i,": try-error ",args$block," ",names(DATA)[i])
      print(s)
      if(any(diff(abs(cond)) > 2)) { # do we trust cond result?
        message("trait ",i,": possible issue with data found - results get *stronger* after conditioning")
        result=finemap.signals(d, method="single", return.pp=TRUE)
        meth="single, sfail+condfail"
      } else {
        result=finemap.signals(d, method="cond", return.pp=TRUE)
        meth="cond, sfail"
      }
    } else {
      result=s
      meth="susie"
    }
  }
  if(!file.exists(tmpfiles[i])) {
    save(result,prop_pos,file=tmpfiles[i])
    message("trait ",i,": saved\t",tmpfiles[i],"\n\t",date())
  }
