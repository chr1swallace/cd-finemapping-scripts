## qR.rb -p skylake-himem -t "06:00:00" -c 4 -r -y 1-4 -n which -j export
library(randomFunctions)
library(devtools)
library(magrittr)
devtools::load_all("/home/cew54/RP/coloc")
source("dirs.rb")
library(pbapply)
## source("count_cores.R")
## library(parallel); options(mc.cores=ncores)
options(warn=2)

## find all susie files
d=file.path(DIR,"byblock_bytrait")
args=getArgs(defaults=list(which=1),numeric="which")

WANTED=list( "CD_DeLange_28067908_1-hg38.tsv.gz")

## add blocks
block_limits=fread("blocks.txt")
setnames(block_limits,c("block","chr","start","end"))

wanted=WANTED[[args$which]]

## for(wanted in WANTED[-c(1)]) {
## mclapply(WANTED,function(wanted) {

## debug
## files=file.path(d,"chr3_block32",paste0(wanted,".RData")); i=1

files=list.files(d,pattern=wanted,full=TRUE,recursive=TRUE)
blocks_fm=dirname(files) %>% basename()
names(files)=blocks_fm
blocks_single=setdiff(block_limits$block, blocks_fm)
  ## remove chr6_block23-25 (MHC)
  blocks_fm %<>% grep("chr6_block2[345]", ., value=TRUE, invert=TRUE)
  blocks_single %<>% grep("chr6_block2[345]", ., value=TRUE, invert=TRUE)
message("dataset ",wanted,"\tfinemapped files found: ",length(files))
message("dataset ",wanted,"\tsingle cv blocks: ",length(blocks_single))

## also load summary data for these snps
read_block=function(block) {
  message("block ",block)
  (load(f=file.path(DIR,"byblock",paste0("datasingle_",block,".RData")))) # DATA, mix of orig, and imputation
  ret=with(DATA[[wanted]],
           data.table(block=block,snp=snp,beta=beta,varbeta=varbeta))
  if("r2pred" %in% names(DATA[[wanted]]))
    ret$r2pred=DATA[[wanted]]$r2pred
 (load(files[block])) # result
  s=result # susie result if susie passed QC, or data.table of single otherwise
  if(class(s)=="susie") {
  s.null=is.null(s)
  ## add check on sld and set s=NULL if too big
  if(max(s$sld^2) > 0.5 || !length(s$sets$cs))
    s.null=TRUE
  message("\t",if(s.null) { "null" } else { "susie" })

  ## set s=NULL if p values for cred sets too big
  if(!s.null) {
    cs=s$sets$cs %>% lapply(., names) %>% unlist()
    m=match(cs,DATA[[wanted]]$snp)
    p=(DATA[[wanted]]$beta/sqrt(DATA[[wanted]]$varbeta))[m] %>% abs() %>% pnorm(., lower=FALSE)
    if(max(p) > 0.01)
      s.null=TRUE
    s_dt=data.table(snp=names(s$pip),pip=s$pip)
    sets=s$sets
    ## sets0=susie_get_cs(s)
    for(i in s$sets$cs_index)
      s_dt[[paste0("pip_set_",i)]]=s$alpha[ i, s_dt$snp]
    for(i in seq_along(s$sets$cs_index))
      s_dt[[paste0("set_",s$sets$cs_index[1])]]=s_dt$snp %in% sets$cs[[i]]
    sets=grep("set_",colnames(ret),value=TRUE)
    for(v in sets)
      s_dt[[v]][is.na(s_dt[[v]])] = FALSE
    s_dt[is.na(pip),pip:=FALSE]
    ret=merge(ret,s_dt,all.x=TRUE,by="snp")
  }
  }
  tmp=DATA[[wanted]][c("beta","varbeta","snp","N","s","type")]
  single=suppressWarnings(finemap.abf(tmp))
  s_dt=data.table(snp=single$snp,single.pp=single$SNP.PP)
  ret=merge(ret,s_dt,all.x=TRUE,by="snp")
  ## debug
  ## ret=copy(ret)
  ## ret[,z:=beta/sqrt(varbeta)]; ret[,p:=pnorm(-abs(z))]; ret[,bp:=as.numeric(sub("1:","",snp))]
  ## ggplot(ret, aes(x=bp,y=-log10(p),col=set_1)) + geom_point()
  ret
}

results=lapply(blocks_fm, read_block) %>% rbindlist(., fill=TRUE)

## single for everything
data=fread(cmd="zcat /home/cew54/share2/02-Processed/CD_DeLange_28067908_1-hg38.tsv.gz")
setkey(data, CHR38, BP38)
head(data)
data[,snp:=paste(CHR38,BP38,sep=":")][,beta:=hm_BETA][,varbeta:=SE^2]
SINGLE=vector("list",length(blocks_single))
names(SINGLE)=blocks_single
for(b in blocks_single) {
  i=which(block_limits$block==b)
  d=data[ CHR38==block_limits$chr[i] & BP38 >= block_limits$start[i] & BP38 <= block_limits$end[i] ]
  dlist=as.list(d[!is.na(beta) & !is.na(varbeta) & varbeta>0,.(snp,beta,varbeta)])
  dlist$N=12194 + 28072
  dlist$s=12194/(12194 + 28072)
  dlist$type="cc"
  ret=with(dlist,
           data.table(block=b,snp=snp,beta=beta,varbeta=varbeta))
  single=suppressWarnings(finemap.abf(dlist))
  s_dt=data.table(snp=single$snp,single.pp=single$SNP.PP)
  SINGLE[[b]]=merge(ret,s_dt,all.x=TRUE,by="snp")
}
FM=vector("list",length(blocks_fm))
names(FM)=blocks_fm
for(b in blocks_fm) {
  i=which(block_limits$block==b)
  d=data[ CHR38==block_limits$chr[i] & BP38 >= block_limits$start[i] & BP38 <= block_limits$end[i] ]
  dlist=as.list(d[!is.na(beta) & !is.na(varbeta) & varbeta>0,.(snp,beta,varbeta)])
  dlist$N=12194 + 28072
  dlist$s=12194/(12194 + 28072)
  dlist$type="cc"
  ret=with(dlist,
           data.table(block=b,snp=snp,beta=beta,varbeta=varbeta))
  fm=suppressWarnings(finemap.abf(dlist))
  s_dt=data.table(snp=fm$snp,single.pp=fm$SNP.PP)
  FM[[b]]=merge(ret,s_dt,all.x=TRUE,by="snp")
}


sresults=rbindlist(SINGLE)
fmresults=rbindlist(FM)

save(sresults, results, fmresults, file="mikhail-cd-temp.RData")
if(!interactive())
  q("no")

################################################################################

## start here

(load("mikhail-cd-temp.RData"))

## lresults=fread("ibd_for_mikhail.csv")
head(results)
head(sresults)
head(fmresults)

m=merge(results, fmresults[,.(block,snp,single.pp2=single.pp)], by=c("block","snp"), all.x=TRUE)
dim(results)
dim(fmresults)
dim(m)
head(m)
m[is.na(single.pp), single.pp:=single.pp2][,single.pp:=NULL]

results2=rbind(results,sresults,fill=TRUE)

results2=results2[rev(order(snp,r2pred))]
  results2=results2[!duplicated(snp)] # keep non-imputed, ie where r2pred=NA
message("dataset ",wanted,"\tfiltered data rows: ",nrow(results2))
tmp=results2[,grep("^set_",names(results2),value=TRUE,invert=TRUE), with=FALSE]
fwrite(tmp, file="cd_for_mikhail.csv")
system("gzip cd_for_mikhail.csv")

bb=tmp[,.(all(is.na(pip_set_6)),all(is.na(pip_set_8)),all(is.na(pip_set_9)),all(is.na(pip_set_10))),
       by="block"]
bb=bb[,-1] %>% as.matrix()
cor(bb)
