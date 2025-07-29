#!/usr/bin/env Rscript
## LONGVAR cfDNAme age prediction

## Load necessary libraries
packages <- c("data.table", "readxl", "GenomicRanges", "methylKit", "genomation", "ggsci", "impute", "ggplot2", "stringr", "dplyr", "tidyr", "ggpubr", "rstatix", "ggpmisc", "caret", "e1071", "glmnet")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
metafh <- args[1]
covfiledir <- args[2]
outdir <- args[3]
captureAllCpGs <- args[4]

heytoday <- Sys.Date()

ref_m <- fread(captureAllCpGs,header=FALSE,sep="\t",data.table=FALSE)
colnames(ref_m) <- c("chr","start","end","flag")
ref_m <- ref_m[!(ref_m$chr %in% c("chrX")),]
rownames(ref_m) <- paste0(ref_m$chr,"_",ref_m$start,"_",ref_m$end)
ref_m_gr <- makeGRangesFromDataFrame(ref_m,keep.extra.columns=FALSE,ignore.strand=TRUE)

metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)
T1metainfo <- subset(metainfo, Year=="T1")
T2metainfo <- subset(metainfo, Year=="T2")


ReadCov2Meobj <- function(metainfo, covfiledir, ref_m_gr, fraction_of_na_regions=0.05, fraction_of_na_samples=0.1, keepallloci="yes", keepallsample="yes") {
    gccodes <- metainfo$GCcode
    covfhs <- paste0(covfiledir, "/", gccodes, "_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
    covfhsl <- as.list(covfhs)

    sampleidsl <- as.list(metainfo$Samplename)
    samtreatment <- ifelse(metainfo$Sex=="Female", 0, 1)
    batchobj <- methRead(location=covfhsl, sample.id=sampleidsl, pipeline="bismarkCoverage", assembly='hg38', context="CpG", treatment=samtreatment, mincov=10)

    regionCounts_all <- regionCounts(batchobj, ref_m_gr)
    filtered_cpgCounts_all <- filterByCoverage(regionCounts_all, lo.count = 0, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    meth_all_cpgs <- methylKit::unite(filtered_cpgCounts_all, destrand = FALSE, min.per.group = 1L)

    m <- methylKit::getData(meth_all_cpgs)
    all_cpgs <- paste(m$chr, m$start, m$end, sep = "_")
    meth_all_cpgs.perc <- percMethylation(meth_all_cpgs)
    rownames(meth_all_cpgs.perc) <- all_cpgs
    meth_all_cpgs.perc = unique(meth_all_cpgs.perc)

    min_cov_locus = 50
    min_cov_sample = 10
    fraction_of_na_regions = fraction_of_na_regions
    fraction_of_na_samples = fraction_of_na_samples
    covs = m[,grepl('coverage', colnames(m))]
    covs = apply(covs, 2, as.numeric)
    mean_covs_per_sample = apply(covs, 2, function(x){mean(x[!is.na(x)])})
    mean_covs_per_region = apply(covs, 1, function(x){mean(x[!is.na(x)])})
    frequency_of_na_per_region = apply(covs, 1, function(x){length(x[is.na(x)])})
    frequency_of_na_per_sample = apply(covs, 2, function(x){length(x[is.na(x)])})
    names(mean_covs_per_region) = rownames(meth_all_cpgs.perc)
    names(frequency_of_na_per_region) = rownames(meth_all_cpgs.perc)
    names(frequency_of_na_per_sample) = colnames(meth_all_cpgs.perc)
    names(mean_covs_per_sample) = colnames(meth_all_cpgs.perc)
    samples = colnames(meth_all_cpgs.perc)
    kept_samples = samples[(frequency_of_na_per_sample <= nrow(meth_all_cpgs.perc)*fraction_of_na_samples) & (mean_covs_per_sample >= min_cov_sample)]
    Covfailed_samples = samples[!samples%in%kept_samples]
    metainfo$Covfailed <- ifelse(metainfo$Samplename %in% kept_samples, "FALSE","TRUE")

    loci = rownames(meth_all_cpgs.perc)

    kept_loci = loci[(frequency_of_na_per_region <= nrow(meth_all_cpgs.perc)*fraction_of_na_regions) & (mean_covs_per_region >= min_cov_locus)]
    Covfailed_loci = loci[!loci%in%kept_loci]

    print(paste0(length(Covfailed_loci), ' loci failed, ', length(kept_loci), ' loci passed.'))
    print(Covfailed_loci)
    print(paste0(length(Covfailed_samples), ' samples failed, ', length(kept_samples), ' samples passed.'))
    print(Covfailed_samples)

    #meth_regions_perc_ft = meth_all_cpgs.perc[na.omit(kept_loci),]
    if (keepallloci == "yes") {
        meth_regions_perc_ft = meth_all_cpgs.perc
    } else {
        meth_regions_perc_ft = meth_all_cpgs.perc[na.omit(kept_loci),]
    }
    
    if (keepallsample == "yes") {
        meth_regions_perc_ft2 = meth_regions_perc_ft
    } else {
        meth_regions_perc_ft2 = meth_regions_perc_ft[,na.omit(kept_samples)]
    }

    imputed = impute.knn(data = meth_regions_perc_ft2, k = 10)
    meth_all_cpgs_imputed  = imputed$data
    res = meth_all_cpgs_imputed
    return(res)
}

#####Removing CpGs are poorly covered across all samples
Allsam <- ReadCov2Meobj(metainfo, covfiledir, ref_m_gr, fraction_of_na_regions=0.001, fraction_of_na_samples=0.1, keepallloci="no", keepallsample="yes")
ref_m2 <- ref_m[rownames(ref_m) %in% rownames(Allsam),]
ref_m_gr2 <- makeGRangesFromDataFrame(ref_m2,keep.extra.columns=FALSE,ignore.strand=TRUE)

#####Selecting CpGs without distinguishing Female and Male
T1res <- ReadCov2Meobj(T1metainfo, covfiledir, ref_m_gr2, fraction_of_na_regions=0.01, fraction_of_na_samples=0.1, keepallloci="no", keepallsample="yes")


#######Spearman correlation between Age and each CpGs in T1 set
GetTrainSetRho <- function(trainmetadt, trainobj, trainset_name) {
    Tage <- trainmetadt[trainmetadt$Samplename==colnames(trainobj),]$Age
    Tsex <- trainmetadt[trainmetadt$Samplename==colnames(trainobj),]$Sex
    Tsamname <- trainmetadt[trainmetadt$Samplename==colnames(trainobj),]$Samplename

    trainobj <- trainobj[ ! apply(trainobj, 1, sd, na.rm=TRUE) < 0.01,]

    Tagescor <- cor(t(trainobj), Tage, method="spearman")
    Tagescor <- data.frame(Tagescor)
    Tagescor$CpG <- rownames(Tagescor)
    Tagescor <- Tagescor[order(Tagescor$Tagescor),]
    Tagescor$CpG <- factor(Tagescor$CpG,levels=c(Tagescor$CpG))
    
    heytoday <- Sys.Date()
    outfig1 <- paste0(outdir,"/LONGVAR.", trainset_name,"_", nrow(trainmetadt), "sam.Age_spearmancor_",nrow(Tagescor),"CpGs.", heytoday, ".pdf")
    pdf(outfig1,height=4.5,width=8)
    p1 <- ggplot(Tagescor,aes(x=CpG,y=Tagescor))+geom_point(fill="black",shape=20,aes(size=abs(Tagescor)>0.5))+geom_hline(yintercept=0,color="gray")+geom_hline(yintercept=0.5,color="red3",linetype="dashed")+geom_hline(yintercept=-0.5,color="red3",linetype="dashed")+geom_hline(yintercept=0.3,color="red2",linetype="dashed")+geom_hline(yintercept=-0.3,color="red2",linetype="dashed")+geom_hline(yintercept=0.1,color="red1",linetype="dashed")+geom_hline(yintercept=-0.1,color="red1",linetype="dashed")+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=12,color="black"),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.text=element_text(size=11),axis.title=element_text(size=14,color="black"),legend.position="none",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="CpG site",y="Spearman correlation")+scale_y_continuous(limits=c(-0.71,0.71),breaks=c(-0.7,-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7))+expand_limits(x= c(min(as.numeric(Tagescor$CpG))-200, max(as.numeric(Tagescor$CpG))+200))+scale_size_manual(values=c(0.5,0.9))
    print(p1)
    dev.off()

    smaxcpg <- Tagescor[Tagescor$Tagescor==max(Tagescor$Tagescor),]$CpG
    smincpg <- Tagescor[Tagescor$Tagescor==min(Tagescor$Tagescor),]$CpG

    res_smaxcpg <- trainobj[rownames(trainobj)==smaxcpg,]
    res_smaxcpgdt <- data.frame(res_smaxcpg)
    colnames(res_smaxcpgdt) <- "Meper"
    res_smaxcpgdt$Samplename <- rownames(res_smaxcpgdt)
    maxdt <- merge(res_smaxcpgdt,trainmetadt,by=c("Samplename"))
    maxdt$Sex <- factor(maxdt$Sex,levels=c("Female","Male"))

    outfig2 <- paste0(outdir, "/LONGVAR.", trainset_name,"_", nrow(trainmetadt), "sam.Age_spearmancor_CpGs.highestCor.", heytoday, ".pdf")
    xupper <- max(maxdt$Meper)+1
    xlower <- min(maxdt$Meper,0)
    pdf(outfig2,height=4.5,width=5.4)
    p2 <- ggplot(maxdt,aes(x=Meper,y=Age))+geom_point(aes(fill=Sex),shape=21,size=1.5)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Methylation (%)",y="Age")+scale_fill_manual(values=c("darkblue","darkred"),drop=FALSE)+ggtitle(smaxcpg)+xlim(xlower,xupper)+geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="spearman")
    print(p2)
    dev.off()

    res_smincpg <- trainobj[rownames(trainobj)==smincpg,]
    res_smincpgdt <- data.frame(res_smincpg)
    colnames(res_smincpgdt) <- "Meper"
    res_smincpgdt$Samplename <- rownames(res_smincpgdt)
    mindt <- merge(res_smincpgdt,trainmetadt,by=c("Samplename"))
    mindt$Sex <- factor(mindt$Sex,levels=c("Female","Male"))

    outfig2 <- paste0(outdir, "/LONGVAR.", trainset_name,"_", nrow(trainmetadt), "sam.Age_spearmancor_CpGs.lowestCor.", heytoday, ".pdf")
    xupper <- max(mindt$Meper)+1
    xlower <- min(mindt$Meper)-1
    pdf(outfig2,height=4.5,width=5.4)
    p2 <- ggplot(mindt,aes(x=Meper,y=Age))+geom_point(aes(fill=Sex),shape=21,size=1.5)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Methylation (%)",y="Age")+scale_fill_manual(values=c("darkblue","darkred"),drop=FALSE)+ggtitle(smincpg)+xlim(xlower,xupper)+geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+stat_correlation(label.x="right",use_label(c("R", "P")),method="spearman")
    print(p2)
    dev.off()

    return(Tagescor)
}

T1agescor <- GetTrainSetRho(trainmetadt=T1metainfo, trainobj=T1res, trainset_name="T1")


######select with correlation and do ELN
T1age <- T1metainfo[T1metainfo$Samplename==colnames(T1res),]$Age
T1sex <- T1metainfo[T1metainfo$Samplename==colnames(T1res),]$Sex
T1samname <- T1metainfo[T1metainfo$Samplename==colnames(T1res),]$Samplename

strongcordt <- T1agescor[T1agescor$Tagescor>=0.1 | T1agescor$Tagescor<=-0.1,]
strongCpGs <- strongcordt$CpG

res_cor <- T1res[rownames(T1res) %in% strongCpGs,]
agemat <- as.matrix(t(res_cor))
Sexfactor <- matrix(ifelse(T1sex=="Female",0,1))
colnames(Sexfactor) <- c("Sex")
agemat1 <- cbind(Sexfactor,agemat)

set.seed(88886)
fitpenalty <- c(0, rep(1, ncol(agemat)))

cvfit <- cv.glmnet(agemat1, T1age, penalty.factor=fitpenalty, alpha=0.85, type.measure="mse", nfolds = 10)
prelambda <- round(cvfit$lambda.min, 4)
print(prelambda)

FeatureGLMLoo <- function(ftmat, fitpenalty, outcome, alpha, lambda) {
    alpha <- as.numeric(alpha)
    n_train <- nrow(ftmat)
    allpredval <- c()

    loopreds <- sapply(1:n_train, function(x) {
        traindt <- ftmat[-x,]
        trainage <- outcome[-x]
        #wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
	    testage <- outcome[x]
        testdt <- t(as.data.frame(ftmat[x,]))
	
        rglm.model <- glmnet(traindt, trainage, penalty.factor=fitpenalty, alpha=alpha,lambda=lambda)
        rglm.preds.res <- predict(rglm.model, testdt, type="response", s=lambda)
        #rglm.preds.coef <- predict(rglm.model, testdt, type="coefficient", s=lambda)

        allpredval <- c(allpredval, rglm.preds.res)
        return(allpredval)
    })

    outdt <- as.data.frame(loopreds)
    colnames(outdt) <- c("Methylation_age")
    outdt$Chronological_age <- outcome

    return(outdt)
}

T1pred <- FeatureGLMLoo(ftmat=agemat1,fitpenalty=fitpenalty, outcome=T1age, alpha=0.85, lambda=prelambda)
T1pred$Sex <- T1sex
T1pred$Samplename <- T1samname
cor(T1pred$Methylation_age, T1pred$Chronological_age)

T1pred$Delta <- T1pred$Methylation_age-T1pred$Chronological_age
T1pred_MAE <- round(sum(abs(T1pred$Methylation_age-T1pred$Chronological_age))/nrow(T1pred),4)
T1pred_MAE

T1pred_MAE_F <- round(sum(abs(T1pred[T1pred$Sex=="Female",]$Methylation_age-T1pred[T1pred$Sex=="Female",]$Chronological_age))/nrow(T1pred[T1pred$Sex=="Female",]),4)
T1pred_MAE_M <- round(sum(abs(T1pred[T1pred$Sex=="Male",]$Methylation_age-T1pred[T1pred$Sex=="Male",]$Chronological_age))/nrow(T1pred[T1pred$Sex=="Male",]),4)
print(paste0("MAE Female: ", T1pred_MAE_F))
print(paste0("MAE Male: ", T1pred_MAE_M))

T1pred_MAE_label <- deparse(bquote(MAE : ~ .(T1pred_MAE)))

outfig3 <- paste0(outdir,"/LONGVAR.T1_",nrow(T1pred),"sam.ELN.LOO_prediction.",heytoday,".pdf")
pdf(outfig3,height=4.6,width=5.8)
p3 <- ggplot(T1pred,aes(x=Chronological_age,y=Methylation_age))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_abline(slope=1, intercept=0,color="gray",linetype="dotted")+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Chronological age",y="Predicted methylation age")+scale_fill_manual(values=c("darkblue","darkred"))+geom_smooth(method=lm, se=FALSE,color="black", linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+scale_y_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+annotate("text", x=25, y=92, parse=TRUE, label=T1pred_MAE_label,size=4)
print(p3)
dev.off()

outfig31 <- paste0(outdir,"/LONGVAR.T1_",nrow(T1pred),"sam.ELN.LOO_prediction.DeltaAge.",heytoday,".pdf")
pdf(outfig31,height=4.6,width=6)
p31 <- ggplot(T1pred,aes(x=Chronological_age,y=Delta))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_hline(yintercept=0,color="gray",linetype="dotted")+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Chronological age",y=expression(Delta ~ age))+scale_fill_manual(values=c("darkblue","darkred"))+scale_color_manual(values=c("darkblue","darkred"))+geom_smooth(aes(color=Sex),method=lm, se=FALSE,linetype="dashed")+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))
print(p31)
dev.off()


T1finalrgml <- glmnet(agemat1, T1age, penalty.factor=fitpenalty, alpha=0.85, lambda=prelambda)
savedmodel <- paste0(outdir, "/LONGVAR.T1.",nrow(T1pred),"sam.ELN.Unpenalize_Sex.alpha0.85.rds")
saveRDS(T1finalrgml, savedmodel)

ref_m3 <- ref_m2[rownames(ref_m2) %in% strongCpGs,]
ref_m_gr3 <- makeGRangesFromDataFrame(ref_m3,keep.extra.columns=FALSE,ignore.strand=TRUE)

T2res <- ReadCov2Meobj(T2metainfo, covfiledir, ref_m_gr3, fraction_of_na_regions=0.01, fraction_of_na_samples=0.1, keepallloci="yes", keepallsample="yes")
#T2res <- ReadCov2Meobj(T2metainfo, covfiledir, ref_m_gr2, fraction_of_na_regions=0.01, fraction_of_na_samples=0.1, keepallloci="yes", keepallsample="yes")
T2age <- T2metainfo[T2metainfo$Samplename==colnames(T2res),]$Age
T2sex <- T2metainfo[T2metainfo$Samplename==colnames(T2res),]$Sex
T2samname <- T2metainfo[T2metainfo$Samplename==colnames(T2res),]$Samplename

T2rawmat <- T2res[rownames(T2res) %in% strongCpGs,]
T2agemat <- as.matrix(t(T2res))
Sexfactor <- matrix(ifelse(T2sex=="Female",0,1))
colnames(Sexfactor) <- c("Sex")
T2agemat2 <- cbind(Sexfactor,T2agemat)

T2rglm.pred.res <- predict(T1finalrgml, T2agemat2, type="response", s=prelambda)
T2pred <- as.data.frame(T2rglm.pred.res)
colnames(T2pred) <- c("Methylation_age")
T2pred$Chronological_age <- T2age
T2pred$Sex <- T2sex
T2pred$Samplename <- T2samname

T2rglm.pred.coef <- predict(T1finalrgml, T2agemat, type="coefficient", s=prelambda)
coefdt <- data.frame(T2rglm.pred.coef@Dimnames[[1]][T2rglm.pred.coef@i + 1], T2rglm.pred.coef@x)

T2pred$Delta <- T2pred$Methylation_age-T2pred$Chronological_age
T2pred_MAE <- round(sum(abs(T2pred$Methylation_age-T2pred$Chronological_age))/nrow(T2pred),4)

T2pred_MAE_F <- round(sum(abs(T2pred[T2pred$Sex=="Female",]$Methylation_age-T2pred[T2pred$Sex=="Female",]$Chronological_age))/nrow(T2pred[T2pred$Sex=="Female",]),4)
T2pred_MAE_M <- round(sum(abs(T2pred[T2pred$Sex=="Male",]$Methylation_age-T2pred[T2pred$Sex=="Male",]$Chronological_age))/nrow(T2pred[T2pred$Sex=="Male",]),4)
print(paste0("MAE Female: ", T2pred_MAE_F))
print(paste0("MAE Male: ", T2pred_MAE_M))

T2pred_MAE_label <- deparse(bquote(MAE : ~ .(T2pred_MAE)))

outfig4 <- paste0(outdir,"/LONGVAR.T2_", nrow(T2pred), "sam.ELN.T1model_prediction.", heytoday ,".pdf")
pdf(outfig4,height=4.6,width=5.8)
p4 <- ggplot(T2pred,aes(x=Chronological_age,y=Methylation_age))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_abline(slope=1, intercept=0,color="gray",linetype="dotted")+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Chronological age",y="Predicted methylation age")+scale_fill_manual(values=c("darkblue","darkred"))+geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+scale_y_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+annotate("text", x=25, y=92, parse=TRUE, label=T2pred_MAE_label,size=4)
print(p4)
dev.off()

outfig41 <- paste0(outdir,"/LONGVAR.T2_",nrow(T2pred),"sam.ELN.T1model_prediction.DeltaAge.",heytoday,".pdf")
pdf(outfig41,height=4.6,width=6)
p41 <- ggplot(T2pred,aes(x=Chronological_age,y=Delta))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_hline(yintercept=0,color="gray",linetype="dotted")+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Chronological age",y=expression(Delta ~ age))+scale_fill_manual(values=c("darkblue","darkred"))+scale_color_manual(values=c("darkblue","darkred"))+geom_smooth(aes(color=Sex),method=lm, se=FALSE,linetype="dashed")+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))
print(p41)
dev.off()

T1pred$Subject <- T1metainfo[T1metainfo$Samplename==colnames(T1res),]$Subject
T2pred$Subject <- T2metainfo[T2metainfo$Samplename==colnames(T2res),]$Subject
model1match <- merge(T1pred,T2pred,by=c("Subject","Sex"))
model1match$Diff <- model1match$Methylation_age.y-model1match$Methylation_age.x


outfig5 <- paste0(outdir,"/LONGVAR.T1_T2.matchsample.ELN.T1model_prediction.",heytoday,".pdf")
pdf(outfig5,height=4.6,width=5.8)
p5 <- ggplot(model1match,aes(x=`Methylation_age.x`, y=`Methylation_age.y`))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_abline(slope=1, intercept=1,color="gray",linetype="dotted")+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),plot.title=element_text(hjust=0.5),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=14),axis.title.y=element_text(size=14,color="black"),axis.title.x=element_text(size=14,color="black"),legend.position="right",legend.title=element_blank(),legend.background=element_blank(),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+scale_fill_manual(values=c("darkblue","darkred"))+scale_color_manual(values=c("darkblue","darkred"))+geom_smooth(aes(color=Sex),method=lm, se=FALSE,linetype="dashed")+labs(x="T1 predicted methylation age",y="T2 predicted methylation age")+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+scale_y_continuous(limits=c(11,101),breaks=seq(21,91,by=10))
print(p5)
dev.off()

outfig6 <- paste0(outdir,"/LONGVAR.T1_T2_delta.matchsample.ELN.T1model_prediction.",heytoday,".pdf")
pdf(outfig6,height=4.6,width=5.8)
p6 <- ggplot(model1match,aes(x=`Delta.x`, y=`Delta.y`))+geom_point(aes(fill=Sex),shape=21,size=1.5)+geom_abline(slope=1, intercept=1,color="gray",linetype="dotted")+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),plot.title=element_text(hjust=0.5),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=14),axis.title.y=element_text(size=14,color="black"),axis.title.x=element_text(size=14,color="black"),legend.position="right",legend.title=element_blank(),legend.background=element_blank(),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+scale_fill_manual(values=c("darkblue","darkred"))+scale_color_manual(values=c("darkblue","darkred"))+labs(x=expression(T1 ~ Delta ~ age),y=expression(T2 ~ Delta ~ age))+scale_x_continuous(limits=c(-22,16),breaks=seq(-20,15,by=5))+scale_y_continuous(limits=c(-22,16),breaks=seq(-20,15,by=5))
print(p6)
dev.off()

model1pred <- rbind(T1pred,T2pred)
model1pred <- model1pred[,c("Methylation_age","Samplename","Delta")]
upmetainfo <- merge(metainfo,model1pred,by=c("Samplename"))

agemefh <- paste0(outdir,"/LONGVAR.metainfo.merged.size.estmole.ageclock.",heytoday,".tsv")
write.table(upmetainfo,agemefh,quote=FALSE,sep="\t",row.names=FALSE)

commontheme <- function() {
    theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
      , panel.grid.major=element_blank(),panel.grid.minor=element_blank()
      , axis.text.y=element_text(size=14,color="black")
      , axis.text.x=element_text(size=14,color="black")
      , legend.text=element_text(size=12)
      , axis.title.x=element_text(size=16,color="black",vjust=-0.5)
      , axis.title.y=element_text(size=16,color="black",vjust=1.5)
      , strip.text=element_text(size=14,face="bold")
      , strip.background=element_rect(fill="white")
      , plot.title=element_text(hjust=0.5))
}

plot_ageme_prediction <- function() {
    ageme <- fread(agemefh,header=TRUE,sep="\t",data.table=FALSE)
    ageme$Year <- factor(ageme$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    ageme$Sex <- factor(ageme$Sex, levels=c("Female","Male"))
    t1dt <- subset(ageme,Year=="Year1")
    t2dt <- subset(ageme,Year=="Year2")

    T1pred_MAE <- round(sum(abs(t1dt$Methylation_age-t1dt$Age))/nrow(t1dt),4)
    T1pred_MAE_label <- deparse(bquote(MAE : ~ .(T1pred_MAE)))

    outfig1 <- paste0(outdir,"/LONGVAR.Year1_",nrow(t1dt),"sam.ELN.LOO_prediction.",heytoday,".pdf")
    pdf(outfig1,height=5,width=5)
    p1 <- ggplot(t1dt,aes(x=Age,y=Methylation_age))+geom_point(aes(color=Sex),size=2,shape=16)+geom_abline(slope=1, intercept=0,color="gray",linetype="dotted")+commontheme()+theme(legend.position=c(0.8,0.2))+labs(x="Chronological age",y="Predicted methylation age")+scale_color_manual(values=c("#EA7317","#2364AA"))+geom_smooth(method=lm, se=FALSE,color="black", linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=5)+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+scale_y_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+annotate("text", x=25, y=92, parse=TRUE, label=T1pred_MAE_label,size=5)
    print(p1)
    dev.off()

    T2pred_MAE <- round(sum(abs(t2dt$Methylation_age-t2dt$Age))/nrow(t2dt),4)
    T2pred_MAE_label <- deparse(bquote(MAE : ~ .(T2pred_MAE)))

    outfig2 <- paste0(outdir,"/LONGVAR.Year2_", nrow(t2dt), "sam.ELN.T1model_prediction.", heytoday ,".pdf")
    pdf(outfig2,height=5,width=5)
    p2 <- ggplot(t2dt,aes(x=Age,y=Methylation_age))+geom_point(aes(color=Sex),size=2,shape=16)+geom_abline(slope=1, intercept=0,color="gray",linetype="dotted")+commontheme()+theme(legend.position=c(0.8,0.2))+labs(x="Chronological age",y="Predicted methylation age")+scale_color_manual(values=c("#EA7317","#2364AA"))+geom_smooth(method=lm, se=FALSE,color="black",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=5)+scale_x_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+scale_y_continuous(limits=c(10,100),breaks=seq(20,90,by=10))+annotate("text", x=25, y=92, parse=TRUE, label=T2pred_MAE_label,size=5)
    print(p2)
    dev.off()
}
