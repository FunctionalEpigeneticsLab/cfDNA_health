#LONGVAR
# Deconvolution of CpG methylation data using EpiDISH and NNLS

## Load necessary libraries
packages <- c("data.table", "tidyr", "ggplot2", "cowplot", "ggpmisc","readxl", "GenomicRanges", "methylKit", "genomation", "ggsci", "impute", "EpiDISH", "nnls", "EMeth")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
    install.packages(new_packages)
}

## Load all required packages
invisible(lapply(packages, library, character.only = TRUE))

# Store today's date for output file naming
heytoday <- Sys.Date()

## Parse command line arguments (or set defaults for interactive use)
args <- commandArgs(trailingOnly = TRUE)
metafh <- args[1]
covfiledir <- args[2]
outdir <- args[3]
refavg <- args[4]


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

get_sample_mat <- function(metafh, covfiledir, refavg) {
    # Read sample metadata
    metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)
    gccodes <- metainfo$GCcode
    covfhs <- paste0(covfiledir, "/", gccodes, "_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
    covfhsl <- as.list(covfhs)
    sampleidsl <- as.list(metainfo$Samplename)
    samtreatment <- ifelse(metainfo$Year=="T1", 0, 1)
    # Read methylation data
    batchobj <- methRead(location=covfhsl, sample.id=sampleidsl, pipeline="bismarkCoverage", assembly='hg38', context="CpG", treatment=samtreatment, mincov=10)

    # Read reference methylation matrix
    ref_m <- fread(refavg,header=TRUE,sep="\t",data.table=FALSE)
    rownames(ref_m) <- paste0(ref_m$chr,"_",ref_m$start,"_",ref_m$end)
    ref_m_gr <- makeGRangesFromDataFrame(ref_m,keep.extra.columns=FALSE,ignore.strand=TRUE)
    regionCounts_all <- regionCounts(batchobj, ref_m_gr)
    filtered_cpgCounts_all <- filterByCoverage(regionCounts_all, lo.count = 0, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    meth_all_cpgs <- methylKit::unite(filtered_cpgCounts_all, destrand = FALSE, min.per.group = 1L)

    m <- methylKit::getData(meth_all_cpgs)
    all_cpgs <- paste(m$chr, m$start, m$end, sep = "_")
    meth_all_cpgs.perc <- percMethylation(meth_all_cpgs)
    rownames(meth_all_cpgs.perc) <- all_cpgs
    meth_all_cpgs.perc = unique(meth_all_cpgs.perc)

    # Check loci and samples by coverage and missingness
    min_cov_locus = 30
    min_cov_sample = 10
    fraction_of_na_regions = 0.05
    fraction_of_na_sample = 0.25
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
    kept_samples = samples[(frequency_of_na_per_sample < nrow(meth_all_cpgs.perc)*fraction_of_na_sample) & (mean_covs_per_sample > min_cov_sample)]
    Covfailed_samples = samples[!samples%in%kept_samples]
    metainfo$Covfailed <- ifelse(metainfo$Samplename %in% kept_samples, "FALSE","TRUE")

    loci = rownames(meth_all_cpgs.perc)

    kept_loci = loci[(frequency_of_na_per_region < nrow(meth_all_cpgs.perc)*fraction_of_na_regions) & (mean_covs_per_region > min_cov_locus)]
    Covfailed_loci = loci[!loci%in%kept_loci]

    print(paste0(length(Covfailed_loci), ' loci failed, ', length(kept_loci), ' loci passed.'))
    print(Covfailed_loci)
    print(paste0(length(Covfailed_samples), ' samples failed, ', length(kept_samples), ' samples passed.'))
    print(Covfailed_samples)

    # Impute missing values
    meth_regions_perc_ft = meth_all_cpgs.perc
    imputed = impute.knn(data = meth_regions_perc_ft, k = 10)
    meth_all_cpgs_imputed  = imputed$data

    return(meth_all_cpgs_imputed)
}

get_epidish_deconv_res <- function(metafh, covfiledir, refavg) {
    # Read sample metadata
    metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)

    # Get sample matrix and reference methylation matrix
    meth_all_cpgs_imputed <- get_sample_mat(metafh, covfiledir, refavg)

    ref_m <- fread(refavg,header=TRUE,sep="\t",data.table=FALSE)
    rownames(ref_m) <- paste0(ref_m$chr,"_",ref_m$start,"_",ref_m$end)
    ref_m_filter = ref_m[rownames(meth_all_cpgs_imputed),-which(names(ref_m) %in% c("chr","start","end"))]

    res = meth_all_cpgs_imputed

    # Deconvolution using EpiDISH
    predicted_cell_proportions_epidish <- epidish(as.matrix(res), as.matrix(ref_m_filter))$estF
    predicted_cell_proportions_epidish <- round(100*predicted_cell_proportions_epidish,2)
    epidishavg <- colMeans(predicted_cell_proportions_epidish)
    epidishavg <- data.frame(epidishavg)
    epidishavg$Celltype <- rownames(epidishavg)
    write.table(epidishavg, "./LONGVAR.deconv.epidish.REFsam.set8.avgacrosssamples.tsv", quote=FALSE, row.names=FALSE,sep="\t")

    colnames(predicted_cell_proportions_epidish) <- colnames(ref_m_filter)
    predicted_cell_proportions_epidish <- as.data.frame(predicted_cell_proportions_epidish)
    predicted_cell_proportions_epidish$Samplename <- colnames(res)
    epidish_res <- merge(metainfo, predicted_cell_proportions_epidish, by=c("Samplename"))
    #write.table(epidish_res, "./LONGVAR.deconv.epidish.REFsam.set8.tsv",row.names=FALSE, quote=FALSE, sep="\t")

    predicted_cell_proportions_epidish$method <- "EpiDISH"
    return(predicted_cell_proportions_epidish)
}

get_nnls_deconv_res <- function(metafh, covfiledir, refavg) {
    # Read sample metadata
    metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)

    # Get sample matrix and reference methylation matrix
    meth_all_cpgs_imputed <- get_sample_mat(metafh, covfiledir, refavg)

    # Read reference methylation matrix
    ref_m <- fread(refavg,header=TRUE,sep="\t",data.table=FALSE)
    rownames(ref_m) <- paste0(ref_m$chr,"_",ref_m$start,"_",ref_m$end)
    ref_m_filter = ref_m[rownames(meth_all_cpgs_imputed),-which(names(ref_m) %in% c("chr","start","end"))]
    res = meth_all_cpgs_imputed

    set.seed(9)
    nnls_sol <- function(res, ref_m_filter) { 
        X <- matrix(0, nrow=dim(res)[2], ncol=dim(ref_m_filter)[2]) 
        for(i in 1:ncol(res))
        X[i,] <- coef(nnls::nnls(as.matrix(ref_m_filter),res[,i]))
        X
    }

    X_nnls <- nnls_sol(res,ref_m_filter)
    colnames(X_nnls) <- colnames(ref_m_filter)
    rownames(X_nnls) <- colnames(res)
    X_nnls <- as.data.frame(X_nnls)
    X_nnls$Samplename <- colnames(res)
    nnls_res <- merge(metainfo, X_nnls, by=c("Samplename"))

    X_nnls$method <- "NNLS"
    return(X_nnls)
}

get_emeth_deconv_res <- function(metafh, covfiledir, refavg) {
    # Read sample metadata
    metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)

    # Get sample matrix and reference methylation matrix
    meth_all_cpgs_imputed <- get_sample_mat(metafh, covfiledir, refavg)

    # Read reference methylation matrix
    ref_m <- fread(refavg,header=TRUE,sep="\t",data.table=FALSE)
    rownames(ref_m) <- paste0(ref_m$chr,"_",ref_m$start,"_",ref_m$end)
    ref_m_filter = ref_m[rownames(meth_all_cpgs_imputed),-which(names(ref_m) %in% c("chr","start","end"))]
    res = meth_all_cpgs_imputed

    # Deconvolution using EMeth
    deconv = 'emeth_laplace'
    maximum_nu = 0
    maximum_iter = 100
    predicted_cell_proportions_emeth <- matrix(NA, ncol = dim(res)[2], nrow = dim(ref_m_filter)[2])
    set.seed(9)

    em <- cv.emeth(Y = as.matrix(res[rownames(ref_m_filter),]),eta = c(rep(0, dim(as.matrix(res[rownames(ref_m_filter),]))[2])), mu = as.matrix(ref_m_filter), aber = FALSE, V = 'c',init = "default", family = "laplace", folds=10, nu=maximum_nu, maxiter = maximum_iter, verbose = TRUE)

    predicted_cell_proportions_emeth <- round(100*em[[1]]$rho,2)
    colnames(predicted_cell_proportions_emeth) <- colnames(ref_m_filter)
    predicted_cell_proportions_emeth <- as.data.frame(predicted_cell_proportions_emeth)
    predicted_cell_proportions_emeth$Samplename <- colnames(res)
    emeth_res <- merge(metainfo, predicted_cell_proportions_emeth, by=c("Samplename"))

    predicted_cell_proportions_emeth$method <- "EMeth"
    return(predicted_cell_proportions_emeth)
}

compare_deconv_res <- function(metafh, covfiledir, refavg) {
    # Read sample metadata
    metainfo <- fread(metafh, header=TRUE, sep="\t", data.table=FALSE)
    submeta <- metainfo[,c("SubjectID","Samplename","Year")]

    # Compare EpiDISH and NNLS results
    epidish_res <- get_epidish_deconv_res(metafh, covfiledir, refavg)
    nnls_res <- get_nnls_deconv_res(metafh, covfiledir, refavg)

    # combine results
    combined_res <- rbind(epidish_res, nnls_res)

    combined_res$MEP <- combined_res$`Erythrocyte-progenitors`+ combined_res$`Megakaryocytes`

    # plot comparison
    majorcelltype <- c("Granulocytes","Monocytes","MEP","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes")
    subcomb <- combined_res[,colnames(combined_res) %in% c(majorcelltype,"Samplename","method")]
    combined_res_long <- gather(subcomb, Celltype, Proportion, `B-cells`:`Vascular-endothelium`)
    dt <- merge(submeta, combined_res_long, by=c("Samplename"))
    dt1 <- subset(dt, Year=="T1")
    dt1$Celltype <- factor(dt1$Celltype, levels=majorcelltype,labels=c("granulocytes","monocytes","MEP","NK cells","B cells","T CD4","VE","hepatocytes"))
    dt1$method <- factor(dt1$method, levels=c("EpiDISH", "NNLS"))

    outfig <- paste0(outdir, "/LONGVAR.deconv.T1sam.comparison.epidish.nnls.REFsam.set8.", heytoday, ".pdf")
    pdf(outfig, width=9.8, height=4.6)
    p0 <- ggplot(dt1, aes(x=method, y=Proportion))+geom_line(aes(group=SubjectID),size=0.2,color="gray")+geom_boxplot(aes(fill=method),width=0.4,alpha=0.5,outlier.shape=NA)+geom_point(aes(color=method),size=0.5, fill="white")+commontheme()+theme(axis.text.x=element_text(size=10,color="black",angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(size=10,color="black"),strip.placement="outside",strip.text=element_text(size=10,face="bold"),strip.background=element_rect(color="white"), axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="none",plot.title=element_text(hjust=0.5))+scale_fill_manual(values=c("black","gray50"))+scale_color_manual(values=c("black","gray50"))+labs(x="",y="proportion (%)")+facet_wrap(~Celltype,scales="free_y",nrow=1)+scale_y_continuous(expand=expansion(mult=c(0.085, .1)))
    print(p0)
    dev.off()

    # plot correlation
    epidishdt <- merge(epidish_res, submeta, by=c("Samplename"))
    nnlsdt <- merge(nnls_res, submeta, by=c("Samplename"))
    epidishdt1 <- subset(epidishdt, Year=="T1")
    nnlsdt1 <- subset(nnlsdt, Year=="T1")
    epidishdt1$MEP <- epidishdt1$`Erythrocyte-progenitors` + epidishdt1$`Megakaryocytes`
    nnlsdt1$MEP <- nnlsdt1$`Erythrocyte-progenitors` + nnlsdt1$`Megakaryocytes`
    
    outfig2 <- paste0(outdir, "/LONGVAR.deconv.T1sam.comparison.epidish.nnls.REFsam.set8.correlation.", heytoday, ".pdf")
    pdf(outfig2, width=9.6, height=9.6)
    plot_list <- vector("list", length = length(majorcelltype))
    for (i in seq_along(majorcelltype)) {
        celltype <- majorcelltype[i]
        epidt <- epidishdt1[,c("Samplename", celltype)]
        nnlst <- nnlsdt1[,c("Samplename", celltype)]
        colnames(epidt) <- c("Samplename", "EpiDISH")
        colnames(nnlst) <- c("Samplename", "NNLS")
        merged_dt <- merge(epidt, nnlst, by = "Samplename")
        p2 <- ggplot(merged_dt, aes(x=EpiDISH, y=NNLS))+geom_point(color="black",size=1,shape=16)+commontheme()+theme(axis.title.x=element_text(size=14,color="black",vjust=-0.5),axis.title.y=element_text(size=14,color="black",vjust=1.5))+labs(title=celltype, x="EpiDISH (%)",y="NNLS (%)")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=4.5)
        plot_list[[i]] <- p2
    }
    pp2 <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 3, align = "hv", axis = "tblr")
    print(pp2)
    dev.off()

    # calculate average difference between methods
    epidish_avg <- colMeans(epidishdt1[,majorcelltype], na.rm=TRUE)
    nnls_avg <- colMeans(nnlsdt1[,majorcelltype], na.rm=TRUE)
    avg_diff <- epidish_avg - nnls_avg
    avg_diff <- data.frame(Celltype=names(avg_diff), AvgDiff=avg_diff)
    avg_diff$Celltype <- factor(avg_diff$Celltype, levels=majorcelltype, labels=c("granulocytes","monocytes","MEP","NK cells","B cells","T CD4","VE","hepatocytes"))
}
