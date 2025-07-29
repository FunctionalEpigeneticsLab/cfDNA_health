#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort (size profile and 2mer end motif)

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "reshape2", "ggpubr", "rstatix", "ggpmisc", "janitor", "circlize", "ComplexHeatmap")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
sizestatsfh <- args[2]
sizeprofilefh <- args[3]
kmer <- "2mer"
motifh <- args[4]
sizemotiffh <- args[5]
outdir <- args[6]

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

get_inter_intra_cor <- function(deconvresfh, sizeprofilefh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))

    sizedt0 <- fread(sizeprofilefh,header=TRUE,sep="\t",data.table=FALSE)
    sizedt0[is.na(sizedt0)] <- 0

    sizedt <- subset(sizedt0, Size>=50 & Size<=450)

    sizemat <- sizedt[,! (colnames(sizedt) %in% c("Size"))]
    sizecor <- cor(sizemat,method="kendall")

    flattenCorrMatrix <- function(cormat) {
        ut <- upper.tri(cormat)
        data.frame(Sample1=rownames(cormat)[row(cormat)[ut]], Sample2=rownames(cormat)[col(cormat)[ut]], Kcor=(cormat)[ut])
    }

    lsizemat <- flattenCorrMatrix(sizecor)

    smetainfo <- metainfo[,c("GCcode","SubjectID","Year","Time","Age","Sex","BMI")]
    smetainfo1 <- smetainfo
    smetainfo2 <- smetainfo
    colnames(smetainfo1) <- c("Sample1","SubjectID1","Year1","Time1","Age1","Sex1","BMI1")
    colnames(smetainfo2) <- c("Sample2","SubjectID2","Year2","Time2","Age2","Sex2","BMI2")

    mlsizemat <- merge(lsizemat, smetainfo1, by=c("Sample1"),all.x=TRUE)
    mlsizemat <- merge(mlsizemat, smetainfo2, by=c("Sample2"), all.x=TRUE)

    mlsizemat <- mlsizemat %>% mutate(Age_group1=cut(Age1,breaks=c(19,29,39,49,59,69,100),labels=c("20s","30s","40s","50s","60s","70s+")))
    mlsizemat <- mlsizemat %>% mutate(Age_group2=cut(Age2-1,breaks=c(19,29,39,49,59,69,100),labels=c("20s","30s","40s","50s","60s","70s+")))


    summary(mlsizemat[mlsizemat$SubjectID1==mlsizemat$SubjectID2,]$Kcor)
    summary(mlsizemat[mlsizemat$SubjectID1!=mlsizemat$SubjectID2,]$Kcor)

    wilcox.test(mlsizemat[mlsizemat$SubjectID1==mlsizemat$SubjectID2,]$Kcor,mlsizemat[mlsizemat$SubjectID1!=mlsizemat$SubjectID2,]$Kcor, alternative="greater")

    summary(mlsizemat[(mlsizemat$SubjectID1!=mlsizemat$SubjectID2) & (mlsizemat$Age_group1==mlsizemat$Age_group2),]$Kcor)
    summary(mlsizemat[(mlsizemat$SubjectID1!=mlsizemat$SubjectID2) & (mlsizemat$Age_group1!=mlsizemat$Age_group2),]$Kcor)

    summary(mlsizemat[(mlsizemat$SubjectID1!=mlsizemat$SubjectID2) & (mlsizemat$Sex1==mlsizemat$Sex2),]$Kcor)
    summary(mlsizemat[(mlsizemat$SubjectID1!=mlsizemat$SubjectID2) & (mlsizemat$Sex1!=mlsizemat$Sex2),]$Kcor)
}

plot_sizeprofile_heatmap <- function(deconvresfh, sizeprofilefh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))
    smetainfo <- metainfo[,c("GCcode","SubjectID","Year","Time","Age","Sex","BMI")]

    sizedt0 <- fread(sizeprofilefh,header=TRUE,sep="\t",data.table=FALSE)
    sizedt0[is.na(sizedt0)] <- 0

    sizedt <- subset(sizedt0, Size>=50 & Size<=450)
    longsizedt <- gather(sizedt, GCcode, Freq, -Size)
    wsizedt <- spread(longsizedt, Size, Freq)

    wsizemat <- wsizedt[,!(colnames(wsizedt) %in% c("GCcode"))]
    res.pca <- prcomp(wsizemat, scale = TRUE)
    summary(res.pca)
    
    wsizemetadt <- merge(wsizedt, smetainfo, by=c("GCcode"))
    wsizemetadt1 <- subset(wsizemetadt, Year=="Year1")
    wsizemetadt2 <- subset(wsizemetadt, Year=="Year2")

    mat1 <- as.matrix(wsizemetadt1[,! (colnames(wsizemetadt1) %in% c("SubjectID","GCcode","Year","Time","Age","Sex","BMI"))])
    row.names(mat1) <- wsizemetadt1$SubjectID

    col_fun1 = colorRamp2(c(20, 30, 40, 50, 60, 70, 90), c("#f7fcf5","#e5f5e0", "#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45"))
    col_fun2 = colorRamp2(c(16.5, 37), c("#FCFBFD","#3F007D"))

    panel_fun = function(index, nm) {
        grid.rect(gp=gpar(fill=c("#737373")))
        grid.text(paste0(index+49,"-",index+149, " bp"), 0.5, 0.5, gp=gpar(col="white",fontsize = 12))
    }
    align_to = list("50-150"=1:100, "151-250"=101:200, "251-350"=201:300, "351-450"=301:400)
    ha <- rowAnnotation(Age=wsizemetadt1$Age, Sex=wsizemetadt1$Sex, BMI=wsizemetadt1$BMI,simple_anno_size=unit(1, "cm"),col=list(Age=col_fun1, Sex=c("Female"="#EA7317","Male"="#2364AA"), BMI=col_fun2))
    va <- HeatmapAnnotation(foo=anno_block(align_to=align_to, panel_fun = panel_fun),bar=anno_mark(at=c(1,96,118, 138, 231, 281, 371),labels=c("50bp","145bp","167bp","187bp","280bp","330bp","420bp"), which="column", side="bottom"))
    va2 <- HeatmapAnnotation(Density=anno_density(mat1,type="violin",gp=gpar(col="black")))
    #va2 <- HeatmapAnnotation(Density=anno_lines(mat1))

    outfig0 <- paste0(outdir,"/LONGVAR.sizeprofile.Year1.heatmap.",heytoday,".pdf")
    pdf(outfig0,height=7,width=11)
    p0 <- ComplexHeatmap::Heatmap(mat1,name="mat",right_annotation=ha,bottom_annotation=va,top_annotation=va2,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,column_names_gp=gpar(fontsize=7))
    print(p0)
    dev.off()

    mat2 <- as.matrix(wsizemetadt2[,! (colnames(wsizemetadt2) %in% c("SubjectID","GCcode","Year","Time","Age","Sex","BMI"))])
    row.names(mat2) <- wsizemetadt2$SubjectID
    va2 <- HeatmapAnnotation(Density=anno_density(mat2,type="violin",gp=gpar(col="black")))
    outfig1 <- paste0(outdir,"/LONGVAR.sizeprofile.Year2.heatmap.",heytoday,".pdf")
    pdf(outfig1,height=7,width=11)
    p1 <- ComplexHeatmap::Heatmap(mat2,name="mat",right_annotation=ha,bottom_annotation=va,top_annotation=va2,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,column_names_gp=gpar(fontsize=7))
    print(p1)
    dev.off()
}

plot_sizeprofile_corr <- function(deconvresfh, sizeprofilefh, sizestatsfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))

    sizedt0 <- fread(sizeprofilefh,header=TRUE,sep="\t",data.table=FALSE)
    sizedt0[is.na(sizedt0)] <- 0

    sizedt <- subset(sizedt0, Size>=50 & Size<=450)
    longsizedt <- gather(sizedt, GCcode, Freq, -Size)
    wsizedt <- spread(longsizedt, Size, Freq)

    smetainfo <- metainfo[,c("GCcode","Year")]
    wsizemetadt <- merge(wsizedt, smetainfo, by=c("GCcode"))
    wsizemetadt1 <- subset(wsizemetadt, Year=="Year1")
    wsizemetadt2 <- subset(wsizemetadt, Year=="Year2")

    wsizemetamat1 <- wsizemetadt1[,!(colnames(wsizemetadt1) %in% c("GCcode","Year"))]
    sizeselfcor1 <- cor(wsizemetamat1)
    sizeselfspcor1 <- cor(wsizemetamat1,method="pearson")

    panel_fun = function(index, nm) {
        grid.rect(gp=gpar(fill=c("#737373")))
        grid.text(paste0(index+49,"-",index+149, " bp"), 0.5, 0.5, gp=gpar(col="white",fontsize = 12))
    }
    align_to = list("50-150"=1:100, "151-250"=101:200, "251-350"=201:300, "351-450"=301:400)
    va <- HeatmapAnnotation(bar=anno_mark(at=c(1,96,118, 138, 231, 281, 371),labels=c("50bp","145bp","167bp","187bp","280bp","330bp","420bp"), which="column", side="bottom"))

    panel_fun1 = function(index, nm) {
        grid.rect(gp=gpar(fill=c("#737373")))
        grid.text(paste0(index+49,"-",index+149, " bp"), 0.5, 0.5, gp=gpar(col="white",fontsize = 12),rot=270)
    }
    ha <- rowAnnotation(bar1=anno_mark(at=c(1,96,118, 138, 231, 281, 371),labels=c("50bp","145bp","167bp","187bp","280bp","330bp","420bp"), which="row", side="left"))
    col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#f7fbff", "#c6dbef", "#6baed6","#2171b5", "#08306b"))

    #mysizegroup <- c(rep(1,97),rep(2,41),rep(3,146),rep(4,117))
    #50-145, 146-166, 167-187, 188-279, 280-329, 330-420, 421-250bp
    mysizegroup <- c(rep(1,96),rep(2,21),rep(3,21),rep(4,92),rep(5,50),rep(6,91),rep(7,30))

    outfig0 <- paste0(outdir,"/LONGVAR.sizeprofile.Year1.sizepearcor.",heytoday,".pdf")
    pdf(outfig0,height=8,width=8.4)
    p0 <- ComplexHeatmap::Heatmap(sizeselfspcor1,name="pearson\ncorr",left_annotation=ha, bottom_annotation=va,cluster_rows=FALSE,cluster_columns=FALSE,row_split=mysizegroup,row_title=c("sub-\nnucleosome","","mono-\nnucleosome","","","di-\nnucleosome",""),row_title_gp=gpar(fontsize=11, fontface = "bold"),column_split=mysizegroup,column_title=c("sub-\nnucleosome","","mono-\nnucleosome","","","di-\nnucleosome",""),column_title_gp=gpar(fontsize=11, fontface = "bold"),column_title_side=c("bottom"),show_row_names=FALSE,show_column_names=FALSE,column_names_gp=gpar(fontsize=7),col = col_fun)
    print(p0)
    dev.off()

    wsizemetamat2 <- wsizemetadt2[,!(colnames(wsizemetadt2) %in% c("GCcode","Year"))]
    sizeselfcor2 <- cor(wsizemetamat2)
    sizeselfspcor2 <- cor(wsizemetamat2,method="pearson")
    va0 <- HeatmapAnnotation(Density=anno_density(wsizemetamat2,gp=gpar(col="black")))

    outfig2 <- paste0(outdir,"/LONGVAR.sizeprofile.Year2.sizepearcor.",heytoday,".pdf")
    pdf(outfig2,height=8,width=8.4)
    p0 <- ComplexHeatmap::Heatmap(sizeselfspcor2,name="pearson\ncorr",left_annotation=ha, bottom_annotation=va,cluster_rows=FALSE,cluster_columns=FALSE,row_split=mysizegroup,row_title=c("sub-\nnucleosome","","mono-\nnucleosome","","","di-\nnucleosome",""),row_title_gp=gpar(fontsize=11, fontface = "bold"),column_split=mysizegroup,column_title=c("sub-\nnucleosome","","mono-\nnucleosome","","","di-\nnucleosome",""),column_title_gp=gpar(fontsize=11, fontface = "bold"),column_title_side=c("bottom"),show_row_names=FALSE,show_column_names=FALSE,column_names_gp=gpar(fontsize=7),col = col_fun)
    print(p0)
    dev.off()

    sizestats <- fread(sizestatsfh,header=TRUE,sep="\t",data.table=FALSE)
    dt0 <- merge(metainfo, sizestats, by=c("GCcode"))
    dt1 <- subset(dt0, Year=="Year1")
    dt2 <- subset(dt0, Year=="Year2")

    outfig0 <- paste0(outdir,"/LONGVAR.sizerange.Year1.corr.",heytoday,".pdf")
    pdf(outfig0,height=4.2,width=4.2)
    p0 <- ggplot(dt1,aes(x=Sizearea_sub2, y=Sizearea_sub))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(146-166bp)",y="Cumulative frequencies\n(50-145bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt1,aes(x=Sizearea_mo, y=Sizearea_sub))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(50-145bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt1,aes(x=Sizearea_mo, y=Sizearea_sub2))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(146-166bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt1,aes(x=Sizearea_mo, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt1$Sizearea_di)+0.01)
    print(p0)
    p0 <- ggplot(dt1,aes(x=Sizearea_sub, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(50-145bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt1$Sizearea_di)+0.01)
    print(p0)
    p0 <- ggplot(dt1,aes(x=Sizearea_sub2, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(146-166bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt1$Sizearea_di)+0.01)
    print(p0)
    dev.off()

    outfig0 <- paste0(outdir,"/LONGVAR.sizerange.Year2.corr.",heytoday,".pdf")
    pdf(outfig0,height=4.2,width=4.2)
    p0 <- ggplot(dt2,aes(x=Sizearea_sub2, y=Sizearea_sub))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(146-166bp)",y="Cumulative frequencies\n(50-145bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt2,aes(x=Sizearea_mo, y=Sizearea_sub))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(50-145bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt2,aes(x=Sizearea_mo, y=Sizearea_sub2))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(146-166bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")
    print(p0)
    p0 <- ggplot(dt2,aes(x=Sizearea_mo, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(167-187bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt2$Sizearea_di)+0.01)
    print(p0)
    p0 <- ggplot(dt2,aes(x=Sizearea_sub, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(50-145bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt2$Sizearea_di)+0.01)
    print(p0)
    p0 <- ggplot(dt2,aes(x=Sizearea_sub2, y=Sizearea_di))+geom_point(color="black",fill="white",size=1.5,shape=21)+commontheme()+labs(x="Cumulative frequencies\n(146-166bp)",y="Cumulative frequencies\n(330-420bp)")+geom_smooth(color="darkblue",method=lm, se=FALSE,linetype="solid")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson")+ylim(0,max(dt2$Sizearea_di)+0.01)
    print(p0)
    dev.off()
}

plot_endmotif <- function(deconvresfh, motifh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))
    smetainfo <- metainfo[,c("GCcode","Year")]

    motifdt <- fread(motifh,header=TRUE,sep="\t",data.table=FALSE)
    mdt <- merge(smetainfo,motifdt,by=c("GCcode"))
    lmdt <- gather(mdt, endmotif, motifreq, AA:TT)
    lmdt1 <- subset(lmdt, Year=="Year1")
    lmdt2 <- subset(lmdt, Year=="Year2")

    outfig0 <- paste0(outdir,"/LONGVAR.endmotif.dotboxplot.",heytoday,".pdf")
    pdf(outfig0,height=5.2,width=9.2)
    p1 <- ggplot(lmdt1, aes(x=endmotif, y=motifreq))+geom_boxplot(color="black",width=0.7,outlier.shape=NA)+geom_dotplot(color="black",binwidth=0.1,binaxis='y', stackdir='center', dotsize=0.6,position=position_dodge(0.5),alpha=0.5)+commontheme()+theme(legend.position="none")+labs(title="Year1",x="Dinucleotide end motif",y="Motif frequency (%)")
    print(p1)
    p2 <- ggplot(lmdt2, aes(x=endmotif, y=motifreq))+geom_boxplot(color="black",width=0.7,outlier.shape=NA)+geom_dotplot(color="black",binwidth=0.1,binaxis='y', stackdir='center', dotsize=0.6,position=position_dodge(0.5),alpha=0.5)+commontheme()+theme(legend.position="none")+labs(title="Year2",x="Dinucleotide end motif",y="Motif frequency (%)")
    print(p2)
    dev.off()

    ###Sex difference
    library(ggh4x)
    smetainfo <- metainfo[,c("GCcode","Year","Sex","Age","BMI")]
    mdt <- merge(smetainfo, motifdt, by=c("GCcode"))
    sexmotif <- c("AC","CC","TG","TT")
    sexdt <- mdt[,colnames(mdt) %in% c("GCcode","Year","Sex","Age","BMI",sexmotif)]
    lsexdt <- gather(sexdt, endmotif, motifreq, AC:TT)
    sexscales_y <- list(
       `AC` = scale_y_continuous(limits = c(3.5, 8), breaks=c(4, 5, 6)),
       `CC` = scale_y_continuous(limits = c(10.5, 21), breaks=c(12, 14, 16)),
       `TG` = scale_y_continuous(limits = c(7.5, 14), breaks=c(8, 9, 10)),
       `TT` = scale_y_continuous(limits = c(1.5,6), breaks=c(2,3,4)))
    outfig1 <- paste0(outdir,"/LONGVAR.endmotif.Sexdiff.",heytoday,".pdf")
    pdf(outfig1, width=5.2, height=7.2)
    stat.test <- lsexdt %>% group_by(Year, endmotif) %>% wilcox_test(motifreq ~ Sex, p.adjust.method="bonferroni") %>% add_xy_position(x = "Sex")
    p1 <- ggplot(lsexdt, aes(x=Sex, y=motifreq))+geom_boxplot(lwd=1, aes(color=Sex),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Sex),position=position_jitter(0.2))+commontheme()+theme(plot.title = element_text(size=16,face="bold"),legend.position="none")+labs(x="",y="Motif frequency (%)")+stat_pvalue_manual(stat.test,label="p",size=4,step.increase=0.01,vjust=0,tip.length=0)+scale_color_manual(values=c("#EA7317","#2364AA"))+facet_grid(endmotif~Year,scales="free_y")+scale_y_continuous(expand=expansion(mult = c(0.08, 0.1)))+ggh4x::facetted_pos_scales(y=sexscales_y)
    print(p1)
    dev.off()

    ###corr BMI
    bmimotif <- c("AC","CG","GG","TG")
    bmidt <- mdt[,colnames(mdt) %in% c("GCcode","Year","Sex","Age","BMI",bmimotif)]
    lbmidt <- gather(bmidt, endmotif, motifreq, AC:TG)
    outfig2 <- paste0(outdir,"/LONGVAR.endmotif.corBMI.",heytoday,".pdf")
    pdf(outfig2, width=5.2, height=7.2)
    p2 <- ggplot(lbmidt, aes(x=BMI,y=motifreq))+geom_point(color="black",size=1.5,shape=16)+commontheme()+theme(legend.position="none")+labs(x="BMI",y="Motif frequency (%)")+scale_y_continuous(expand=expansion(mult = c(0.1, 0.2)))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="right",label.y="top",use_label(c("R", "P")),method="pearson",size=4)+facet_grid(endmotif~Year, scales="free_y")
    print(p2)
    dev.off()

    ###corr Age
    agemotif <- c("AG","AT","TG","TT")
    agemdt <- mdt[,colnames(mdt) %in% c("GCcode","Year","Sex","Age","BMI",agemotif)]
    lagemdt <- gather(agemdt, endmotif, motifreq, AG:TT)
    outfig3 <- paste0(outdir,"/LONGVAR.endmotif.corAge.",heytoday,".pdf")
    pdf(outfig3, width=5.2, height=7.2)
    p3 <- ggplot(lagemdt, aes(x=Age,y=motifreq))+geom_point(color="black",size=1.5,shape=16)+commontheme()+theme(legend.position="none")+labs(x="Age",y="Motif frequency (%)")+scale_y_continuous(expand=expansion(mult = c(0.1, 0.2)))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="right",label.y="top",use_label(c("R", "P")),method="pearson",size=4)+facet_grid(endmotif~Year, scales="free_y")
    print(p3)
    dev.off()
}

GetCorrMat <- function(smdt1,mmdt1) {
    cormat <- data.frame()
    for (mysize in 50:450) {
      for (mymotif in unique(mmdt1$Motif)) {
         print(paste0("Size: ", mysize, ";Motif: ", mymotif))
         subsmdt1 <- subset(smdt1, Size==mysize)
         submmdt1 <- subset(mmdt1, Motif==mymotif)
         subdt <- merge(subsmdt1,submmdt1,by=c("GCcode"))
         mycor <- cor(subdt$Size_freq,subdt$Motif_freq, method="spearman")
         samcor <- data.frame(Size=mysize, Motif=mymotif, Scorr=mycor)
         cormat <- rbind(cormat, samcor)
      }
    }
    return(cormat)
}

GetConcMotifCorr <- function(motifcormetadt,inYear) {
  myconcmotifcor <- data.frame(Motif=character(),Corr=numeric())
  for (i in 4:ncol(motifcormetadt)) {
      curcorr <- cor(motifcormetadt$Adj_Unique_Molecules,motifcormetadt[,i],method="spearman")
      curdt <- data.frame(Motif=colnames(motifcormetadt)[i],Corr=curcorr)
      myconcmotifcor <- rbind(myconcmotifcor,curdt)
  }
  myconcmotifcor$Year <- inYear
  return(myconcmotifcor)
}

plotsizecormotif <- function(mdt,metadt,inYear,outdir,kmer) {
    mdt1 <- subset(mdt,Year==inYear)
    smdt1 <- mdt1 %>% group_by(GCcode, Size) %>% summarise(Size_freq=sum(Size_Motif_Freq))
    mmdt1 <- mdt1 %>% group_by(GCcode, Motif) %>% summarise(Motif_freq=sum(Size_Motif_Freq))
    totabund <- sum(mmdt1$Motif_freq)
    rankmotif1 <- mmdt1 %>% group_by(Motif) %>% summarise(Abundance=sum(Motif_freq)/totabund) %>% mutate(rank=rank(-Abundance))
    orderedmotif1 <- rankmotif1[order(rankmotif1$rank),]$Motif

    sizemotifcormat1 <- GetCorrMat(smdt1,mmdt1)
    wmmdt1 <- spread(mmdt1, Motif, Motif_freq)
    concmetainfo <- metadt[,c("GCcode","Year","Adj_Unique_Molecules")]
    motifcormetadt1 <- merge(concmetainfo, wmmdt1, by=c("GCcode"))
    t1mcor <- GetConcMotifCorr(motifcormetadt1,inYear)

    sizemotifcordt1 <- spread(sizemotifcormat1, Size, Scorr)
    sizemotifcordt1 <- merge(rankmotif1, sizemotifcordt1, by=c("Motif"))
    concsizemotifcordt1 <- merge(t1mcor, sizemotifcordt1, by=c("Motif"))
    #concsizemotifcordt1 <- concsizemotifcordt1[match(orderedmotif1,concsizemotifcordt1$Motif),]
    smcormat1 <- as.matrix(concsizemotifcordt1[,! (colnames(concsizemotifcordt1) %in% c("Motif","Corr","Year","Abundance","rank"))])
    row.names(smcormat1) <- concsizemotifcordt1$Motif

    mincorval <- min(t1mcor$Corr)
    maxcorval <- max(t1mcor$Corr)
    corval2 <- mincorval+(maxcorval-mincorval)/6
    corval3 <- mincorval+(maxcorval-mincorval)/3
    corval4 <- maxcorval-(maxcorval-mincorval)/3
    corval5 <- maxcorval-(maxcorval-mincorval)/6
    col_fun0 = colorRamp2(c(mincorval,corval2,corval3,0,corval4,corval5,maxcorval),c("#542788","#998ec3","#d8daeb","#f7f7f7","#fee0b6","#f1a340","#b35806"))

    minabdval <- min(concsizemotifcordt1$Abundance)
    maxabdval <- max(concsizemotifcordt1$Abundance)
    midabdval <- minabdval+(maxabdval-minabdval)/2
    col_fun1 = colorRamp2(c(minabdval, midabdval, maxabdval), c("#4d9221","#f7f7f7","#c51b7d"))
    panel_fun = function(index, nm) {
        grid.rect(gp=gpar(fill=c("#737373")))
        grid.text(paste0(nm, " bp"), 0.5, 0.5, gp=gpar(col="white",fontsize = 12))
    }
    align_to = list("50-150"=1:101, "151-250"=102:201, "251-350"=202:301, "351-450"=302:401)
    va <- HeatmapAnnotation(foo=anno_block(align_to=align_to, panel_fun = panel_fun),bar=anno_mark(at=c(1,96,118, 138, 231, 281, 371),labels=c("50bp","145bp","167bp","187bp","280bp","330bp","420bp"), which="column", side="bottom"))
    #ha <- rowAnnotation(CorrConc=concsizemotifcordt1$Corr, Abundance=concsizemotifcordt1$Abundance,foo=anno_mark(at=c(2,1,4,64,21,7,29),labels=c("CCCA","CCTG","CCAG","TAAA","TGAG","CCCT","CCAA")),simple_anno_size=unit(0.8, "cm"),col=list(CorrConc=col_fun0,Abundance=col_fun1))
    #ha <- rowAnnotation(CorrConc=concsizemotifcordt1$Corr, Abundance=concsizemotifcordt1$Abundance,simple_anno_size=unit(0.8, "cm"),col=list(CorrConc=col_fun0,Abundance=col_fun1),annotation_legend_param=list(CorrConc=list(at=c(round(min(concsizemotifcordt1$Corr),2), 0, round(max(concsizemotifcordt1$Corr),2)))))
    ha <- rowAnnotation(Abundance=concsizemotifcordt1$Abundance,simple_anno_size=unit(0.8, "cm"),col=list(Abundance=col_fun1))

    outfig0 <- paste0(outdir,"/LONGVAR.size_motif.",kmer,".profile.",inYear,".scorr.heatmap.",heytoday,".pdf")
    pdf(outfig0,height=7,width=10.2)
    p0 <- ComplexHeatmap::Heatmap(smcormat1,name="corr_mat",heatmap_legend_param=list(at=c(round(min(smcormat1),2), 0, round(max(smcormat1),2))),cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE, right_annotation=ha, bottom_annotation=va, column_names_gp=gpar(fontsize=12, fontface="bold"),row_names_gp=gpar(fontsize=12,fontface="bold"))
    print(p0)
    dev.off()
}

plot_cor_size_motif(sizemotiffh, deconvresfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    smetadt <- metainfo[,colnames(metainfo) %in% c("GCcode","Year")]

    sizemotif <- fread(sizemotiffh,header=TRUE,sep="\t",data.table=FALSE)
    mdt <- merge(sizemotif,smetadt,by=c("GCcode"),all.x=TRUE)

    mdt1 <- subset(mdt,Year=="T1")
    mdt2 <- subset(mdt,Year=="T2")

    plotsizecormotif(mdt1,metainfo,"T1",outdir,kmer)
    plotsizecormotif(mdt2,metainfo,"T2",outdir,kmer)
}
