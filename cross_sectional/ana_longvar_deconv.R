#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort (deconvolution results)

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "reshape2", "ggpubr", "rstatix", "ggpmisc", "scales", "tidyr", "plyr", "GGally")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
wbcfh <- args[2]
agemefh <- args[3]
outdir <- args[4]


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

compare_wbc_deconv <- function(deconvresfh, wbcfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    
    wbc <- fread(wbcfh,header=TRUE,sep="\t",data.table=FALSE)
    wbc$PLymphocytes <- wbc$`T-cells-CD4`+wbc$`T-cells-CD8`+wbc$`B-cells`+wbc$`NK-cells`
    subwbc <- wbc[,c("SubjectID","Samplename","Age","Sex","BMI","Granulocytes","Neutrophils_Per","PLymphocytes","Lymphocytes_Per","Monocytes.x","Monocytes_Per")]

    wbc$Granulocytes_blood <- 100*wbc$Granulocytes/(wbc$Granulocytes+wbc$PLymphocytes+wbc$Monocytes.x)
    wbc$PLymphocytes_blood <- 100*wbc$PLymphocytes/(wbc$Granulocytes+wbc$PLymphocytes+wbc$Monocytes.x)
    wbc$Monocytes_blood <- 100*wbc$Monocytes.x/(wbc$Granulocytes+wbc$PLymphocytes+wbc$Monocytes.x)

    outfig0 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.relblood.corWBCprop.",heytoday,".pdf")
    pdf(outfig0, height=4.5, width=4.5)
    xylower <- min(wbc$Granulocytes_blood, wbc$Neutrophils_Per)-2
    xyupper <- max(wbc$Granulocytes_blood, wbc$Neutrophils_Per)+5
    p0 <- ggplot(wbc, aes(x=Granulocytes_blood, y=Neutrophils_Per))+geom_abline(intercept=0,linetype="dashed",color="gray")+geom_point(color="black",size=3,shape=16)+commontheme()+labs(x="% Granulocytes in plasma blood cells",y="% Neutrophils in WBC")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=4.5)+xlim(xylower, xyupper)+ylim(xylower, xyupper)
    print(p0)
    xylower <- min(wbc$PLymphocytes_blood, wbc$Lymphocytes_Per)-2
    xyupper <- max(wbc$PLymphocytes_blood, wbc$Lymphocytes_Per)+2
    p0 <- ggplot(wbc, aes(x=PLymphocytes_blood, y=Lymphocytes_Per))+geom_abline(intercept=0,linetype="dashed",color="gray")+geom_point(color="black",size=3,shape=16)+commontheme()+labs(x="% Lymphocytes in plasma blood cells",y="% Lymphocytes in WBC")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=4.5)+xlim(xylower, xyupper)+ylim(xylower, xyupper)
    print(p0)
    xylower <- min(wbc$Monocytes_blood, wbc$Monocytes_Per)-2
    xyupper <- max(wbc$Monocytes_blood, wbc$Monocytes_Per)+2
    p0 <- ggplot(wbc, aes(x=Monocytes_blood, y=Monocytes_Per))+geom_abline(intercept=0,linetype="dashed",color="gray")+geom_point(color="black",size=3,shape=16)+commontheme()+labs(x="% Monocytes in plasma blood cells",y="% Monocytes in WBC")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")),method="pearson",size=4.5)+xlim(xylower, xyupper)+ylim(xylower, xyupper)
    print(p0)
    dev.off()

    subwbc <- wbc[,c("SubjectID","Samplename","Age","Sex","BMI","Granulocytes_blood","Neutrophils_Per","PLymphocytes_blood","Lymphocytes_Per","Monocytes_blood","Monocytes_Per")]
    longsubwbc <- melt(subwbc, id.vars=c("SubjectID","Samplename","Age","Sex","BMI"),value.name="Proportion",variable.name="Celltype")

    outfig0 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.relblood.WBCprop_boxplot.",heytoday,".pdf")
    pdf(outfig0,height=4.8,width=4)
    longsubwbc0 <- subset(longsubwbc,Celltype=="Granulocytes_blood" | Celltype=="Neutrophils_Per")
    longsubwbc0$Celltype <- factor(longsubwbc0$Celltype)
    stat.test0 <- longsubwbc0 %>% pairwise_wilcox_test(Proportion ~ Celltype, paired=TRUE, p.adjust.method="bonferroni") %>% add_xy_position(x = "Celltype")
    p0 <- ggplot(longsubwbc0, aes(x=Celltype,y=Proportion))+geom_boxplot(width=0.5,position=position_dodge(1))+geom_point(shape=21, size=1.5)+geom_line(aes(group=SubjectID),color="gray")+commontheme()+theme(legend.position="none")+labs(x="",y="Proportion (%)")+scale_x_discrete(breaks=c("Granulocytes_blood","Neutrophils_Per"),labels=c("Plasma\ngranulocytes","WBC\nneutrophils"))+stat_pvalue_manual(stat.test0,label="p",size=4.5,step.increase=0.01,vjust=0,tip.length=0)+scale_y_continuous(expand = c(0.1, 0.25))
    print(p0)
    longsubwbc1 <- subset(longsubwbc,Celltype=="PLymphocytes_blood" | Celltype=="Lymphocytes_Per")
    longsubwbc1$Celltype <- factor(longsubwbc1$Celltype)
    stat.test1 <- longsubwbc1 %>% pairwise_wilcox_test(Proportion ~ Celltype,p.adjust.method="bonferroni") %>% add_xy_position(x = "Celltype")
    p1 <- ggplot(longsubwbc1, aes(x=Celltype,y=Proportion))+geom_boxplot(width=0.5,position=position_dodge(1))+geom_point(shape=21, size=1.5)+geom_line(aes(group=SubjectID),color="gray")+commontheme()+theme(legend.position="none")+labs(x="",y="Proportion (%)")+scale_x_discrete(breaks=c("PLymphocytes_blood","Lymphocytes_Per"),labels=c("Plasma\nlymphocytes","WBC\nlymphocytes"))+stat_pvalue_manual(stat.test1,label="p",size=4.5,step.increase=0.01,vjust=0,tip.length=0)+scale_y_continuous(expand = c(0.1, 0.25))
    print(p1)
    longsubwbc2 <- subset(longsubwbc,Celltype=="Monocytes_blood" | Celltype=="Monocytes_Per")
    longsubwbc2$Celltype <- factor(longsubwbc2$Celltype)
    stat.test2 <- longsubwbc2 %>% pairwise_wilcox_test(Proportion ~ Celltype,p.adjust.method="bonferroni") %>% add_xy_position(x = "Celltype")
    p2 <- ggplot(longsubwbc2, aes(x=Celltype,y=Proportion))+geom_boxplot(width=0.5,position=position_dodge(1))+geom_point(shape=21, size=1.5)+geom_line(aes(group=SubjectID),color="gray")+commontheme()+theme(legend.position="none")+labs(x="",y="Proportion (%)")+scale_x_discrete(breaks=c("Monocytes_blood","Monocytes_Per"),labels=c("Plasma\nmonocytes","WBC\nmonocytes"))+stat_pvalue_manual(stat.test2,label="p",size=4.5,step.increase=0.01,vjust=0,tip.length=0)+scale_y_continuous(expand = c(0.1, 0.25))
    print(p2)
    dev.off()
}

eval_cross_deconv <- function(deconvresfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))

    dt1 <- subset(metainfo,Year=="Year1")
    dt2 <- subset(metainfo,Year=="Year2")
    
    selectcelltype <- c("Granulocytes","Monocytes","Megakaryocytes","Erythrocyte-progenitors","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes")

    cv <- function(x) {
        cval <- sd(x)/mean(x)
        cval <- round(cval, digits=2)
        return(cval)
    }
    cvmajorcell1 <- data.frame(apply(dt1[,colnames(dt1) %in% selectcelltype],2,cv))
    colnames(cvmajorcell1) <- "CV"
    cvmajorcell1$Celltype <- rownames(cvmajorcell1)
    pmax <- function(x) 0.85*max(x)
    cvmajorcell1$Proportion <- apply(dt1[,colnames(dt1) %in% selectcelltype],2,pmax)
    cvmajorcell1$Celltype <- factor(cvmajorcell1$Celltype, levels=selectcelltype)
    cvmajorcell2 <- apply(dt2[,colnames(dt2) %in% selectcelltype],2,cv)
    cvmajorcell2 <- data.frame(cvmajorcell2)
    colnames(cvmajorcell2) <- "CV"
    cvmajorcell2$Celltype <- rownames(cvmajorcell2)
    cvmajorcell2$Proportion <- apply(dt2[,colnames(dt2) %in% selectcelltype],2,pmax)
    cvmajorcell2$Celltype <- factor(cvmajorcell2$Celltype, levels=selectcelltype)


    sdt1 <- dt1[,c("SubjectID","Samplename","Year","Time","Age","Sex","Height","Weight","BMI",selectcelltype)]
    sdt1long <- melt(sdt1, id.vars=c("SubjectID","Samplename","Year","Time","Age","Sex","Height","Weight","BMI"),value.name="Proportion",variable.name="Celltype")

    sdt2 <- dt2[,c("SubjectID","Samplename","Year","Time","Age","Sex","Height","Weight","BMI",selectcelltype)]
    sdt2long <- melt(sdt2, id.vars=c("SubjectID","Samplename","Year","Time","Age","Sex","Height","Weight","BMI"),value.name="Proportion",variable.name="Celltype")

    ###boxplot
    outfig0 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.Majorcellprop_boxplot.",heytoday,".pdf")
    pdf(outfig0,height=4.5,width=5.8)
    p0 <- ggplot(sdt1long, aes(x=Celltype,y=Proportion))+geom_boxplot(color="black",width=0.6,outlier.shape=1)+commontheme()+theme(axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="none")+labs(title="Year1",x="Cell type",y="Proportion (%)")
    print(p0)
    p0 <- ggplot(sdt2long, aes(x=Celltype,y=Proportion))+geom_boxplot(color="black",width=0.6,outlier.shape=1)+commontheme()+theme(axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position="none")+labs(title="Year2",x="Cell type",y="Proportion (%)")
    print(p0)
    dev.off()

    ###density
    manucelltype <- factor(c("Granulocytes","Monocytes","Megakaryocytes","Erythrocyte-progenitors","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes","Epithelial cells","Others"))
    cols <- rainbow(length(manucelltype))
    names(cols) <- manucelltype

    majorcols <- cols[c("Granulocytes","Monocytes","Megakaryocytes","Erythrocyte-progenitors","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes")]
    names(majorcols) <- selectcelltype

    outfig1 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.Majorcellprop_density.",heytoday,".pdf")
    pdf(outfig1,height=7.2,width=9.6)
    p0 <- ggplot(sdt1long, aes(x=Proportion,color=Celltype))+geom_histogram(aes(y=..density..),fill="white")+geom_density(linewidth=1)+commontheme()+theme(legend.position="none")+labs(title="Year1",x="Proportion (%)",y="Density")+scale_color_manual(values=majorcols)+geom_text(data=cvmajorcell1, x=Inf, hjust=1.2, y=Inf, vjust=1.5, color="black", size=4.5, aes(label=paste0("CV: ", CV)))+facet_wrap(~Celltype,scales="free",labeller = label_wrap_gen())
    print(p0)
    p0 <- ggplot(sdt2long, aes(x=Proportion,color=Celltype))+geom_histogram(aes(y=..density..),fill="white")+geom_density(linewidth=1)+commontheme()+theme(legend.position="none")+labs(title="Year2",x="Proportion (%)",y="Density")+scale_color_manual(values=majorcols)+geom_text(data=cvmajorcell2, x=Inf, hjust=1.2, y=Inf, vjust=1.5, color="black", size=4.5, aes(label=paste0("CV: ", CV)))+facet_wrap(~Celltype,scales="free",labeller = label_wrap_gen())
    print(p0)
    dev.off()

    ###Sex difference
    sexcelltype <- c("Megakaryocytes","Erythrocyte-progenitors","T-cells-CD4")
    outfig3 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.Year1.specelltype.Sexdiff.",heytoday,".pdf")
    pdf(outfig3,height=4.8,width=4)
    for (i in sexcelltype) {
        print(i)
        subdt <- subset(sdt1long,Celltype==i)
        stat.test <- subdt %>% group_by(Year) %>% wilcox_test(Proportion ~ Sex, p.adjust.method="bonferroni") %>% add_xy_position(x = "Sex")
        p3 <- ggplot(subdt, aes(x=Sex, y=Proportion))+geom_boxplot(lwd=1, aes(color=Sex),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Sex),position=position_jitter(0.2))+commontheme()+theme(plot.title = element_text(size=16,face="bold"),legend.position="none")+labs(title=i, x="",y="Proportion (%)")+scale_y_continuous(expand=expansion(mult = c(0.05, 0.1)))+stat_pvalue_manual(stat.test,label="p",size=4.5,step.increase=0.01,vjust=0,tip.length=0)+scale_color_manual(values=c("#EA7317","#2364AA"))
        print(p3)
    }
    dev.off()

    outfig4 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.Year2.specelltype.Sexdiff.",heytoday,".pdf")
    pdf(outfig4,height=4.8,width=4)
    for (i in sexcelltype) {
        print(i)
        subdt <- subset(sdt2long,Celltype==i)
        stat.test <- subdt %>% group_by(Year) %>% wilcox_test(Proportion ~ Sex, p.adjust.method="bonferroni") %>% add_xy_position(x = "Sex")
        p4 <- ggplot(subdt, aes(x=Sex, y=Proportion))+geom_boxplot(lwd=1, aes(color=Sex),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Sex),position=position_jitter(0.2))+commontheme()+theme(plot.title = element_text(size=16,face="bold"),legend.position="none")+labs(title=i, x="",y="Proportion (%)")+scale_y_continuous(expand=expansion(mult = c(0.05, 0.1)))+stat_pvalue_manual(stat.test,label="p",size=4.5,step.increase=0.01,vjust=0,tip.length=0)+scale_color_manual(values=c("#EA7317","#2364AA"))
        print(p4)
    }
    dev.off()

    ###corr BMI
    outfig5 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.BMICorrMegakaryocytes.",heytoday,".pdf")
    pdf(outfig5,height=4,width=4.8)
    megasubdt1 <- subset(sdt1long, Celltype=="Megakaryocytes")
    p51 <- ggplot(megasubdt1,aes(x=BMI, y=Proportion))+geom_point(color="black",size=2,shape=16)+commontheme()+labs(title="Year1",x="BMI",y="Megakaryocytes Proportion (%)")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(use_label(c("R", "P")),method="pearson",size=4.5,label.x="right",label.y="bottom")
    print(p51)
    megasubdt2 <- subset(sdt2long, Celltype=="Megakaryocytes")
    p52 <- ggplot(megasubdt2,aes(x=BMI, y=Proportion))+geom_point(color="black",size=2,shape=16)+commontheme()+labs(title="Year2",x="BMI",y="Megakaryocytes Proportion (%)")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(use_label(c("R", "P")),method="pearson",size=4.5,label.x="right",label.y="bottom")
    print(p52)
    dev.off()

    ###corr Age
    outfig6 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.AgeCorrMonocytes.",heytoday,".pdf")
    pdf(outfig6,height=4,width=4.8)
    monosubdt1 <- subset(sdt1long, Celltype=="Monocytes")
    p61 <- ggplot(monosubdt1,aes(x=Age, y=Proportion))+geom_point(color="black",size=2,shape=16)+commontheme()+labs(title="Year1",x="Age",y="Monocytes Proportion (%)")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(use_label(c("R", "P")),method="pearson",size=4.5,label.x="left",label.y="top")
    print(p61)
    monosubdt2 <- subset(sdt2long, Celltype=="Monocytes")
    p62 <- ggplot(monosubdt2,aes(x=Age, y=Proportion))+geom_point(color="black",size=2,shape=16)+commontheme()+labs(title="Year2",x="Age",y="Monocytes Proportion (%)")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(use_label(c("R", "P")),method="pearson",size=4.5,label.x="left",label.y="top")
    print(p62)
    dev.off()

    ###Celltype correlation
    corcelltype <- c("Granulocytes","Megakaryocytes","Erythrocyte-progenitors")
    cellcomb <- t(combn(corcelltype,2))
    outfig1 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.deconv.Celltypecor.",heytoday,".pdf")
    pdf(outfig1,height=5,width=5)
    for (i in 1:nrow(cellcomb)) {
        ct1 <- cellcomb[i,][1]
        ct2 <- cellcomb[i,][2]
        print(cellcomb[i,])
        subdt <- metainfo[,c("SubjectID","Samplename","Year","Time","Age","Sex",ct1,ct2)]
        subdt1 <- subset(subdt, Year=="Year1")
        p1 <- ggplot(subdt1,aes(x = subdt1[,ct1], y = subdt1[,ct2]))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year1",x=paste0(ct1," (%)"),y=paste0(ct2," (%)"))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_cor(aes(label=after_stat(paste0("italic(R)~","`=`~","'", r, "'","~italic(P)~","`=`~","'",p, "'"))),label.x.npc=0.4,size=4.5)
        print(p1)
        subdt2 <- subset(subdt, Year=="Year2")
        p2 <- ggplot(subdt2,aes(x = subdt2[,ct1], y = subdt2[,ct2]))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year2",x=paste0(ct1," (%)"),y=paste0(ct2," (%)"))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_cor(aes(label=after_stat(paste0("italic(R)~","`=`~","'", r, "'","~italic(P)~","`=`~","'",p, "'"))),label.x.npc=0.4,size=4.5)
        print(p2)
    }
    dev.off()

  
    celldt1 <- dt1[,colnames(dt1) %in% selectcelltype]
    colnames(celldt1) <- c("B","Gra","Hep","MK","EP","Mon","NK","T-CD4","VE")
    celldt2 <- dt2[,colnames(dt2) %in% selectcelltype]
    colnames(celldt2) <- c("B","Gra","Hep","MK","EP","Mon","NK","T-CD4","VE")
    #cellcomb <- t(combn(majorcells,2))

    my_fn <- function(data, mapping, ...){
      p <- ggplot(data = data, mapping = mapping) + geom_point(size=0.5,alpha=0.7) + geom_smooth(method=lm, fill="blue", color="blue", ...)
      p
    }

    outfig1 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.deconv.majorCelltypecor.pairwise.",heytoday,".pdf")
    pdf(outfig1,height=10.8,width=10.8)
    pp1 <- ggpairs(celldt1, columns = 1:9, lower=list(continuous=my_fn),title="Year1")+commontheme()+theme(axis.text.x=element_text(size=11,color="black"),axis.text.y=element_text(size=11,color="black"))
    print(pp1)
    pp2 <- ggpairs(celldt2, columns = 1:9, lower=list(continuous=my_fn),title="Year2")+commontheme()+theme(axis.text.x=element_text(size=11,color="black"),axis.text.y=element_text(size=11,color="black"))
    print(pp2)
    dev.off()
}

eval_inter_intra_deconv <- function(deconvresfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))
    
    ###Inter/intra variabilities
    #use the same age group for T1 and T2
    interdt <- metainfo
    interdt$Age <- ifelse(interdt$Year=="Year2",interdt$Age-1, interdt$Age)
    interstrata <- interdt %>% mutate(Age_group=cut(Age,breaks=c(19,29,39,49,59,69,100),labels=c("20s","30s","40s","50s","60s","70s+")))
    cordt <- interstrata[,c("Granulocytes","Monocytes","Megakaryocytes","Erythrocyte-progenitors","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes")]
    rownames(cordt) <- paste0(interstrata$SubjectID,"_",interstrata$Year,"_",interstrata$Age_group,"_",interstrata$Sex)
    tcordt <- t(cordt)
    cormat <- cor(tcordt,method="kendall")

    flattenCorrMatrix <- function(cormat) {
      ut <- upper.tri(cormat)
      data.frame(Sample1=rownames(cormat)[row(cormat)[ut]], Sample2=rownames(cormat)[col(cormat)[ut]], kcor=(cormat)[ut])
    }

    lcormat <- flattenCorrMatrix(cormat)
    rcormat <- separate(data=lcormat, col=Sample1, sep="_", into = c("SubjectID1", "Year1", "Agegroup1", "Sex1"))
    rcormat <- separate(data=rcormat, col=Sample2, sep="_", into = c("SubjectID2", "Year2", "Agegroup2", "Sex2"))

    l <- rep(list(0:1), 4)
    mycomb <- expand.grid(l)
    allrcormat <- data.frame(SubjectID1=character(),Year1=character(),Agegroup1=character(),Sex1=character(),SubjectID2=character(),Year2=character(),Agegroup2=character(),Sex2=character(),kcor=numeric(),Group=character())

    for (z in 1:nrow(mycomb)) {
        curcomb <- mycomb[z,]
        if (sum(curcomb) != 0) {
            curoperators <- ifelse(curcomb==0,'==','!=')
            #print(curoperators)
            eval1 <- curoperators[1]
            eval2 <- curoperators[2]
            eval3 <- curoperators[3]
            eval4 <- curoperators[4]
            subrcormat <- subset(rcormat, get(eval1)(rcormat$SubjectID1, rcormat$SubjectID2) & get(eval2)(rcormat$Year1, rcormat$Year2) & get(eval3)(rcormat$Agegroup1, rcormat$Agegroup2) & get(eval4)(rcormat$Sex1, rcormat$Sex2))
            if (dim(subrcormat)[1] != 0) {
                print(dim(subrcormat))
                #print(head(subrcormat))
                curnames <- ifelse(curcomb==0,'Same','Diff')
                groupname <- paste0(curnames[1], "_subject;", curnames[2], "_Year;", curnames[3], "_agegroup;", curnames[4], "_sex")
                print(groupname)
                print(summary(subrcormat$kcor))
                subrcormat$Group <- groupname
            }
            allrcormat <- rbind(allrcormat, subrcormat)
        }
    }

    allrcormat$Group1 <- ifelse(allrcormat$Group=="Same_subject;Diff_Year;Same_agegroup;Same_sex", "Intra individuals","Inter individuals")
    suballrcormat <- subset(allrcormat,Group=="Diff_subject;Diff_Year;Diff_agegroup;Same_sex" | Group=="Diff_subject;Diff_Year;Same_agegroup;Diff_sex")
    suballrcormat$Group2 <- ifelse(suballrcormat$Group=="Diff_subject;Diff_Year;Diff_agegroup;Same_sex","Inter individuals among same sex","Inter individuals among same age group")

    mu <- ddply(allrcormat, "Group1", summarise, grp.mean=round(mean(kcor),2))
    mu2 <- ddply(suballrcormat, "Group2", summarise, grp.mean=round(mean(kcor),2))
    outfig1 <- paste0(outdir,"/LONGVAR.deconv.epidish.REFsam.set8.kendallcorcoef.interintraVAR.",heytoday,".pdf")
    pdf(outfig1,height=4.5,width=5.8)
    p1 <- ggplot(allrcormat,aes(x=kcor,color=Group1,fill=Group1))+geom_histogram(aes(y=..density..), binwidth=0.06, alpha=0.6,position="identity")+geom_vline(data=mu, aes(xintercept=grp.mean, color=Group1),linetype="dashed")+geom_text(data=mu, aes(label=grp.mean, x=grp.mean, y=3.7), size=3.5, color="black")+commontheme()+theme(legend.title=element_blank(),legend.position=c(0.25,0.85))+labs(x="Correlation coefficient (Kendall)",y="Density")+scale_color_manual(values=c("#9970ab","#5aae61"))+scale_fill_manual(values=c("#9970ab","#5aae61"))
    print(p1)
    dev.off()
}
