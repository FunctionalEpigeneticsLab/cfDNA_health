#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort (methylation, size, motif)

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "ggpubr", "rstatix", "ggpmisc", "janitor", "circlize", "ComplexHeatmap")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))
heytoday <- Sys.Date()

args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
roimefh0 <- args[2]
roimefh1 <- args[3]
gt70motifh <- args[4]
lt30motifh <- args[5]
gt70sizestatfh <- args[6]
lt70sizestatfh <- args[7]
outdir <- args[8]

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

PlotRefROIme <- function(deconvresfh ,roime0, roime1, outdir) {
    metafh <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metafh$Year <- factor(metafh$Year, levels=c("T1","T2"),labels=c("Year1","Year2"))
    roime0 <- fread(roimefh0,header=TRUE,sep="\t",data.table=FALSE)
    roime1 <- fread(roimefh1,header=TRUE,sep="\t",data.table=FALSE)
    roime0dt <- roime0[,c("Chrom","Start","End","RefMeper","SamMeper","GCcode")]
    roime1dt <- roime1[,c("Chrom","Start","End","RefMeper","SamMeper","GCcode")]
    wroime0dt <- spread(roime0dt, GCcode, SamMeper)
    wroime1dt <- spread(roime1dt, GCcode, SamMeper)
    wroime0dt <- wroime0dt[,colnames(wroime0dt) %in% c("Chrom","Start","End","RefMeper",metafh$GCcode)]
    wroime1dt <- wroime1dt[,colnames(wroime1dt) %in% c("Chrom","Start","End","RefMeper",metafh$GCcode)]
    fwroime0dt <- wroime0dt %>% drop_na()
    fwroime1dt <- wroime1dt %>% drop_na()
    
    fwroime0mat <- as.matrix(fwroime0dt[,!(colnames(fwroime0dt) %in% c("Chrom","Start","End"))])
    fwroime1mat <- as.matrix(fwroime1dt[,!(colnames(fwroime1dt) %in% c("Chrom","Start","End"))])

    refCGsite0 <- fwroime0dt[,c("Chrom","Start","End")]
    refCGsite1 <- fwroime1dt[,c("Chrom","Start","End")]
    froime0 <- merge(refCGsite0, roime0, by=c("Chrom","Start","End"), all.x=TRUE)
    froime1 <- merge(refCGsite1, roime1, by=c("Chrom","Start","End"), all.x=TRUE)
    roime0dt <- froime0 %>% group_by(GCcode) %>% summarise(totun=sum(Uncount),totme=sum(Mecount)) %>% mutate(Acrossmc=totme/(totme+totun)) %>% select(GCcode, Acrossmc)
    roime1dt <- froime1 %>% group_by(GCcode) %>% summarise(totun=sum(Uncount),totme=sum(Mecount)) %>% mutate(Acrossmc=totme/(totme+totun)) %>% select(GCcode, Acrossmc)
    mdt0 <- merge(metafh, roime0dt, by=c("GCcode"), all.x=TRUE)
    mdt1 <- merge(metafh, roime1dt, by=c("GCcode"), all.x=TRUE)

    outfig0 <- paste0(outdir,"/LONGVAR.all.roime.acrossme.sam2ref.bar.",heytoday,".pdf")
    pdf(outfig0,height=6.8,width=11.2)
    p0 <- ggplot(mdt0, aes(x=SubjectID, y=Acrossmc))+geom_hline(yintercept=mean(wroime0dt$RefMeper),color="red")+geom_bar(stat="identity",width=0.3, fill="darkred")+commontheme()+scale_x_discrete(labels=NULL, breaks=NULL)+labs(x="Sample",y="Methylation level across highly methylated CpGs (2953)")+geom_text(x=20, y=mean(wroime0dt$RefMeper), color="black", label="Reference methylation level")+facet_grid(Year~.)
    print(p0)
    p1 <- ggplot(mdt1, aes(x=SubjectID, y=Acrossmc))+geom_hline(yintercept=mean(wroime1dt$RefMeper),color="red")+geom_bar(stat="identity",width=0.3, fill="darkblue")+commontheme()+scale_x_discrete(labels=NULL, breaks=NULL)+labs(x="Sample",y="Methylation level across lowly methylated CpGs (2752)")+geom_text(x=20, y=mean(wroime1dt$RefMeper), color="black", label="Reference methylation level")+facet_grid(Year~.)
    print(p1)
    dev.off()
}

AnaROIendmotif <- function(deconvresfh, gt70motifh, lt30motifh, outdir){
    metafh <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metafh$Year <- factor(metafh$Year, levels=c("T1","T2"),labels=c("Year1","Year2"))
    gt70motifdt <- fread(gt70motifh,header=TRUE,sep="\t",data.table=FALSE)
    lt30motifdt <- fread(lt30motifh,header=TRUE,sep="\t",data.table=FALSE)

    lgt70motif <- gather(gt70motifdt, endmotif, motifreq, AA:TT)
    lgt70motif$Mestate <- "Highly_methylated"
    llt30motif <- gather(lt30motifdt, endmotif, motifreq, AA:TT)
    llt30motif$Mestate <- "Lowly_methylated"

    smotifdt <- rbind(lgt70motif, llt30motif)
    mdt <- merge(metafh, smotifdt, by=c("GCcode"))
    mdt$Mestate <- factor(mdt$Mestate, levels=c("Highly_methylated","Lowly_methylated"))
    mdt$endmotif <- factor(mdt$endmotif)
    mdt1 <- subset(mdt, Year=="Year1")
    mdt2 <- subset(mdt, Year=="Year2")

    outfig0 <- paste0(outdir,"/LONGVAR.all.roime_gt70vslt30.endmotif.",heytoday,".pdf")
    pdf(outfig0,height=5.6,width=11.2)
    stat.test1 <- mdt1 %>% group_by(endmotif) %>% t_test(motifreq ~ Mestate, paired=TRUE) %>% adjust_pvalue (method="bonferroni") %>% add_xy_position(x="endmotif")
    p1 <- ggplot(mdt1, aes(x=endmotif, y=motifreq, color=Mestate))+geom_boxplot(aes(color=Mestate),position=position_dodge(0.5),width=0.7,outlier.shape=NA)+geom_dotplot(aes(fill=Mestate),binwidth=0.1,binaxis='y', stackdir='center', dotsize=0.6,position=position_dodge(0.5),alpha=0.5)+commontheme()+theme(legend.position="top")+stat_pvalue_manual(stat.test1,label="p.adj",size=3.5,vjust=0,tip.length=0)+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+labs(title="Year1",x="Motif",y="Motif frequency (%)")
    print(p1)
    stat.test2 <- mdt2 %>% group_by(endmotif) %>% t_test(motifreq ~ Mestate, paired=TRUE) %>% adjust_pvalue (method="bonferroni") %>% add_xy_position(x="endmotif")
    p2 <- ggplot(mdt2, aes(x=endmotif, y=motifreq, color=Mestate))+geom_boxplot(aes(color=Mestate),position=position_dodge(0.5),width=0.7,outlier.shape=NA)+geom_dotplot(aes(fill=Mestate),binwidth=0.1,binaxis='y', stackdir='center', dotsize=0.6,position=position_dodge(0.5),alpha=0.5)+commontheme()+theme(legend.position="top")+labs(title="Year2",x="Motif",y="Motif frequency (%)")+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+stat_pvalue_manual(stat.test2,label="p.adj",size=3.5,vjust=0,tip.length=0)
    print(p2)
    dev.off()

    selectedmotifs <- c("AT","CA","CC","CG","GC","GG","TA","TG")
    smdt1 <- mdt1[mdt1$endmotif %in% selectedmotifs,]
    smdt1$endmotif <- factor(smdt1$endmotif, levels=selectedmotifs)
    smdt2 <- mdt2[mdt2$endmotif %in% selectedmotifs,]
    smdt2$endmotif <- factor(smdt2$endmotif, levels=selectedmotifs)
    outfig1 <- paste0(outdir,"/LONGVAR.all.roime_gt70vslt30.selectedendmotif.paired.",heytoday,".pdf")
    pdf(outfig1,height=5,width=9.6)
    sub.stat.test1 <- stat.test1 %>% filter(endmotif %in% selectedmotifs) %>% add_xy_position(x="endmotif")
    p1 <- ggplot(smdt1, aes(x=Mestate, y=motifreq, group=SubjectID))+geom_line(size=0.2, color="gray",position=position_jitter(w=0.2, h=0, seed=1))+geom_point(aes(color=Mestate),size=0.6,position=position_jitter(w=0.2, h=0, seed=1))+commontheme()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="top")+stat_pvalue_manual(sub.stat.test1,label="p.adj",size=4,vjust=0,tip.length=0)+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+labs(x="",y="Motif frequency (%)")+facet_wrap(~endmotif,nrow=1)+scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+guides(color=guide_legend(override.aes = list(size=2))) 
    print(p1)
    sub.stat.test2 <- stat.test2 %>% filter(endmotif %in% selectedmotifs) %>% add_xy_position(x="endmotif")
    p2 <- ggplot(smdt2, aes(x=Mestate, y=motifreq, group=SubjectID))+geom_line(size=0.2, color="gray", position=position_jitter(w=0.2, h=0, seed=1))+geom_point(aes(color=Mestate),size=0.6,position=position_jitter(w=0.2, h=0, seed=1))+commontheme()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="top")+stat_pvalue_manual(sub.stat.test2,label="p.adj",size=4,vjust=0,tip.length=0)+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+labs(x="",y="Motif frequency (%)")+facet_wrap(~endmotif, nrow=1)+scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+guides(color=guide_legend(override.aes = list(size=2)))
    print(p2)
    dev.off()

    ## calculate mean differences of selected motifs
    hldt <- merge(lgt70motif, llt30motif, by=c("GCcode","endmotif"))
    subhldt0 <- hldt[hldt$endmotif %in% selectedmotifs,]
    subhldt <- merge(metafh[,c("GCcode","SubjectID","Year")], subhldt0, by=c("GCcode"))
    subhldt1 <- subset(subhldt, Year=="Year1")
    subhldt2 <- subset(subhldt, Year=="Year2")
    meandiff1 <- subhldt1 %>% group_by(endmotif) %>% summarise(mean_diff=mean(motifreq.x-motifreq.y, na.rm=TRUE))
    meandiff2 <- subhldt2 %>% group_by(endmotif) %>% summarise(mean_diff=mean(motifreq.x-motifreq.y, na.rm=TRUE))
}


AnaROIsizestat <- function(deconvresfh, gt70sizestatfh, lt30sizestatfh, outdir) {
    metafh <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metafh$Year <- factor(metafh$Year, levels=c("T1","T2"),labels=c("Year1","Year2"))
    gt70sizestat <- fread(gt70sizestatfh,header=TRUE,sep="\t",data.table=FALSE)
    lt30sizestat <- fread(lt30sizestatfh,header=TRUE,sep="\t",data.table=FALSE)
    gt70sizestat <- gt70sizestat[,c("GCcode","Sizearea_sub","Sizearea_mo","Sizearea_inter","Sizearea_di")]
    lt30sizestat <- lt30sizestat[,c("GCcode","Sizearea_sub","Sizearea_mo","Sizearea_inter","Sizearea_di")]

    lgt70sizestat <- gather(gt70sizestat, Sizearea, cumfreq, Sizearea_sub:Sizearea_di)
    lgt70sizestat$Mestate <- "Highly_methylated"
    llt30sizestat <- gather(lt30sizestat, Sizearea, cumfreq, Sizearea_sub:Sizearea_di)
    llt30sizestat$Mestate <- "Lowly_methylated"

    lsizestat <- rbind(lgt70sizestat, llt30sizestat)
    mdt <- merge(metafh, lsizestat, by=c("GCcode"))
    mdt$Mestate <- factor(mdt$Mestate, levels=c("Highly_methylated","Lowly_methylated"))
    mdt$Sizearea <- factor(mdt$Sizearea, levels=c("Sizearea_sub","Sizearea_mo","Sizearea_inter","Sizearea_di"), labels=c("50-145bp","167-187bp","200-329bp","330-420bp"))
    mdt1 <- subset(mdt, Year=="Year1")
    mdt2 <- subset(mdt, Year=="Year2")
    mdt1$cumfreq <- mdt1$cumfreq*100
    mdt2$cumfreq <- mdt2$cumfreq*100

    outfig0 <- paste0(outdir,"/LONGVAR.all.roime_gt70vslt30.sizestats.",heytoday,".pdf")
    pdf(outfig0,height=5,width=7.2)
    stat.test1 <- mdt1 %>% group_by(Sizearea) %>% t_test(cumfreq ~ Mestate, paired=TRUE) %>% adjust_pvalue(method="bonferroni") %>% add_xy_position(x="Mestate")
    p1 <- ggplot(mdt1, aes(x=Mestate, y=cumfreq, group=SubjectID))+geom_line(size=0.3, color="gray", position=position_jitter(seed=9,width=0.1, height=0))+geom_point(aes(color=Mestate),size=0.7,position=position_jitter(seed=9, width=0.1, height=0))+commontheme()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="top")+stat_pvalue_manual(stat.test1,label="p.adj",size=4,vjust=0,tip.length=0)+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+labs(x="Size range",y="Cumulative size frequencies (%)")+facet_wrap(~Sizearea, nrow=1)+scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+guides(color=guide_legend(override.aes = list(size=2)))
    print(p1)
    stat.test2 <- mdt2 %>% group_by(Sizearea) %>% t_test(cumfreq ~ Mestate, paired=TRUE) %>% adjust_pvalue(method="bonferroni") %>% add_xy_position(x="Mestate")
    p2 <- ggplot(mdt2, aes(x=Mestate, y=cumfreq, group=SubjectID))+geom_line(size=0.3, color="gray", position=position_jitter(seed=9,width=0.1, height=0))+geom_point(aes(color=Mestate),size=0.7,position=position_jitter(seed=9,width=0.1, height=0))+commontheme()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="top")+stat_pvalue_manual(stat.test2,label="p.adj",size=4,vjust=0,tip.length=0)+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))+labs(x="Size range",y="Cumulative size frequencies (%)")+facet_wrap(~Sizearea, nrow=1)+scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+guides(color=guide_legend(override.aes = list(size=2)))
    print(p2)
    dev.off()

    ## calculate mean differences of size statistics
    hldt <- merge(lgt70sizestat, llt30sizestat, by=c("GCcode","Sizearea"))
    subhldt0 <- merge(metafh[,c("GCcode","SubjectID","Year")], hldt, by=c("GCcode"))
    subhldt1 <- subset(subhldt0, Year=="Year1")
    subhldt2 <- subset(subhldt0, Year=="Year2")
    meandiff1 <- subhldt1 %>% group_by(Sizearea) %>% summarise(mean_diff=mean(100*cumfreq.x-100*cumfreq.y, na.rm=TRUE))
    meandiff2 <- subhldt2 %>% group_by(Sizearea) %>% summarise(mean_diff=mean(100*cumfreq.x-100*cumfreq.y, na.rm=TRUE))
}

ReadinPooledSize <- function(fh, inYear, roiname) {
    dt <- fread(fh, header=FALSE, sep="\t")
    dt <- as.data.frame(table(dt[,V1]))
    colnames(dt) <- c("Size","Count")
    dt$Freq <- 100*dt$Count/sum(dt$Count)
    dt$Year <- inYear
    dt$MEROI <- roiname
    return(dt)
}

t1gt70fh <- "pooled.T1.meroi.gt70.size.txt.gz"
t1lt30fh <- "pooled.T1.meroi.lt30.size.txt.gz"
t2gt70fh <- "pooled.T2.meroi.gt70.size.txt.gz"
t2lt30fh <- "pooled.T2.meroi.lt30.size.txt.gz"

PlotPooledSize <- function(outdir) {
    t1h <- ReadinPooledSize(t1gt70fh, "T1", "Highly_methylated")
    t1l <- ReadinPooledSize(t1lt30fh, "T1", "Lowly_methylated")
    t2h <- ReadinPooledSize(t2gt70fh, "T2", "Highly_methylated")
    t2l <- ReadinPooledSize(t2lt30fh, "T2", "Lowly_methylated")

    t1 <- rbind(t1h, t1l)
    t2 <- rbind(t2h, t2l)

    t1$MEROI <- factor(t1$MEROI, levels=c("Highly_methylated","Lowly_methylated"))
    t2$MEROI <- factor(t2$MEROI, levels=c("Highly_methylated","Lowly_methylated"))
    t1$Size <- as.numeric(as.character(t1$Size))
    t2$Size <- as.numeric(as.character(t2$Size))

    outfig0 <- paste0(outdir,"/LONGVAR.all.roime_gt70vslt30.sizeprofile.",heytoday,".pdf")
    pdf(outfig0,height=4,width=5.8)
    p1 <- ggplot(t1,aes(x=Size,y=Freq,color=MEROI))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=0.8)+commontheme()+theme(legend.title=element_blank(),legend.position=c(0.75,0.8))+labs(x="Size (bp)",y="Frequency (%)")+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")+scale_color_manual(values=c("darkred","darkblue"))
    print(p1)
    p2 <- ggplot(t2,aes(x=Size,y=Freq,color=MEROI))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=0.8)+commontheme()+theme(legend.title=element_blank(),legend.position=c(0.75,0.8))+labs(x="Size (bp)",y="Frequency (%)")+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")+scale_color_manual(values=c("darkred","darkblue"))
    print(p2)
    dev.off()
    
    ###size difference
    t1diff <- merge(t1l, t1h, by=c("Size","Year"))
    t2diff <- merge(t2l, t2h, by=c("Size","Year"))
    t1diff$Size <- as.numeric(as.character(t1diff$Size))
    t2diff$Size <- as.numeric(as.character(t2diff$Size))
    t1diff$Freqdiff <- t1diff$Freq.x - t1diff$Freq.y
    t2diff$Freqdiff <- t2diff$Freq.x - t2diff$Freq.y

    outfig0 <- paste0(outdir,"/LONGVAR.all.roime_gt70vslt30.sizeprofile.freqdiff.",heytoday,".pdf")
    pdf(outfig0,height=4,width=5.8)
    p1 <- ggplot(t1diff,aes(x=Size,y=Freqdiff))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_hline(yintercept=0,linetype="dashed",color = "gray", linewidth=1)+geom_line(color="black",linewidth=1)+commontheme()+theme(legend.title=element_blank(),legend.position=c(0.75,0.8))+labs(x="Size (bp)",y="Frequency difference (%)\nlowly to highly methylation")+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")+scale_color_manual(values=c("darkred","darkblue"))+scale_y_continuous(expand=expansion(mult = c(0.1, 0.1)))
    print(p1)
    p2 <- ggplot(t2diff,aes(x=Size,y=Freqdiff))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_hline(yintercept=0,linetype="dashed",color = "gray", linewidth=1)+geom_line(color="black",linewidth=1)+commontheme()+theme(legend.title=element_blank(),legend.position=c(0.75,0.8))+labs(x="Size (bp)",y="Frequency difference (%)\nlowly to highly methylation")+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")+scale_color_manual(values=c("darkred","darkblue"))+scale_y_continuous(expand=expansion(mult = c(0.1, 0.1)))
    print(p2)
    dev.off()
}

cutfh0 <- "LONGVAR.cleavageprofile.allsample.Majorcell_gt70.tsv"
cutfh1 <- "LONGVAR.cleavageprofile.allsample.Majorcell_lt30.tsv"

PlotCleavageCorrConc <- function(deconvresfh, outdir) {
    cut0 <- fread(cutfh0,header=TRUE,sep="\t",data.table=FALSE)
    cut1 <- fread(cutfh1,header=TRUE,sep="\t",data.table=FALSE)

    metafh <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)

    mdt0 <- merge(metafh, cut0, by=c("GCcode"))
    mdt01 <- subset(mdt0, Year=="T1")
    mdt02 <- subset(mdt0, Year=="T2")

    mdt1 <- merge(metafh, cut1, by=c("GCcode"))
    mdt11 <- subset(mdt1, Year=="T1")
    mdt12 <- subset(mdt1, Year=="T2")

    outfig1 <- paste0(outdir,"/LONGVAR.all.roime_gt70.pos6cut.corr.logGE.",heytoday,".pdf")
    pdf(outfig1,height=5,width=5)
    p0 <- ggplot(mdt01, aes(x=`6`, y=Adj_Unique_Molecules))+geom_point(color="black",size=2,alpha=0.7,shape=16)+geom_smooth(method="lm",color="darkblue", se=FALSE, linetype="dashed")+commontheme()+labs(title="Year1",x="Cleavage % at position C",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+stat_cor(method="spearman",label.y=0.92,size=5)+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))
    print(p0)
    p0 <- ggplot(mdt02, aes(x=`6`, y=Adj_Unique_Molecules))+geom_point(color="black",size=2,alpha=0.7,shape=16)+geom_smooth(method="lm",color="darkblue", se=FALSE, linetype="dashed")+commontheme()+labs(title="Year2",x="Cleavage % at position C",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+stat_cor(method="spearman",label.y=0.92,size=5)+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))
    print(p0)
    dev.off()

    outfig1 <- paste0(outdir,"/LONGVAR.all.roime_lt30.pos6cut.corr.logGE.",heytoday,".pdf")
    pdf(outfig1,height=5,width=5)
    p0 <- ggplot(mdt11, aes(x=`6`, y=Adj_Unique_Molecules))+geom_point(color="black",size=2,alpha=0.7,shape=16)+geom_smooth(method="lm",color="darkblue", se=FALSE, linetype="dashed")+commontheme()+labs(title="Year1",x="Cleavage % at position C",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+stat_cor(method="spearman",label.y=0.92,size=5)+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))
    print(p0)
    p0 <- ggplot(mdt12, aes(x=`6`, y=Adj_Unique_Molecules))+geom_point(color="black",size=2,alpha=0.7,shape=16)+geom_smooth(method="lm",color="darkblue", se=FALSE, linetype="dashed")+commontheme()+labs(title="Year2",x="Cleavage % at position C",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+stat_cor(method="spearman",label.y=0.92,size=5)+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))
    print(p0)
    dev.off()
}
