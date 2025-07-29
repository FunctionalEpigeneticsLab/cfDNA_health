#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "reshape2", "ggpubr", "rstatix", "ggpmisc", "scales")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
qubitfh <- args[2]
outdir <- args[3]


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

eval_qubit_spikein_conc <- function(deconvresfh, qubitfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    qubit <- fread(qubitfh, header=TRUE, sep="\t", data.table=FALSE)
    colnames(qubit) <- c("Libraryname","Samplename","GCcode","Qubit_concentration")
    mdt <- merge(metainfo, qubit, by=c("GCcode"))

    outfig1 <- paste0(outdir,"/LONGVAR.estimated_adj_molecules_cor_qubit.",heytoday,".pdf")
    pdf(outfig1,height=4.5,width=4.5)
    p1 <- ggplot(mdt,aes(x=Adj_Unique_Molecules, y=Qubit_concentration))+geom_point(fill="black",size=2,shape=21)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=14),axis.title=element_text(size=16,color="black"),legend.position="none",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Estimated cfDNA concentration\n(Spike-in to sample ratio)",y="Qubit (ng/µl)")+geom_smooth(method=lm,se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="left",use_label(c("R", "P")))
    print(p1)
    dev.off()
}

eval_cross_conc <- function(deconvresfh, outdir) {
    metainfo <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metainfo$Year <- factor(metainfo$Year, levels=c("T1", "T2"), labels=c("Year1", "Year2"))
    metainfo$Time <- factor(metainfo$Time, levels=c("Morning","Afternoon"))
    metainfo$Sex <- factor(metainfo$Sex, levels=c("Female", "Male"))
    metainfo$log10_Adj_Unique_Molecules <- log10(metainfo$Adj_Unique_Molecules)

    GEdt0 <- metainfo[,c("SubjectID","Year","Sex","Adj_Unique_Molecules")]
    wideGEdt0 <- spread(GEdt0, Year, Adj_Unique_Molecules)

    ##
    outfig0 <- paste0(outdir, "/LONGVAR.estimated_adj_molecules.Y1vsY2.",heytoday,".pdf")
    pdf(outfig0, height=5, width=5)
    p0 <- ggplot(wideGEdt0, aes(x=Year1,y=Year2))+geom_point(shape=16,size=2,color="black")+geom_abline(intercept=0,linetype="dashed",color="gray")+commontheme()+labs(x="Estimated cfDNA concentration at Year1", y="Estimated cfDNA concentration at Year2")+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_cor(label.x=50, label.y=950, size=4.5) +xlim(15, 1015)+ylim(15, 1015)
    print(p0)
    dev.off()

    outfig1 <- paste0(outdir,"/LONGVAR.estimated_adj_molecules.density.",heytoday,".pdf")
    pdf(outfig1,height=4,width=5.2)
    p1 <- ggplot(metainfo, aes(x=Adj_Unique_Molecules))+geom_histogram(aes(y=..density.. ,color=Year,fill=Year), alpha=0.5, position="identity")+geom_density(aes(color=Year),linewidth=1)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position=c(0.85,0.8),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Estimated cfDNA concentration\n(Spike-in to sample ratio)",y="Density")+scale_color_manual(values=c("gray","royalblue4"))+scale_fill_manual(values=c("gray","royalblue4"))
    print(p1)
    p1 <- ggplot(metainfo, aes(x=Adj_Unique_Molecules))+geom_histogram(aes(y=..density.. ,color=Year,fill=Year), alpha=0.5, position="identity")+geom_density(aes(color=Year),linewidth=1)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=16,color="black"),legend.position=c(0.85,0.8),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)",y="Density")+scale_color_manual(values=c("gray","royalblue4"))+scale_fill_manual(values=c("gray","royalblue4"))+scale_x_log10()
    print(p1)
    dev.off()

    ###Sex difference
    T1metainfo <- metainfo[metainfo$Year=="Year1",]
    T2metainfo <- metainfo[metainfo$Year=="Year2",]
    outfig2 <- paste0(outdir, "/LONGVAR.estimated_adj_molecules.Sex.",heytoday,".pdf")
    pdf(outfig2,height=4.8,width=4)
    stat.test <- T1metainfo %>% group_by(Year) %>% wilcox_test(log10_Adj_Unique_Molecules ~ Sex) %>% add_xy_position(x = "Sex")
    p2 <- ggplot(T1metainfo, aes(x=Sex,y=Adj_Unique_Molecules))+geom_boxplot(aes(color=Sex),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Sex),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),legend.position="none")+labs(title="Year1",x="",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_color_manual(values=c("#EA7317","#2364AA"))+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000),expand=c(0.05, 0.1))+stat_pvalue_manual(stat.test,label="p",size=4.5,step.increase=0.02,vjust=0,tip.length=0)
    print(p2)
    stat.test <- T2metainfo %>% group_by(Year) %>% wilcox_test(log10_Adj_Unique_Molecules ~ Sex) %>% add_xy_position(x = "Sex")
    p2 <- ggplot(T2metainfo, aes(x=Sex,y=Adj_Unique_Molecules))+geom_boxplot(aes(color=Sex),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Sex),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),legend.position="none")+labs(title="Year2",x="",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_color_manual(values=c("#EA7317","#2364AA"))+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000),expand=c(0.05, 0.1))+stat_pvalue_manual(stat.test,label="p",size=4.5,step.increase=0.02,vjust=0,tip.length=0)
    print(p2)
    dev.off()


    ###Corr BMI
    outfig2 <- paste0(outdir, "/LONGVAR.estimated_adj_molecules.corBMI.", heytoday, ".pdf")
    pdf(outfig2,height=5,width=5)
    p2 <- ggplot(T1metainfo, aes(x=BMI,y=Adj_Unique_Molecules))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year1",x="BMI",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="right",label.y="bottom",use_label(c("R", "P")),method="pearson",size=4.5)
    print(p2)
    p2 <- ggplot(T2metainfo, aes(x=BMI,y=Adj_Unique_Molecules))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year2",x="BMI",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_correlation(label.x="right",label.y="bottom",use_label(c("R", "P")),method="pearson",size=4.5)
    print(p2)
    dev.off()

    ###Corr Age
    outfig2 <- paste0(outdir, "/LONGVAR.estimated_adj_molecules.corAge.", heytoday, ".pdf")
    pdf(outfig2,height=5,width=5)
    p2 <- ggplot(T1metainfo, aes(x=Age,y=Adj_Unique_Molecules))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year1",x="Age",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_cor(aes(label=after_stat(paste0("italic(R)~","`=`~","'", r, "'","~italic(P)~","`=`~","'",p, "'"))),label.x.npc = "center",label.y.npc = "bottom",size=4.5)
    print(p2)
    p2 <- ggplot(T2metainfo, aes(x=Age,y=Adj_Unique_Molecules))+geom_point(color="black",size=2,shape=16)+commontheme()+theme(legend.position="none")+labs(title="Year2",x="Age",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_y_log10(breaks=c(20,50,100,200,400,700,1000),labels=c(20,50,100,200,400,700,1000))+geom_smooth(method=lm, se=FALSE,color="darkblue",linetype="dashed")+stat_cor(aes(label=after_stat(paste0("italic(R)~","`=`~","'", r, "'","~italic(P)~","`=`~","'",p, "'"))),label.x.npc = "center",label.y.npc = "bottom",size=4.5)
    print(p2)
    dev.off()
}
