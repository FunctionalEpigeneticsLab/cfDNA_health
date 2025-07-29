#!/usr/bin/env Rscript
## This script is used to generate plots for the Diurnal cohort

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "reshape2", "ggpubr", "rstatix", "ggpmisc", "scales", "cosinor")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
sizecountfh <- args[2]
sizestatsfh <- args[3]
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

eval_diurnal_concentration <- function(deconvresfh, outdir) {
    dt <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    dt$Day <- factor(dt$Day, levels=c("Day1", "Day2", "Day3"))
    dt$Time <- factor(dt$Time, levels=c("Morning","Afternoon","Evening"))
    dt$Sex <- factor(dt$Sex, levels=c("F", "M"), labels=c("Female","Male"))

    dt %>% select(Samplename,SubjectID,Day,Time,Adj_Unique_Molecules) %>% group_by(Day,Time) %>% summarize(meanconc=mean(Adj_Unique_Molecules))

    GEdt <- dt[,c("Samplename","GCcode","SubjectID","Day","Time","Age","Sex","Adj_Unique_Molecules")]
    GEdt$logGE <- log10(GEdt$Adj_Unique_Molecules)
    GEdt %>% select(Samplename,SubjectID,Day,Time,logGE) %>% group_by(Day,Time) %>% summarize(logmeanconc=mean(logGE))
    
    ###statistical test
    sGEdt <- GEdt[,c("SubjectID","Time","Day","logGE")]
    sGEdt %>% group_by(Time, Day) %>% shapiro_test(logGE)
    sGEdt %>% group_by(Time) %>% levene_test(logGE ~ Day)

    ######Two-way mixed ANOVA test
    res.aov <- anova_test(data = sGEdt, dv = logGE, wid = SubjectID, within = c(Day,Time))
    get_anova_table(res.aov)

    ######Effect of time at each level of day
    one.way2 <- sGEdt %>% group_by(Day) %>% anova_test(dv = logGE, wid = SubjectID, within = Time) %>% get_anova_table() %>% adjust_pvalue(method = "bonferroni")
    print(one.way2)

    pwc2 <- sGEdt %>% group_by(Day) %>% pairwise_t_test(logGE ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% select(-df, -statistic, -p) # Remove details
    pwc2

    ###plot boxplot
    outfig1 <- paste0(outdir, "/Diurnal.log10_estimated_adj_molecules.grouped_time.", heytoday, ".pdf")
    pdf(outfig1,height=5.6,width=4.5)
    p1 <- ggplot(GEdt, aes(x=Time,y=Adj_Unique_Molecules))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+labs(x="",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+annotate(geom="text", x=2, y=365, label=paste0("Repeated~measures~ANOVA:italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)+scale_y_log10(breaks=c(20,50,80,120,150,210,270,330),labels=c(20,50,80,120,150,210,270,330),expand = c(0.01, 0.08))
    print(p1)
    dev.off()

    outfig2 <- paste0(outdir, "/Diurnal.log10_estimated_adj_molecules.stratify_day.", heytoday, ".pdf")
    pdf(outfig2,height=4.2,width=7.8)
    stat.test1 <- GEdt %>% group_by(Day) %>% pairwise_t_test(logGE ~ Time, paired=TRUE, p.adjust.method="bonferroni") %>% add_xy_position(x = "Time")
    p1 <- ggplot(GEdt, aes(x=Time,y=Adj_Unique_Molecules))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+labs(x="",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+scale_y_log10(breaks=c(20,50,80,120,150,210,270,330),labels=c(20,50,80,120,150,210,270,330),expand = c(0.01, 0.08))+stat_pvalue_manual(stat.test1,label="p.adj",size=3.5,step.increase=0.02,vjust=0,tip.length=0)+facet_grid(.~Day)
    print(p1)
    dev.off()

    ###plot individual
    GEdt$Series <- paste0(GEdt$Day,"_",GEdt$Time)
    GEdt$Series <- factor(GEdt$Series, levels=c("Day1_Morning","Day1_Afternoon","Day1_Evening","Day2_Morning","Day2_Afternoon","Day2_Evening","Day3_Morning","Day3_Afternoon","Day3_Evening"))
    GEdt$SubjectID <- factor(GEdt$SubjectID, levels=c("1DIU","2DIU","3DIU","4DIU","5DIU","6DIU","7DIU","8DIU","9DIU","10DIU","11DIU","12DIU","13DIU","14DIU","15DIU","16DIU"))

    outfig3 <- paste0(outdir, "/Diurnal.log10_estimated_adj_molecules.individual.", heytoday, ".pdf")
    pdf(outfig3,height=8.2,width=10.2)
    p2 <- ggplot(GEdt,aes(x=Series,y=Adj_Unique_Molecules,group=1))+geom_line(color="black",linewidth=0.8, position=position_dodge(0.3))+geom_point(position=position_dodge(0.3), size=2, shape=21, fill="white")+commontheme()+theme(axis.text.y=element_text(size=10,color="black"),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=11,color="black"),legend.text=element_text(size=12),axis.title=element_text(size=14,color="black"),legend.position="right",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="",y="Estimated cfDNA concentration\n(Spike-in to sample ratio, log10)")+facet_wrap(~SubjectID)+scale_y_log10(breaks=c(20,50,80,150,210,330),labels=c(20,50,80,150,210,330),expand = c(0.01, 0.08))
    print(p2)
    dev.off()
}

eval_diurnal_poolsize <- function(sizecountfh, outdir) {
    fh <- fread(sizecountfh,header=TRUE,sep="\t",data.table=FALSE)
    timedt <- fh %>% group_by(Time,Size) %>% summarise(totcount = sum(Count)) %>% group_by(Time) %>% mutate(Freq=100*totcount/sum(totcount)) %>% select(Time, Size, Freq)
    timedt$Time <- factor(timedt$Time,levels=c("Morning","Afternoon","Evening"))
    timedt <- timedt[timedt$Size>=50 & timedt$Size <=450,]
    longtimedt <- spread(timedt, Time, Freq)
    longtimedt$M2A <- longtimedt$Morning-longtimedt$Afternoon
    longtimedt$M2E <- longtimedt$Morning-longtimedt$Evening
    longtimedt1 <- longtimedt[,c("Size","M2A","M2E")]
    timedt2 <- gather(longtimedt1,compgroup,freqdiff,M2A:M2E)
    timedt2$compgroup <- factor(timedt2$compgroup,levels=c("M2A","M2E"),labels=c("Morning-Afternoon","Morning-Evening"))

    outfig <- paste0(outdir, "/Diurnal.pooled.sizeprofile.time.",heytoday,".pdf")
    pdf(outfig,width=5.8,height=4)
    p0 <- ggplot(timedt,aes(x=Size,y=Freq,color=Time))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=0.8,alpha=0.7)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=14),axis.title=element_text(size=16,color="black"),legend.title=element_blank(),legend.background=element_blank(),legend.position=c(0.85,0.85),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency (%)")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
    print(p0)
    p1 <- ggplot(timedt2,aes(x=Size,y=freqdiff,color=compgroup))+geom_hline(yintercept=0,color="#EFC000",linewidth=1)+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=1)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),legend.text=element_text(size=14),axis.title=element_text(size=16,color="black"),legend.title=element_blank(),legend.background=element_blank(),legend.position=c(0.75,0.8),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency difference (%)")+scale_color_manual(values=c("darkred","royalblue4"))+scale_x_continuous(limits=c(50, 450),breaks=seq(50,450,50))+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
    print(p1)
    dev.off()
}

eval_diurnal_sizestats <- function(deconvresfh, sizestatsfh, outdir) {
    metainfo <- fread(deconvresfh, header=TRUE, sep="\t", data.table=FALSE)
    sizestats <- fread(sizestatsfh, header=TRUE, sep="\t", data.table=FALSE)
    dt0 <- merge(metainfo, sizestats, by=c("GCcode"))
    dt0$LogAdjMol <- log10(dt0$Adj_Unique_Molecules)
    dt0$Subject <- as.factor(dt0$Subject)
    dt0$Time <- factor(dt0$Time,levels=c("Morning","Afternoon","Evening"))
    dt0$Day <- factor(dt0$Day,levels=c("Day1","Day2","Day3"))
    dt0$Sizearea_mo <- dt0$Sizearea_mo*100

    sFSdt <- dt0[,c("SubjectID","Time","Day","Sizearea_mo")]
    sFSdt %>% group_by(Time, Day) %>% shapiro_test(Sizearea_mo)
    sFSdt %>% group_by(Time) %>% levene_test(Sizearea_mo ~ Day)
    res.aov <- anova_test(data = sFSdt, dv = Sizearea_mo, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    one.way2 <- sFSdt %>% group_by(Day) %>% anova_test(dv = Sizearea_mo, wid = SubjectID, within = Time) %>% get_anova_table() %>% adjust_pvalue(method = "bonferroni")
    one.way2
    pwc2 <- sFSdt %>% group_by(Day) %>% pairwise_t_test(Sizearea_mo ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% select(-df, -statistic, -p) # Remove details
    pwc2

    outfig0 <- paste0(outdir,"/Diurnal.Sizearea_monoNuc.time.",heytoday,".pdf")
    pdf(outfig0, height=5.6,width=4.5)
    p0 <- ggplot(sFSdt, aes(x=Time,y=Sizearea_mo))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(x="",y="Cumulative size frequencies (%)\n(167-187bp)")+annotate(geom="text", x=2, y=max(sFSdt$Sizearea_mo)*1.1, label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p0)
    dev.off()

    ###
    dt0$Sizearea_sub <- dt0$Sizearea_sub*100

    sFSdt <- dt0[,c("SubjectID","Time","Day","Sizearea_sub")]
    sFSdt %>% group_by(Time, Day) %>% shapiro_test(Sizearea_sub)
    sFSdt %>% group_by(Time) %>% levene_test(Sizearea_sub ~ Day)
    res.aov <- anova_test(data = sFSdt, dv = Sizearea_sub, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    one.way2 <- sFSdt %>% group_by(Day) %>% anova_test(dv = Sizearea_sub, wid = SubjectID, within = Time) %>% get_anova_table() %>% adjust_pvalue(method = "bonferroni")
    one.way2
    pwc2 <- sFSdt %>% group_by(Day) %>% pairwise_t_test(Sizearea_sub ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% select(-df, -statistic, -p) # Remove details
    pwc2

    outfig0 <- paste0(outdir,"/Diurnal.Sizearea_subNuc.time.",heytoday,".pdf")
    pdf(outfig0, height=5.6,width=4.5)
    p0 <- ggplot(sFSdt, aes(x=Time,y=Sizearea_sub))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(x="",y="Cumulative size frequencies (%)\n(50-145bp)")+annotate(geom="text", x=2, y=max(sFSdt$Sizearea_sub)*1.1, label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p0)
    dev.off()

    dt0$Sizearea_di <- dt0$Sizearea_di*100

    sFSdt <- dt0[,c("SubjectID","Time","Day","Sizearea_di")]
    sFSdt %>% group_by(Time, Day) %>% shapiro_test(Sizearea_di)
    sFSdt %>% group_by(Time) %>% levene_test(Sizearea_di ~ Day)
    res.aov <- anova_test(data = sFSdt, dv = Sizearea_di, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    one.way2 <- sFSdt %>% group_by(Day) %>% anova_test(dv = Sizearea_di, wid = SubjectID, within = Time) %>% get_anova_table() %>% adjust_pvalue(method = "bonferroni")
    one.way2
    pwc2 <- sFSdt %>% group_by(Day) %>% pairwise_t_test(Sizearea_di ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% select(-df, -statistic, -p) # Remove details
    pwc2

    outfig0 <- paste0(outdir,"/Diurnal.Sizearea_diNuc.time.",heytoday,".pdf")
    pdf(outfig0, height=5.6,width=4.5)
    p0 <- ggplot(sFSdt, aes(x=Time,y=Sizearea_di))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(x="",y="Cumulative size frequencies (%)\n(330-420bp)")+annotate(geom="text", x=2, y=max(sFSdt$Sizearea_di)*1.1, label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p0)
    dev.off()

    ###stratified by day
    sFSdt <- dt0[,c("SubjectID","Time","Day","Sizearea_mo")]
    outfig1 <- paste0(outdir, "/Diurnal.Sizearea_monoNuc.stratify_day.", heytoday, ".pdf")
    pdf(outfig1,height=4.2,width=7.8)
    stat.test1 <- sFSdt %>% group_by(Day) %>% pairwise_t_test(Sizearea_mo ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Time")
    p1 <- ggplot(sFSdt, aes(x=Time,y=Sizearea_mo))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(x="",y="Cumulative size frequencies (%)\n(167-187bp)")+stat_pvalue_manual(stat.test1,label="p.adj",size=3.5,step.increase=0.02,vjust=0,tip.length=0)+facet_grid(.~Day)
    print(p1)
    dev.off()

    sFSdt <- dt0[,c("SubjectID","Time","Day","Sizearea_sub")]
    outfig1 <- paste0(outdir, "/Diurnal.Sizearea_subNuc.stratify_day.", heytoday, ".pdf")
    pdf(outfig1,height=4.2,width=7.8)
    stat.test1 <- sFSdt %>% group_by(Day) %>% pairwise_t_test(Sizearea_sub ~ Time, paired = TRUE, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Time")
    p1 <- ggplot(sFSdt, aes(x=Time,y=Sizearea_sub))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(x="",y="Cumulative size frequencies (%)\n(50-145bp)")+stat_pvalue_manual(stat.test1,label="p.adj",size=3.5,step.increase=0.02,vjust=0,tip.length=0)+facet_grid(.~Day)
    print(p1)
    dev.off()
}

eval_diurnal_deconv <- function(deconvresfh, outdir) {
    dt <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    dt$Day <- factor(dt$Day, levels=c("Day1", "Day2", "Day3"))
    dt$Time <- factor(dt$Time, levels=c("Morning","Afternoon","Evening"))
    dt$Sex <- factor(dt$Sex, levels=c("F", "M"), labels=c("Female","Male"))

    col2drop <- c("Ethnicity","GCcode","Total_Reads","Unique_Reads","Duplication_Rate","CpGs_Methylated_Percentage","CHG_Methylated_Percentage","CHH_Methylated_Percentage","chrX_reads","chrY_reads","chrM_reads","Lambda_reads","pUC_reads","Lambda_Methylated_Percentage","pUC_Methylated_Percentage","Percentage_Aligned_on_Target","Raw_Target_reads","Raw_Autosome_reads","Raw_Lambda_reads","Duplication_Rate_Target","Duplication_Rate_Autosome","Duplication_Rate_Lambda","Abs_pUC_Input","Abs_Lambda_Input")
    dt0 <- dt[,!(colnames(dt) %in% col2drop)]

    dt0long <- melt(dt0, id.vars=c("SubjectID","Samplename","Day","Time","SamplingTime","Age","Sex","Height","Weight","BMI","Mean_Coverage","Mean_Insert_size","Median_Insert_size","Mean_Ontar_Insert_size","Median_Ontar_Insert_size","Mode_Ontar_Insert_size","Sizeratio1","Sizeratio2","ShortFragProp","Raw_Molecules","Unique_Molecules","Adj_Unique_Molecules","Covfailed"),value.name="Proportion",variable.name="Celltype")

    ## major cell types diurnal changes
    majorcelltype <- c("Granulocytes","Monocytes","Megakaryocytes","Erythrocyte-progenitors","NK-cells","B-cells","T-cells-CD4","Vascular-endothelium","Liver-hepatocytes")
    plotdt0 <- dt0long[dt0long$Celltype %in% majorcelltype,]
    plotdt <- plotdt0 %>% group_by(Day, Time, Sex, Celltype) %>% summarise(meanprop=mean(Proportion),sdprop=sd(Proportion))
    plotdt$ts <- paste0(plotdt$Day, "_", plotdt$Time)
    plotdt$ts <- factor(plotdt$ts, levels=c("Day1_Morning","Day1_Afternoon","Day1_Evening","Day2_Morning","Day2_Afternoon","Day2_Evening","Day3_Morning","Day3_Afternoon","Day3_Evening"))
    #plotdt$ts <- as.numeric(as.character(plotdt$ts))
    outfig0 <- paste0(outdir, "/Diurnal.deconv.epidish.REFsam.set8.diurnalchanges.majorcell.", heytoday, ".pdf")
    pdf(outfig0,height=4.5, width=7.2)
    for (ct in majorcelltype) {
        subplotdt <- subset(plotdt, Celltype==ct)
        p0 <- ggplot(subplotdt,aes(x=ts,y=meanprop,group=Sex,color=Sex))+geom_errorbar(aes(ymin=meanprop-sdprop, ymax=meanprop+sdprop),width=0.3, linewidth=0.8,position=position_dodge(0.2))+geom_line(aes(color=Sex),linewidth=1, position=position_dodge(0.2))+geom_point(position=position_dodge(0.2), size=2, shape=21, fill="white")+commontheme()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=11,color="black"))+labs(title=ct,x="",y="proportion (%)")+scale_color_manual(values=c("#EA7317","#2364AA"))
        print(p0)
    }
    dev.off()

    ## specific cell types for repeated ANOVA
    outfig3 <- paste0(outdir, "/Diurnal.deconv.epidish.REFsam.set8.timechanges.grouped_time.", heytoday, ".pdf")
    pdf(outfig3,height=5.6,width=4.5)
    monosub <- subset(dt0long,Celltype=="Monocytes")
    res.aov <- anova_test(data = monosub, dv = Proportion, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    p32 <- ggplot(monosub, aes(x=Time,y=Proportion))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="Monocytes",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+annotate(geom="text", x=2, y=1.05*max(monosub$Proportion), label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p32)

    nksub <- subset(dt0long,Celltype=="NK-cells")
    res.aov <- anova_test(data = nksub, dv = Proportion, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    p30 <- ggplot(nksub, aes(x=Time,y=Proportion))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="NK-cells",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+annotate(geom="text", x=2, y=1.05*max(nksub$Proportion), label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p30)

    hepsub <- subset(dt0long,Celltype=="Liver-hepatocytes")
    res.aov <- anova_test(data = hepsub, dv = Proportion, wid = SubjectID, within = c(Day, Time))
    get_anova_table(res.aov)
    p31 <- ggplot(hepsub, aes(x=Time,y=Proportion))+geom_boxplot(lwd=1, aes(color=Time),width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="Liver-hepatocytes",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+annotate(geom="text", x=2, y=1.05*max(hepsub$Proportion), label=paste("Repeated~measures~~ANOVA:~italic(P)==",res.aov$ANOVA$p[2]),parse=TRUE, color="black",size=4)
    print(p31)
    dev.off()

    ###stratified by day
    #plotdt0_stat <- plotdt0 %>% group_by(Day, Celltype) %>% pairwise_t_test(Proportion ~ Time, paired=TRUE, p.adjust.method="bonferroni") %>% add_xy_position(x = "Time")
    outfig3 <- paste0(outdir, "/Diurnal.deconv.epidish.REFsam.set8.timechanges.stratify_day.", heytoday, ".pdf")
    pdf(outfig3,height=4.2,width=7.8)
    monosub <- subset(dt0long,Celltype=="Monocytes")
    stat.test32 <- monosub %>% group_by(Day) %>% pairwise_t_test(Proportion ~ Time, paired=TRUE,p.adjust.method="bonferroni") %>% add_xy_position(x="Time")
    p32 <- ggplot(monosub, aes(x=Time,y=Proportion))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="Monocytes",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+facet_grid(.~Day)+stat_pvalue_manual(stat.test32,label="p.adj",size=3.5,step.increase=0.05,vjust=0,tip.length=0)
    print(p32)

    nksub <- subset(dt0long,Celltype=="NK-cells")
    stat.test30 <- nksub %>% group_by(Day) %>% pairwise_t_test(Proportion ~ Time, paired=TRUE,p.adjust.method="bonferroni") %>% add_xy_position(x="Time")
    p30 <- ggplot(nksub, aes(x=Time,y=Proportion))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="NK-cells",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+facet_grid(.~Day)+stat_pvalue_manual(stat.test30,label="p.adj",size=3.5,step.increase=0.05,vjust=0,tip.length=0)
    print(p30)

    hepsub <- subset(dt0long,Celltype=="Liver-hepatocytes")
    stat.test31 <- hepsub %>% group_by(Day) %>% pairwise_t_test(Proportion ~ Time, paired=TRUE,p.adjust.method="bonferroni") %>% add_xy_position(x="Time")
    p31 <- ggplot(hepsub, aes(x=Time,y=Proportion))+geom_boxplot(color="black",width=0.6,outlier.shape=NA)+geom_jitter(aes(color=Time),position=position_jitter(0.2))+commontheme()+theme(axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=11,color="black"), legend.position="none")+scale_color_manual(values=c("#EFC000","darkred","royalblue4"))+labs(title="Liver-hepatocytes",x="",y="Cell type proportion (%)")+scale_y_continuous(expand = c(0.05, 0.25))+facet_grid(.~Day)+stat_pvalue_manual(stat.test31,label="p.adj",size=3.5,step.increase=0.05,vjust=0,tip.length=0)
    print(p31)
    dev.off()

    ###cell type diurnal rhythm analysis
    plotdt0$decimaltime <- sapply(plotdt0$SamplingTime, function(x) {
       parts <- strsplit(x, ":")[[1]]
       as.numeric(parts[1]) + as.numeric(parts[2]) / 60
    })
    plotdt0$intday <- as.numeric(as.factor(plotdt0$Day))

    for (ct in unique(plotdt0$Celltype)) {
        subplotdt <- subset(plotdt0, Celltype==ct)
        print(paste("Evaluating diurnal rhythm for cell type:", ct))
        #fit <- cosinor.lm(Proportion ~ time(decimaltime), period = 24, subplotdt, na.action = na.omit)
        #print(summary(fit))
        fit2 <- cosinor.lm(Proportion ~ time(decimaltime) + intday, period = 24, subplotdt, na.action = na.omit)
        print(summary(fit2))
    }

    
    day_map <- setNames(0:2, c("Day1", "Day2", "Day3"))
    plotdt0$day_num <- day_map[as.character(plotdt0$Day)]
    plotdt0$decimal_time <- sapply(plotdt0$SamplingTime, function(x) {
    parts <- strsplit(x, ":")[[1]]
        as.numeric(parts[1]) + as.numeric(parts[2]) / 60
    })
    plotdt0$continuous_time <- plotdt0$day_num * 24 + plotdt0$decimal_time

    outfig0 <- paste0(outdir, "/Diurnal.deconv.epidish.REFsam.set8.diurnalchanges.spec_cosinor.", heytoday, ".pdf")
    cosinorct <- c("Monocytes","NK-cells","Liver-hepatocytes")
    subplotdt0 <- subset(plotdt0, Celltype=="Monocytes")
    fit0 <- cosinor.lm(Proportion ~ time(continuous_time) + intday, period = 24, subplotdt0, na.action = na.omit)
    print(summary(fit0))
    subplotdt1 <- subset(plotdt0, Celltype=="NK-cells")
    fit1 <- cosinor.lm(Proportion ~ time(continuous_time) + intday, period = 24, subplotdt1, na.action = na.omit)
    print(summary(fit1))
    subplotdt2 <- subset(plotdt0, Celltype=="Liver-hepatocytes")
    fit2 <- cosinor.lm(Proportion ~ time(continuous_time) + intday, period = 24, subplotdt2, na.action = na.omit)
    print(summary(fit2))
    pdf(outfig0,height=4, width=7)
    p0 <- ggplot(subplotdt0,aes(x=continuous_time,y=Proportion))+geom_line(aes(group=SubjectID),color="gray",linewidth=0.6)+geom_point(size=1, shape=21, fill="white")+commontheme()+theme(legend.position="none")+labs(title="Monocytes",x="time (h)",y="proportion (%)")+scale_color_manual(values=rainbow(16))+stat_smooth(method="lm", formula=y ~ cos(2 * pi * x / 24) + sin(2 * pi * x / 24), se=FALSE, color="black", linewidth=1)+scale_x_continuous(breaks=seq(8, 72, 8), label=c("Day1\n08:00","Day1\n16:00","Day1\n24:00","Day2\n08:00","Day2\n16:00","Day2\n24:00","Day3\n08:00","Day3\n16:00","Day3\n24:00"))+geom_text(x=10, y=0.95*max(subplotdt0$Proportion), label=paste0("Cosinor fit: period = 24h\n", "amplitude p-value = 0.0056"), size=3.5, color="black", hjust=0)+scale_y_continuous(expand = c(0.05, 0.25))
    print(p0)
    p1 <- ggplot(subplotdt1,aes(x=continuous_time,y=Proportion))+geom_line(aes(group=SubjectID),color="gray",linewidth=0.6)+geom_point(size=1, shape=21, fill="white")+commontheme()+theme(legend.position="none")+labs(title="NK-cells",x="time (h)",y="proportion (%)")+scale_color_manual(values=rainbow(16))+stat_smooth(method="lm", formula=y ~ cos(2 * pi * x / 24) + sin(2 * pi * x / 24), se=FALSE, color="black", linewidth=1)+scale_x_continuous(breaks=seq(8, 72, 8), label=c("Day1\n08:00","Day1\n16:00","Day1\n24:00","Day2\n08:00","Day2\n16:00","Day2\n24:00","Day3\n08:00","Day3\n16:00","Day3\n24:00"))+geom_text(x=10, y=0.95*max(subplotdt1$Proportion), label=paste0("Cosinor fit: period = 24h\n", "amplitude p-value = 0.0003"), size=3.5, color="black", hjust=0)+scale_y_continuous(expand = c(0.05, 0.25))
    print(p1)
    p2 <- ggplot(subplotdt2,aes(x=continuous_time,y=Proportion))+geom_line(aes(group=SubjectID),color="gray",linewidth=0.6)+geom_point(size=1, shape=21, fill="white")+commontheme()+theme(legend.position="none")+labs(title="Liver-hepatocytes",x="time (h)",y="proportion (%)")+scale_color_manual(values=rainbow(16))+stat_smooth(method="lm", formula=y ~ cos(2 * pi * x / 24) + sin(2 * pi * x / 24), se=FALSE, color="black", linewidth=1)+scale_x_continuous(breaks=seq(8, 72, 8), label=c("Day1\n08:00","Day1\n16:00","Day1\n24:00","Day2\n08:00","Day2\n16:00","Day2\n24:00","Day3\n08:00","Day3\n16:00","Day3\n24:00"))+geom_text(x=10, y=0.95*max(subplotdt2$Proportion), label=paste0("Cosinor fit: period = 24h\n", "amplitude p-value = 0.0011"), size=3.5, color="black", hjust=0)+scale_y_continuous(expand = c(0.05, 0.25))
    print(p2)
    dev.off()


}
