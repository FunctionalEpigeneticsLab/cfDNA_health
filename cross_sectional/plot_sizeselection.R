#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort (size selection)

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "ggpubr", "rstatix")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

heytoday <- Sys.Date()
args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
sizeselectresfh <- args[2]
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
      , strip.background=element_rect(fill="white"))
}

Plot_SizeSelection_DeconvRes <- function(deconvresfh, sizeselectres, outdir) {
    metainfo <- fread(deconvresfh, header=TRUE, sep="\t", data.table=FALSE)
    res <- fread(sizeselectres, header=TRUE, sep="\t", data.table=FALSE)
    
    #sizeselectres$Celltype <- factor(sizeselectres$Celltype, levels=c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium"))
    #sizeselectres$Sizerange <- factor(sizeselectres$Sizerange, levels=c("All",">250bp"))
    majorcelltype <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    res$Sizerange <- factor(res$Sizerange, levels=c("All","w_selection"),labels=c("All\nfragments","With\nselection"))
    res$Group <- factor(res$Group, levels=c("<167bp","167-250bp",">250bp"))
    t1res <- subset(res, Year=="T1")
    t2res <- subset(res, Year=="T2")
    
    outfig <- paste0(outdir, "/LONGVAR.deconv.epidish.REFsam.set8.sizeselection.combined.Year1.",heytoday,".pdf")
    pdf(outfig, width=6.4, height=4.2)
    for (ct in majorcelltype) {
        subres <- subset(t1res, Celltype==ct)
        subresdt <- subres %>% group_by(Sizerange,Group) %>% mutate(N=n()) %>% mutate(Nobs=paste0('N=',N))
        stat.test <- subres %>% group_by(Group) %>% t_test(Fraction ~ Sizerange, paired=TRUE) %>% adjust_pvalue(method = "bonferroni") %>% add_xy_position(x="Sizerange")
        
        pp <- ggplot(subresdt, aes(x=Sizerange, y=Fraction))+geom_line(aes(group=Samplename), size=0.2, color="gray")+geom_boxplot(aes(fill=Sizerange),width=0.6, alpha=0.5,outlier.shape=NA)+geom_point(aes(color=Sizerange),fill="white",size=0.5)+geom_text(aes(label=Nobs),y=-Inf,vjust=-1,color="black",size=3)+commontheme()+theme(axis.text.x=element_text(size=10,color="black"),axis.text.y=element_text(size=10,color="black"),strip.background=element_blank(), strip.placement="outside", strip.text=element_text(size=12,face="bold"),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="none",plot.title=element_text(hjust=0.5))+scale_fill_manual(values=c("gray50", "darkred"))+scale_color_manual(values=c("gray50", "darkred"))+labs(title=ct, x="",y="Contribution (%)")+facet_wrap(~Group,scales="free_x",nrow=1)+guides(color=guide_legend(override.aes=list(size=2)))+stat_pvalue_manual(stat.test,label="p.adj",size=3.5,vjust=0,tip.length=0)+scale_y_continuous(expand=expansion(mult=c(0.085, .1)))
        print(pp)
    }
    dev.off()

    outfig <- paste0(outdir, "/LONGVAR.deconv.epidish.REFsam.set8.sizeselection.combined.Year2.",heytoday,".pdf")
    pdf(outfig, width=6.4, height=4.2)
    for (ct in majorcelltype) {
        subres <- subset(t2res, Celltype==ct)
        subresdt <- subres %>% group_by(Sizerange,Group) %>% mutate(N=n()) %>% mutate(Nobs=paste0('N=',N))
        stat.test <- subres %>% group_by(Group) %>% t_test(Fraction ~ Sizerange, paired=TRUE) %>% adjust_pvalue(method = "bonferroni") %>% add_xy_position(x="Sizerange")
        
        pp <- ggplot(subresdt, aes(x=Sizerange, y=Fraction))+geom_line(aes(group=Samplename), size=0.2, color="gray")+geom_boxplot(aes(fill=Sizerange),width=0.6, alpha=0.5,outlier.shape=NA)+geom_point(aes(color=Sizerange),fill="white",size=0.5)+geom_text(aes(label=Nobs),y=-Inf,vjust=-1,color="black",size=3)+commontheme()+theme(axis.text.x=element_text(size=10,color="black"),axis.text.y=element_text(size=10,color="black"),strip.background=element_blank(), strip.placement="outside", strip.text=element_text(size=12,face="bold"),axis.ticks.x=element_blank(),legend.title=element_blank(),legend.position="none",plot.title=element_text(hjust=0.5))+scale_fill_manual(values=c("gray50", "darkred"))+scale_color_manual(values=c("gray50", "darkred"))+labs(title=ct, x="",y="Contribution (%)")+facet_wrap(~Group,scales="free_x",nrow=1)+guides(color=guide_legend(override.aes=list(size=2)))+stat_pvalue_manual(stat.test,label="p.adj",size=3.5,vjust=0,tip.length=0)+scale_y_continuous(expand=expansion(mult=c(0.085, .1)))
        print(pp)
    }
    dev.off()
}
