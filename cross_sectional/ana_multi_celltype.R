#!/usr/bin/env Rscript
## This script is used to generate plots for the cross-sectional cohort (celltype, size, motif)

## Load necessary libraries
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "circlize", "ComplexHeatmap", "dendsort", "CCP", "CCA", "cowplot","pheatmap","viridis")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(packages, library, character.only = TRUE))
heytoday <- Sys.Date()

args <- commandArgs(trailingOnly = TRUE)
deconvresfh <- args[1]
sizemotiffh <- args[2]
motifh <- args[3]
sizestatsfh <- args[4]
outdir <- args[5]
kmer <- "2mer"


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

GetCorrMat <- function(metadt1,smmdt1) {
    cormat <- data.frame()
    for (mycell in majorcelltypes) {
      for (mysizemotif in mycombl) {
         print(paste0("Cell: ", mycell, ";MotifSize: ", mysizemotif))
         submetadt1 <- metadt1[,colnames(metadt1) %in% c("GCcode",mycell)]
         subsmmdt1 <- smmdt1[smmdt1$Motif_Size==mysizemotif,c("GCcode","Size_Motif_Freq")]
         subdt <- merge(submetadt1,subsmmdt1,by=c("GCcode"),all.x=TRUE)
         subdt[is.na(subdt)] <- 0
         mycor <- cor(subdt[,colnames(subdt)==mycell],subdt$Size_Motif_Freq, method="spearman")
         samcor <- data.frame(Cell=mycell, Motif_Size=mysizemotif, Scorr=mycor)
         cormat <- rbind(cormat, samcor)
      }
    }
    return(cormat)
}

plotheatmap2mercluster <- function(metadt, smmdt, inYear) {
    tcor <-  GetCorrMat(metadt,smmdt)
    tcordt <- separate_wider_delim(tcor, cols=Motif_Size, delim = "_", names = c("Motif","Size"))
    tcordt$Size <- factor(tcordt$Size, levels=c(50:450))
    tcordt$Cell <- factor(tcordt$Cell, levels=majorcelltypes)
    #tcordt$Motif <- factor(tcordt$Motif, levels=c("C","G","A","T"))
    wtcordt <- spread(tcordt,Size,Scorr)
    wtcormat <- as.matrix(wtcordt[,!(colnames(wtcordt) %in% c("Cell","Motif"))])
    rownames(wtcormat) <- wtcordt$Motif

    panel_fun = function(index, nm) {
        grid.rect(gp=gpar(fill=c("#737373")))
        grid.text(paste0(nm, " bp"), 0.5, 0.5, gp=gpar(col="white",fontsize = 12))
    }
    align_to = list("50-150"=1:101, "151-250"=102:201, "251-350"=202:301, "351-450"=302:401)
    va <- HeatmapAnnotation(foo=anno_block(align_to=align_to, panel_fun = panel_fun),bar=anno_mark(at=c(1,96,118, 138, 231, 281, 371),labels=c("50bp","145bp","167bp","187bp","280bp","330bp","420bp"), which="column", side="bottom"))
    
    #col_fun = colorRamp2(c(min(wtcormat),min(wtcormat)/2,0,max(wtcormat)/2,max(wtcormat)),c("#0571b0","#92c5de","#f7f7f7","#f4a582","#ca0020"))
    #col_fun = colorRamp2(c(min(wtcormat),0,max(wtcormat)),c("blue","#EEEEEE","red"))
    col_fun = colorRamp2(c(-0.42,0,0.42),c("blue","white","red"))

    for (i in 1:length(majorcelltypes)) {
        subwtcordt <- subset(wtcordt, Cell==majorcelltypes[i])
        subwtcormat <- as.matrix(subwtcordt[,!(colnames(subwtcordt) %in% c("Cell","Motif"))])
        rownames(subwtcormat) <- subwtcordt$Motif
        
        subendnucdt <- data.frame(Motif=subwtcordt$Motif)
        subendnuc <- subendnucdt %>% separate_wider_position(Motif, c(EndNuc1=1, EndNuc2=1))
        ha <- rowAnnotation(EndNuc1=subendnuc$EndNuc1, EndNuc2=subendnuc$EndNuc2, simple_anno_size=unit(0.8, "cm"),col=list(EndNuc1=c("C"="#FFBB00","G"="#3F681C","A"="#FB6542","T"="#375E97"),EndNuc2=c("C"="#FFBB00","G"="#3F681C","A"="#FB6542","T"="#375E97")),show_annotation_name=FALSE)

        row_dend = dendsort(hclust(dist(subwtcormat)))

        if (i == 1) {
            p0 <- Heatmap(subwtcormat,name="cor_mat", right_annotation=ha, heatmap_legend_param=list(at=c(round(min(wtcormat),2), 0, round(max(wtcormat),2))),cluster_columns=FALSE, cluster_rows = row_dend, show_row_names=TRUE,show_column_names=FALSE, row_title=majorcelltypes[i], row_names_gp=gpar(fontsize=6),row_title_rot = 0, col=col_fun)
            ht_list = p0
        } else {
            #p0 <- Heatmap(subwtcormat,name="cor_mat", right_annotation=ha, heatmap_legend_param=list(at=c(round(min(wtcormat),2), 0, round(max(wtcormat),2))),cluster_columns=FALSE, cluster_rows=TRUE, row_dend_reorder=TRUE, show_row_names=TRUE,show_column_names=FALSE,row_names_gp=gpar(fontsize=6),row_split=subwtcordt$Cell,row_title_rot = 0, col=col_fun,  show_heatmap_legend = FALSE)
            p0 <- Heatmap(subwtcormat,name="cor_mat", right_annotation=ha, heatmap_legend_param=list(at=c(round(min(wtcormat),2), 0, round(max(wtcormat),2))),cluster_columns=FALSE, cluster_rows = row_dend, show_row_names=TRUE,show_column_names=FALSE,row_title=majorcelltypes[i], row_names_gp=gpar(fontsize=6),row_title_rot = 0, col=col_fun,  show_heatmap_legend = FALSE)
            ht_list = ht_list %v% p0
        }
    }

    outfig0 <- paste0(outdir, "/LONGVAR.celltype.cor.",kmer,".sizemotif.",inYear,".clustered.heatmap.",heytoday,".pdf")
    pdf(outfig0, height=10.4, width=11.8)
    ht_list = ht_list %v% va
    print(ht_list)
    dev.off()

    #p0 <- Heatmap(wtcormat,name="cor_mat",bottom_annotation=va, right_annotation=ha,cluster_columns=FALSE, cluster_rows=FALSE, show_row_names=TRUE,show_column_names=FALSE,row_names_gp=gpar(fontsize=6),row_split=wtcordt$Cell,row_title_rot = 0,heatmap_legend_param=list(at=c(round(min(wtcormat),2), 0, round(max(wtcormat),2))))
}

plot_sizemotif_celltype_cor <- function(deconvresfh, sizemotiffh, outdir) {
    metadt <- fread(deconvresfh,header=TRUE,sep="\t",data.table=FALSE)
    metadt$MEP <- metadt$Megakaryocytes+metadt$`Erythrocyte-progenitors`

    metadt1 <- subset(metadt, Year=="T1")
    metadt2 <- subset(metadt, Year=="T2")

    majorcelltypes <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    smetadt <- metadt[,colnames(metadt) %in% c("GCcode","Year")]

    sizemotif <- fread(sizemotiffh,header=TRUE,sep="\t",data.table=FALSE)

    mymotifcomb <- expand.grid(c("C","G","T","A"),c("C","G","T","A"))
    mymotifcombl <- unlist(paste0(mymotifcomb$Var1,mymotifcomb$Var2))

    sizemotif$Motif_Size <- paste0(sizemotif$Motif, "_", sizemotif$Size)
    sizemotif <- sizemotif[,c("GCcode","Motif_Size","Size_Motif_Freq")]

    smmdt <- merge(sizemotif,smetadt,by=c("GCcode"),all.x=TRUE)
    smmdt1 <- subset(smmdt, Year=="T1")
    smmdt2 <- subset(smmdt, Year=="T2")

    mycomb <- expand.grid(c("C","G","T","A"),c("C","G","T","A"),c(50:450))
    mycombl <- unlist(paste0(mycomb$Var1,mycomb$Var2,"_",mycomb$Var3))

    plotheatmap2mercluster(metadt=metadt1, smmdt=smmdt1, inYear="T1")
    plotheatmap2mercluster(metadt=metadt2, smmdt=smmdt2, inYear="T2")
}

perform_CCA_celltype_endmotif <- function(deconvresfh, motifh, outdir) {
    motifdt <- fread(motifh, header=TRUE, sep="\t", data.table=FALSE)
    metadt <- fread(deconvresfh, header=TRUE, sep="\t", data.table=FALSE)
    metadt$MEP <- metadt$Megakaryocytes+metadt$`Erythrocyte-progenitors`

    majorcelltype <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    cellmeta <- metadt[,c("GCcode","Year","Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")]

    motif2celldt <- merge(cellmeta, motifdt, by=c("GCcode"))

    mymotifcomb <- expand.grid(c("C","G","T","A"),c("C","G","T","A"))
    mymotifcombl <- unlist(paste0(mymotifcomb$Var1,mymotifcomb$Var2))

    motifmat <- scale(motif2celldt[,colnames(motif2celldt) %in% mymotifcombl])
    cellmat <- scale(motif2celldt[,colnames(motif2celldt) %in% majorcelltype])

    lambdas <- estim.regul(cellmat, motifmat, grid1=seq(0.001,1,length=50),grid2=seq(0.001,1,length=50))
    lambda1 <- lambdas$lambda1
    lambda2 <- lambdas$lambda2

    cca_reg <- rcc(cellmat, motifmat, lambda1=lambda1, lambda2=lambda2)
    ccacors <- data.frame(CCAcor=cca_reg$cor,Dim=factor(c(1:length(cca_reg$cor))))

    outfig30 <- paste0(outdir,"/LONGVAR.Cellfrac_endmotif.2mer.rCCAs.bar.all288sam.",heytoday,".pdf")
    pdf(outfig30,height=4,width=5.8)
    p30 <- ggplot(ccacors, aes(x=Dim, y=CCAcor))+geom_bar(stat="identity",fill="white",color="black",width=0.7)+commontheme()+labs(x="Dimensions", y="Canonical correlations")+geom_text(aes(label=round(CCAcor,digits=2)), vjust=1.6, color="black", size=3.5)+scale_y_continuous(limits=c(0, 0.8), breaks=seq(0,0.8,by=0.15))
    print(p30)
    dev.off()

    set.seed(98765)
    n_permutations <- 1000
    permuted_cor1 <- permuted_cor2 <- permuted_cor3 <- permuted_cor4 <- numeric(n_permutations)

    for (i in 1:n_permutations) {
        shuffled_x <- cellmat[sample(1:nrow(cellmat)), ]
        permuted_result <- rcc(shuffled_x, motifmat, lambda1 = lambda1, lambda2 = lambda2)
        permuted_cor1[i] <- permuted_result$cor[1]
        permuted_cor2[i] <- permuted_result$cor[2]
        permuted_cor3[i] <- permuted_result$cor[3]
        permuted_cor4[i] <- permuted_result$cor[4]
    }
    p_value1 <- mean(permuted_cor1 >= cca_reg$cor[1])
    p_value2 <- mean(permuted_cor2 >= cca_reg$cor[2])
    p_value3 <- mean(permuted_cor3 >= cca_reg$cor[3])
    p_value4 <- mean(permuted_cor4 >= cca_reg$cor[4])
    print(paste("P-value for first canonical correlation:", p_value1))

    cell_loadings <- cca_reg$xcoef
    motif_loadings <- cca_reg$ycoef

    cellloadingdt <- data.frame(cca_reg$xcoef)
    cellloadingdt <- cellloadingdt[,c("X1","X2")]
    colnames(cellloadingdt) <- c("CCA1","CCA2")
    cellloadingdt$Celltype <- rownames(cellloadingdt)
    cellloadingdt$Celltype <- factor(cellloadingdt$Celltype,levels=rev(majorcelltype))
    cellloadingdt$ccolorgrp1 <- ifelse(cellloadingdt$CCA1>0, "Color1", "Color2")
    cellloadingdt$ccolorgrp2 <- ifelse(cellloadingdt$CCA2>0, "Color1", "Color2")

    motifloadingdt <- data.frame(cca_reg$ycoef)
    motifloadingdt <- motifloadingdt[,c("X1","X2")]
    colnames(motifloadingdt) <- c("CCA1","CCA2")
    motifloadingdt$Motif <- rownames(motifloadingdt)
    motifloadingdt$mcolorgrp1 <- ifelse(motifloadingdt$CCA1>0, "Color1", "Color2")
    motifloadingdt$mcolorgrp2 <- ifelse(motifloadingdt$CCA2>0, "Color1", "Color2")
    ylim1 <- max(abs(motifloadingdt$CCA1))+0.1
    ylim2 <- max(abs(motifloadingdt$CCA2))+0.1


    outfig31 <- paste0(outdir,"/LONGVAR.Cellfrac_endmotif.2mer.rCCA_loadings.all288sam.dim1.",heytoday,".pdf")
    pdf(outfig31,height=4.2,width=9)
    p1 <- ggplot(cellloadingdt)+geom_col(aes(CCA1, Celltype, fill=ccolorgrp1),width=0.7,color="black")+commontheme()+theme(axis.text.y=element_text(size=10,color="black"),axis.text.x=element_text(size=12,color="black"),legend.position="none")+labs(x="CCA1", y="Cell type")+scale_fill_manual(values=c("#762a83","#1b7837"))
    p2 <- ggplot(motifloadingdt,aes(x=Motif, y=CCA1, fill=mcolorgrp1))+geom_bar(stat="identity",color="black",width=0.7)+commontheme()+theme(axis.text.x=element_text(size=11,color="black"),axis.text.y=element_text(size=12,color="black"),legend.position="none")+labs(x="2mer end motif", y="CCA1")+scale_fill_manual(values=c("#EFB509","#002C54"))+ylim(-ylim1,ylim1)
    pp <- plot_grid(p1, p2, labels = c('',''), ncol=2, rel_widths=c(1, 1.5))
    print(pp)
    dev.off()

    outfig32 <- paste0(outdir,"/LONGVAR.Cellfrac_endmotif.2mer.rCCA_loadings.all288sam.dim2.",heytoday,".pdf")
    pdf(outfig32,height=4.2,width=9)
    p3 <- ggplot(cellloadingdt)+geom_col(aes(CCA2, Celltype, fill=ccolorgrp2),width=0.7,color="black")+commontheme()+theme(axis.text.y=element_text(size=10,color="black"),axis.text.x=element_text(size=12,color="black"),legend.position="none")+labs(x="CCA2", y="Cell type")+scale_fill_manual(values=c("#762a83","#1b7837"))
    p4 <- ggplot(motifloadingdt,aes(x=Motif, y=CCA2, fill=mcolorgrp2))+geom_bar(stat="identity",color="black",width=0.7)+commontheme()+theme(axis.text.x=element_text(size=11,color="black"),axis.text.y=element_text(size=12,color="black"),legend.position="none")+labs(x="2mer end motif", y="CCA2")+scale_fill_manual(values=c("#EFB509","#002C54"))+ylim(-ylim2,ylim2)
    pp <- plot_grid(p3, p4, labels = c('',''), ncol=2, rel_widths=c(1, 1.5))
    print(pp)
    dev.off()


    U <- as.data.frame(as.matrix(cellmat) %*% cca_reg$xcoef)
    V <- as.data.frame(as.matrix(motifmat) %*% cca_reg$ycoef)

    cca_plot_dt1 <- data.frame(U1 = U[,1], V1 = V[,1])
    cca_plot_dt2 <- data.frame(U2 = U[,2], V2 = V[,2])

    outfig3 <- paste0(outdir,"/LONGVAR.Cellfrac_endmotif.2mer.rCCA_dim.all288sam",heytoday,".pdf")
    pdf(outfig3,height=6,width=6)
    p3 <- ggplot(cca_plot_dt1, aes(x=U1, y=V1))+geom_point(color="black",alpha=0.7)+geom_smooth(method="lm",color="darkblue",se=FALSE, linetype="dashed")+commontheme()+labs(x="Canonical variable 1 from cell fractions", y="Canonical variable 1 from 2mer end motif")
    print(p3)
    p3 <- ggplot(cca_plot_dt2, aes(x=U2, y=V2))+geom_point(color="black",alpha=0.7)+geom_smooth(method="lm",color="darkblue",se=FALSE, linetype="dashed")+commontheme()+labs(x="Canonical variable 2 from cell fractions", y="Canonical variable 2 from 2mer end motif")
    print(p3)
    dev.off()
}

####cell-type-specific fragments
#plot selected markers for cell-type-specific fragments
rawmarkerfh <- "set8.averaged.tsv"
plot_ct_marker <- function(outdir) {
    rawmarker <- fread(rawmarkerfh, header=TRUE, sep="\t", data.table=FALSE)

    rawdt <- rawmarker[,colnames(rawmarker) %in% c("chr","start","end","Granulocytes","Monocytes","Meagkaryocytes","Erythrocyte-progenitors","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")]
    rawdt$MEP <- (rawmarker$Megakaryocytes+rawmarker$`Erythrocyte-progenitors`)/2
    rawdt <- rawdt[,! (colnames(rawdt) %in% c("Megakaryocytes","Erythrocyte-progenitors"))]

    majorcelltypes <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    selmarkerdt <- data.frame()
    for (i in 1:length(majorcelltypes)) {
        print(majorcelltypes[i])
        for (j in 1:nrow(rawdt)) {
            ctdiff1 <- rawdt[j,c(majorcelltypes[i])]-max(rawdt[j,!(colnames(rawdt) %in% c("chr","start","end",majorcelltypes[i]))])
            ctdiff2 <- rawdt[j,c(majorcelltypes[i])]-min(rawdt[j,!(colnames(rawdt) %in% c("chr","start","end",majorcelltypes[i]))])
            
            if (ctdiff1 < -0.6 && ctdiff2 < -0.6 ) {
                addrow <- data.frame(rawdt[j,],Celltype=majorcelltypes[i],Mestate="Hypo",Diff=ctdiff2)
            } else if (ctdiff1 > 0.6 && ctdiff2 > 0.6) {
                addrow <- data.frame(rawdt[j,],Celltype=majorcelltypes[i],Mestate="Hyper",Diff=ctdiff1)
            } else {
                next
            }

            if (nrow(selmarkerdt) == 0) {
                selmarkerdt <- addrow
            } else {
                selmarkerdt <- rbind(selmarkerdt,addrow)
            }
        }
    }

    table(selmarkerdt[,c("Celltype","Mestate")])
    names(selmarkerdt) <- gsub(x = names(selmarkerdt), pattern = "\\.", replacement = "-")

    plotdt <- selmarkerdt
    col_order <- c("chr","start","end","Diff","Mestate","Celltype","Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    plotdt <- plotdt[,col_order]

    majorcelltypes <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    col_fun = colorRamp2(c(0,0.25,0.5,0.75,1),c("#440154FF","#443A83FF","#21908CFF","#8FD744FF","#FDE725FF"))
    ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
    

    for (i in 1:length(majorcelltypes)) {
        subdt <- subset(plotdt, Celltype==majorcelltypes[i] & Mestate=="Hypo")
        submat <- as.matrix(subdt[,!(colnames(subdt) %in% c("chr","start","end","Pval","Diff","Mestate","Prank","DiffRank","Celltype"))])
        rownames(submat) <- paste0(subdt$chr,":",subdt$start,"-",subdt$end)
        ha <- rowAnnotation(Marker=subdt$Mestate, simple_anno_size=unit(0.7, "cm"),col=list(Marker=c("Hypo"="darkblue","Hyper"="darkred")),show_annotation_name=FALSE)
        if (i == 1) {
            p0 <- Heatmap(submat,name="Meth_level",right_annotation=ha,cluster_columns=FALSE, cluster_rows=FALSE, show_row_names=FALSE,show_column_names=TRUE, row_title=majorcelltypes[i], row_title_gp = gpar(fontsize=8),row_names_gp=gpar(fontsize=3.5),row_title_rot = 0, col=col_fun)
            ht_list0 = p0
        } else {
            p0 <- Heatmap(submat,name="Meth_level",right_annotation=ha,cluster_columns=FALSE, cluster_rows=FALSE, show_row_names=FALSE,show_column_names=TRUE, row_title=majorcelltypes[i], row_title_gp = gpar(fontsize=8), row_names_gp=gpar(fontsize=3.5),row_title_rot = 0, col=col_fun,  show_heatmap_legend = FALSE)
            ht_list0 = ht_list0 %v% p0
        }
    }

    for (i in 1:length(majorcelltypes)) {
        subdt <- subset(plotdt, Celltype==majorcelltypes[i] & Mestate=="Hyper")
        if (nrow(subdt) == 0) {
            j=1
        } else {
            submat <- as.matrix(subdt[,!(colnames(subdt) %in% c("chr","start","end","Pval","Diff","Mestate","Prank","DiffRank","Celltype"))])
            rownames(submat) <- paste0(subdt$chr,":",subdt$start,"-",subdt$end)
            ha <- rowAnnotation(Marker=subdt$Mestate, simple_anno_size=unit(0.7, "cm"),col=list(Marker=c("Hypo"="darkblue","Hyper"="darkred")),show_annotation_name=FALSE)
            if (j == 1) {
                p0 <- Heatmap(submat,name="Meth_level",right_annotation=ha,cluster_columns=FALSE, cluster_rows=FALSE, show_row_names=FALSE,show_column_names=TRUE, row_title=majorcelltypes[i], row_title_gp = gpar(fontsize=8),row_names_gp=gpar(fontsize=3.5),row_title_rot = 0, col=col_fun)
                ht_list1 = p0
            } else {
                p0 <- Heatmap(submat,name="Meth_level",right_annotation=ha,cluster_columns=FALSE, cluster_rows=FALSE, show_row_names=FALSE,show_column_names=TRUE, row_title=majorcelltypes[i], row_title_gp = gpar(fontsize=8), row_names_gp=gpar(fontsize=3.5),row_title_rot = 0, col=col_fun,  show_heatmap_legend = FALSE)
                ht_list1 = ht_list1 %v% p0
            }
            j = j + 1
        }
    }
    outfig0 <- paste0(outdir, "/selected.cell_type_specific.",heytoday,".pdf")
    pdf(outfig0, height=8.9, width=6.8)
    ht_list <- ht_list0 %v% ht_list1
    print(ht_list)
    dev.off()
}


indir <- "/cellfragment/pool"
Getsizeprofile <- function(indir,inYear,outdir) {
    refsizefh <- paste0(indir,"/",inYear,"pooled.size.count.tsv")
    infh0 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hypo.bed")
    #infh1 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hypo.bed")
    #infh2 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hyper.bed")
    infh3 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hyper.bed")
    refsize <- fread(refsizefh,header=TRUE,sep="\t",data.table=FALSE)
    refsize$RefFreq <- 100*refsize$PoolCount/sum(refsize$PoolCount)

    pdt0 <- fread(infh0,header=FALSE,sep="\t",data.table=FALSE)
    pdt3 <- fread(infh3,header=FALSE,sep="\t",data.table=FALSE)
    colnames(pdt0) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    colnames(pdt3) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    
    pdt <- pdt0
    fpdt <- pdt %>% group_by(Size, Mestate, Celltype) %>% summarise(Count=n()) %>% group_by(Mestate,Celltype) %>% mutate(Freq=100*Count/sum(Count))
    totcnt <- pdt %>% group_by(Mestate,Celltype) %>% summarise(Count=n())

    majorcelltype <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")
    fpdt$Celltype <- factor(fpdt$Celltype, levels=majorcelltype)
    mytext <- data.frame(Celltype=character(), mylabel=character())
    for (i in majorcelltype) {
        subtotcnt <- subset(totcnt, Celltype==i)
        curtext <- data.frame(Celltype=i, mylabel=paste0("Hypo=",subtotcnt[subtotcnt$Mestate=="Hypo",]$Count), Mestate="Hypo")
        mytext <- rbind(mytext,curtext)
    }
    mytext$Celltype <- factor(mytext$Celltype, levels=majorcelltype)
  
    outfig <- paste0(outdir,"/",inYear,"d60pool.fragments.sizeprofile.hypo0_cellspec.",heytoday,".pdf")
    pdf(outfig,height=7.8,width=8.8)
    p0 <- ggplot(fpdt,aes(x=Size,y=Freq))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(color="darkblue",linewidth=0.8,alpha=0.7)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=10,color="black"),legend.text=element_text(size=11),axis.title.x=element_text(size=14,color="black",vjust=-0.5),axis.title.y=element_text(size=14,color="black",vjust=1.5),legend.position=c(0.85,0.15),strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency (%)",color="Marker")+scale_x_continuous(limits=c(49, 451),breaks=seq(50, 450, by=50))+scale_color_manual(values=c("darkblue","darkred"))+facet_wrap(~Celltype, ncol=3)+geom_text(data=mytext,aes(label=mylabel, x = Inf, y = Inf), color="black",hjust=1.05,vjust=1.5)+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
    print(p0)
    dev.off()
}

Getsizeprofile(indir,inYear="T1",outdir)
Getsizeprofile(indir,inYear="T2",outdir)

Getsizeprofile1 <- function(indir,inYear,outdir) {
    #infh0 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hypo.bed")
    #infh1 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hypo.bed")
    infh2 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hyper.bed")
    infh3 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hyper.bed")
    #pdt0 <- fread(infh0,header=FALSE,sep="\t",data.table=FALSE)
    #pdt1 <- fread(infh1,header=FALSE,sep="\t",data.table=FALSE)
    pdt2 <- fread(infh2,header=FALSE,sep="\t",data.table=FALSE)
    pdt3 <- fread(infh3,header=FALSE,sep="\t",data.table=FALSE)
    
    #colnames(pdt0) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    #colnames(pdt1) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    colnames(pdt2) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    colnames(pdt3) <- c("Chrom","Start","End","Size","Mestate","Celltype")

    pdt2 <- subset(pdt2, Celltype=="Liver-hepatocytes")
    pdt3 <- subset(pdt3, Celltype=="Liver-hepatocytes")
    #pdt0$Group <- "Cellspec"
    #pdt1$Group <- "Unspec"
    pdt2$Group <- "Others"
    pdt3$Group <- "Liver-hepatocytes"
    pdt <- rbind(pdt2,pdt3)
    totcnt <- pdt %>% group_by(Group) %>% summarise(Count=n())
    fpdt <- pdt %>% group_by(Size, Group, Mestate, Celltype) %>% summarise(Count=n()) %>% group_by(Group, Mestate,Celltype) %>% mutate(Freq=100*Count/sum(Count))
    outfig <- paste0(outdir,"/",inYear,"d60pool.fragments.sizeprofile.liver_hyperregion_spec_vs_unspec.",heytoday,".pdf")
    pdf(outfig,height=4,width=5.8)
    p0 <- ggplot(fpdt,aes(x=Size,y=Freq,color=Group))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=0.8,alpha=0.9)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=12,color="black"),legend.text=element_text(size=12),axis.title.x=element_text(size=14,color="black",vjust=-0.5),axis.title.y=element_text(size=14,color="black",vjust=1.5),legend.position="top",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency (%)",color="Cell type")+scale_x_continuous(limits=c(49, 451),breaks=seq(50, 450, by=50))+scale_color_manual(values=c("darkred","gray50"))+geom_text(label=paste0("Liver-hepatocytes: ",totcnt[totcnt$Group=="Liver-hepatocytes",]$Count,"\nOthers: ",totcnt[totcnt$Group!="Liver-hepatocytes",]$Count), x=Inf, y=Inf, color="black",hjust=1.05, vjust=1.5)+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
    print(p0)
    dev.off()
    #return(fpdt)
}


Getsizeprofile2 <- function(indir,inYear,outdir) {
    infh0 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hypo.bed")
    #infh1 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hypo.bed")
    #infh2 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.um_hyper.bed")
    infh3 <- paste0(indir,"/",inYear,"d60pool.fragments.cellspecific.cpgsite.me_hyper.bed")
    pdt0 <- fread(infh0,header=FALSE,sep="\t",data.table=FALSE)
    pdt3 <- fread(infh3,header=FALSE,sep="\t",data.table=FALSE)
    colnames(pdt0) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    colnames(pdt3) <- c("Chrom","Start","End","Size","Mestate","Celltype")
    #pdt0 <- pdt0[pdt0$Size >= 50 & pdt0$Size <=450,]
    
    #cellspecdt <- rbind(pdt0,pdt3)
    #cellspecdt <- cellspecdt[,c("Chrom","Start","End","Size","Celltype")]
    #cellspecdt <- pdt3
    cellspecdt <- pdt0
    majorcelltype <- c("Granulocytes","Monocytes","MEP","NK-cells","T-cells-CD4","B-cells","Liver-hepatocytes","Vascular-endothelium")

    outfig <- paste0(outdir,"/",inYear,"d60pool.fragments.one_vs_others.hypo0_specific.",heytoday,".pdf")
    pdf(outfig,height=4.2,width=5.8)
    for (i in majorcelltype) {
        pdt <- cellspecdt
        pdt$Celltype <- ifelse (pdt$Celltype==i, i, "Others")
        fpdt <- pdt %>% group_by(Size, Celltype) %>% summarise(Count=n()) %>% group_by(Celltype) %>% mutate(Freq=100*Count/sum(Count))
        fpdt$Celltype <- factor(fpdt$Celltype, levels=c(i,"Others"))
        totcnt <- pdt %>% group_by(Celltype) %>% summarise(Count=n())
        totcnt$Celltype <- factor(totcnt$Celltype, levels=c(i,"Others"))

        p0 <- ggplot(fpdt,aes(x=Size,y=Freq,color=Celltype))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_line(linewidth=0.8,alpha=0.9)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=12,color="black"),legend.text=element_text(size=12),axis.title.x=element_text(size=14,color="black",vjust=-0.5),axis.title.y=element_text(size=14,color="black",vjust=1.5),legend.position="top",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency (%)",color="Cell type")+scale_x_continuous(limits=c(49, 451),breaks=seq(50, 450, by=50))+scale_color_manual(values=c("darkred","gray50"))+geom_text(label=paste0(i,": ",totcnt[totcnt$Celltype==i,]$Count,"\nOthers: ",totcnt[totcnt$Celltype!=i,]$Count), x=Inf, y=Inf, color="black",hjust=1.05, vjust=1.5)+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
        print(p0)
    }
    dev.off()

    outfig <- paste0(outdir,"/",inYear,"d60pool.fragments.one_vs_others.hypo0_specific.reldiff.",heytoday,".pdf")
    pdf(outfig,height=4.2,width=5.8)
    for (i in majorcelltype) {
        pdt <- cellspecdt
        pdt$Celltype <- ifelse (pdt$Celltype==i, i, "Others")
        fpdt <- pdt %>% group_by(Size, Celltype) %>% summarise(Count=n()) %>% group_by(Celltype) %>% mutate(Freq=100*Count/sum(Count))
        ifpdt <- subset(fpdt, Celltype==i)
        ofpdt <- subset(fpdt, Celltype=="Others")
        mpdt <- merge(ifpdt,ofpdt,by=c("Size"))
        mpdt$RelFreq <- mpdt$Freq.x-mpdt$Freq.y

        totcnt <- pdt %>% group_by(Celltype) %>% summarise(Count=n())
        totcnt$Celltype <- factor(totcnt$Celltype, levels=c(i,"Others"))

        p0 <- ggplot(mpdt,aes(x=Size,y=RelFreq))+geom_vline(xintercept=167,linetype="dashed",color = "gray", linewidth=1)+geom_vline(xintercept=330,linetype="dashed",color = "gray", linewidth=1)+geom_hline(yintercept=0,linetype="dashed",color = "gray", linewidth=0.5)+geom_line(color="darkred",linewidth=0.8,alpha=0.9)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(size=12,color="black"),axis.text.x=element_text(size=12,color="black"),legend.text=element_text(size=12),axis.title.x=element_text(size=14,color="black",vjust=-0.5),axis.title.y=element_text(size=14,color="black",vjust=1.5),legend.position="none",strip.text=element_text(size=12,face="bold"),strip.background=element_rect(fill="white"))+labs(x="Size (bp)",y="Frequency difference (%)\nCell-type-specific to others")+scale_x_continuous(limits=c(49, 451),breaks=seq(50, 450, by=50))+geom_text(label=paste0(i,": ",totcnt[totcnt$Celltype==i,]$Count,"\nOthers: ",totcnt[totcnt$Celltype!=i,]$Count), x=Inf, y=Inf, color="black",hjust=1.05, vjust=1.5)+annotate(x=167,y=-Inf,label="167bp",geom="text",vjust=-1,color="black")+annotate(x=330,y=-Inf,label="330bp",geom="text",vjust=-1,color="black")
        print(p0)
    }
    dev.off()
}

Getsizeprofile2(indir,inYear="T1",outdir)
Getsizeprofile2(indir,inYear="T2",outdir)
