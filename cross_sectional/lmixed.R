#!/usr/bin/env Rscript
## This script is used to perform mixed effects

## Load necessary libraries
packages <- c("data.table", "gee", "lme4", "lmerTest", "janitor")

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

test_GE_mixedeffects <- function(deconvresfh) {
    deconvres <- fread(deconvresfh, header = T, sep = "\t", data.table = F)

    deconvres$Log_Adj_Unique_Molecules <- log10(deconvres$Adj_Unique_Molecules)
    deconvres$Year <- factor(deconvres$Year, levels=c("T1", "T2"))
    deconvres$Time <- factor(deconvres$Time, levels=c("Morning","Afternoon"))
    deconvres$Sex <- factor(deconvres$Sex, levels=c("Female", "Male"))

    st1 <- lmer(Log_Adj_Unique_Molecules ~ Sex+BMI+Age+(1|SubjectID)+(1|Year)+(1|Time),data=deconvres)
    summary(st1)
    shapiro.test(residuals(st1))
    v <- VarCorr(st1)
    v1 <- as.data.frame(v)
    print("Subject ICC:")
    v1$vcov[1]/(v1$vcov[1]+v1$vcov[2]+v1$vcov[3]+v1$vcov[4])
}

test_celltype_mixedeffects <- function(deconvresfh) {
    deconvres <- fread(deconvresfh, header = T, sep = "\t", data.table = F)

    deconvres$Log_Adj_Unique_Molecules <- log10(deconvres$Adj_Unique_Molecules)
    deconvres$Year <- factor(deconvres$Year, levels=c("T1", "T2"))
    deconvres$Time <- factor(deconvres$Time, levels=c("Morning","Afternoon"))
    
    cdt0 <- deconvres %>% clean_names()
    cselectcelltype <- c("granulocytes","monocytes","megakaryocytes","erythrocyte_progenitors","nk_cells","b_cells","t_cells_cd4","vascular_endothelium","liver_hepatocytes")

    for (inYear in c("T1","T2")) {
        tdt <- subset(cdt0, year == inYear)
        for (testcell in cselectcelltype) {
            print(testcell)
            tmp <- summary(tdt[,testcell])
            print(tmp)
        }
    }

    for (testcell in cselectcelltype) {
        st <- lmer(paste0(testcell," ~ sex+bmi+age+(1|subject_id)+(1|year)+(1|time)"),data=cdt0)
        print(paste0(testcell," model:"))
        print(summary(st))
        print(shapiro.test(residuals(st)))
        v <- VarCorr(st)
        v1 <- as.data.frame(v)
        cicc <- v1$vcov[1]/(v1$vcov[1]+v1$vcov[2]+v1$vcov[3]+v1$vcov[4])
        print(paste0("Subject ICC: ", cicc))
    }

    for (testcell in cselectcelltype) {
        cdt0[,testcell] <- log10(cdt0[,testcell]+0.001)
        st <- lmer(paste0(testcell," ~ sex+bmi+age+(1|subject_id)+(1|year)+(1|time)"),data=cdt0)
        print(paste0(testcell," model:"))
        print(summary(st))
        print(shapiro.test(residuals(st)))
        v <- VarCorr(st)
        v1 <- as.data.frame(v)
        cicc <- v1$vcov[1]/(v1$vcov[1]+v1$vcov[2]+v1$vcov[3]+v1$vcov[4])
        print(paste0("Subject ICC: ", cicc))
    }
}
