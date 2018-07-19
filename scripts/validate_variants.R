library(ggplot2)
library(reshape2)

preprocess.varscan <-  function(df.varscan){
    colnames(df.varscan)[[2]] <- "POS"
    colnames(df.varscan)[[3]] <- "REF"
    colnames(df.varscan)[[4]] <- "ALT"
    df.varscan[,"AD"] <- sapply(df.varscan[,"Cons.Cov.Reads1.Reads2.Freq.P.value"], function(x){
        x <- as.character(x)
        as.numeric(strsplit(x, ":")[[1]][4])
    })
    df.varscan[,"RAD"] <- sapply(df.varscan[,"StrandFilter.R1..R1..R2..R2..pval"], function(x){
        x <- as.character(x)
        as.numeric(strsplit(x, ":")[[1]][6])
    })
    df.varscan[,"DP"] <- sapply(df.varscan[,"Cons.Cov.Reads1.Reads2.Freq.P.value"], function(x){
        x <- as.character(x)
        as.numeric(strsplit(x, ":")[[1]][2])
    })
    ## df.varscan[,"FREQ"] <- sapply(df.varscan[,"Cons.Cov.Reads1.Reads2.Freq.P.value"], function(x){
    ##     x <- as.character(x)
    ##     y <- strsplit(x, ":")[[1]][5]
    ##     as.numeric(gsub("%", "", y))/100
    ## })
    df.varscan[,"FREQ"] <- df.varscan[,"AD"]/df.varscan[,"DP"]
    df.varscan[,"mode"] <- "varscan"
    return(df.varscan)
}

plot.snp <- function(prefix){
    ## Read in and preprocess tsv
    df.ivar <- read.table(paste("../data/variants/",prefix,".ivar.tsv", sep=""), header=TRUE, sep="\t")
    df.ivar[,"mode"] <- "ivar"
    df.varscan <- read.table(paste("../data/variants/",prefix, ".varscan.tsv", sep=""), header=TRUE, sep="\t")
    df.varscan <- preprocess.varscan(df.varscan)
    ## Remove indels from iVar output
    df.ivar.snv <- df.ivar[grep("^-|^\\+", df.ivar$ALT, value = FALSE, invert=TRUE),]
    ## Get only maximum per position because VarScan reports only 1 variant per site(assumes a diploid genome).
    df.ivar.snv.max <- aggregate(FREQ ~ POS, df.ivar.snv, max)
    df.ivar.snv <- merge(df.ivar.snv, df.ivar.snv.max)
    df.ivar.snv <- df.ivar.snv[with(df.ivar.snv, order(POS)),]
    merged.col <- merge(df.ivar.snv[,c("POS", "FREQ", "mode")], df.varscan[,c("POS", "FREQ", "mode")], by="POS", suffixes=c(".ivar", ".varscan"))
    merged.col[,"FREQ.ivar"] <- round(merged.col[,"FREQ.ivar"], 4)
    merged.col[,"FREQ.varscan"] <- round(merged.col[,"FREQ.varscan"], 4)
    rho <- cor(merged.col$FREQ.ivar, merged.col$FREQ.varscan, method="spearman")
    jit.h <- (range(merged.col$`FREQ.varscan`)[2] - range(merged.col$`FREQ.varscan`)[1])/100
    jit.w <- (range(merged.col$`FREQ.ivar`)[2] - range(merged.col$`FREQ.ivar`)[1])/100
    pdf(paste("../plots/variants/",prefix,".snp.pdf", sep = ""))
    print(ggplot(merged.col, aes(x=`FREQ.ivar`, y=`FREQ.varscan`)) + geom_abline(intercept = 0, slope = 1, color="#7e7e7e", size=1, alpha = 0.75, linetype="dashed") + geom_jitter(width = jit.w, height = jit.h) + theme_classic() + ggtitle(paste("Freq of iSNVs called by iVar vs VarScan (Spearman's Rho =", round(rho, 2), ", File =", prefix,")")))
    dev.off()
    return(rho)
}

plot.indel <- function(prefix){
    ## Read in and preprocess tsv
    df.ivar <- read.table(paste("../data/variants/",prefix,".ivar.tsv", sep=""), header=TRUE, sep="\t")
    df.ivar[,"mode"] <- "ivar"
    df.varscan <- read.table(paste("../data/variants/",prefix, ".varscan.indel.tsv", sep=""), header=TRUE, sep="\t")
    df.varscan <- preprocess.varscan(df.varscan)
    ## Remove indels from iVar output
    df.ivar.snv <- df.ivar[grep("^-|^\\+", df.ivar$ALT, value = FALSE),]
    ## Get only maximum per position because VarScan reports only 1 variant per site(assumes a diploid genome).
    df.ivar.snv.max <- aggregate(FREQ ~ POS, df.ivar.snv, max)
    df.ivar.snv <- merge(df.ivar.snv, df.ivar.snv.max)
    df.ivar.snv <- df.ivar.snv[with(df.ivar.snv, order(POS)),]
    merged.col <- merge(df.ivar.snv[,c("POS", "FREQ", "mode")], df.varscan[,c("POS", "FREQ", "mode")], by="POS", suffixes=c(".ivar", ".varscan"))
    merged.col[,"FREQ.ivar"] <- round(merged.col[,"FREQ.ivar"], 4)
    merged.col[,"FREQ.varscan"] <- round(merged.col[,"FREQ.varscan"], 4)
    rho <- cor(merged.col$FREQ.ivar, merged.col$FREQ.varscan, method="spearman")
    jit.h <- (range(merged.col$`FREQ.varscan`)[2] - range(merged.col$`FREQ.varscan`)[1])/100
    jit.w <- (range(merged.col$`FREQ.ivar`)[2] - range(merged.col$`FREQ.ivar`)[1])/100
    pdf(paste("../plots/variants/",prefix,".indel.pdf", sep = ""))
    print(ggplot(merged.col, aes(x=`FREQ.ivar`, y=`FREQ.varscan`)) + geom_abline(intercept = 0, slope = 1, color="#7e7e7e", size=1, alpha = 0.75, linetype="dashed") + geom_jitter(width = jit.w, height = jit.h) + theme_classic() + ggtitle(paste("Freq of iSNVs called by iVar vs VarScan (Spearman's Rho =", round(rho, 2), ", File =", prefix,")")))
    dev.off()
    return(rho)
}

plot.indel("ZI-merge-27_a.sorted");
plot.indel("ZI-merge-26_a.sorted");

plot.snp("BC01.trimmed.sorted");
plot.snp("BC02.trimmed.sorted");
plot.snp("ZI-merge-27_a.sorted");
plot.snp("ZI-merge-26_a.sorted");

plot.indel("ZI-merge-27_a.sorted.2");
plot.indel("ZI-merge-26_a.sorted.2");

plot.snp("BC01.trimmed.sorted.2");
plot.snp("BC02.trimmed.sorted.2");
plot.snp("ZI-merge-27_a.sorted.2");
plot.snp("ZI-merge-26_a.sorted.2");
