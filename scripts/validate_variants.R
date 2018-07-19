library(ggplot2)
library(reshape2)

df.ivar <- read.table("../data/variants/BC01.trimmed.sorted.ivar.tsv", header=TRUE)
df.varscan <- read.table("../data/variants/BC01.trimmed.sorted.varscan.indel.tsv", header=TRUE)
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

df.varscan[,"FREQ"] <- sapply(df.varscan[,"Cons.Cov.Reads1.Reads2.Freq.P.value"], function(x){
    x <- as.character(x)
    y <- strsplit(x, ":")[[1]][5]
    as.numeric(gsub("%", "", y))/100
})

df.varscan[,"mode"] <- "varscan"
df.ivar[,"mode"] <- "ivar"

df.ivar.snv <- df.ivar[grep("^-|^\\+", df.ivar$ALT, value = FALSE, invert=TRUE),]

merged <- rbind(df.ivar.snv[,c("POS", "FREQ", "mode")], df.varscan[,c("POS", "FREQ", "mode")])

## ggplot(merged, aes(x=POS, y=FREQ, color=mode)) + geom_jitter(width = 100, alpha = 0.75)

## Get only maximum per position
df.ivar.snv.max <- aggregate(FREQ ~ POS, df.ivar.snv, max)
df.ivar.snv <- merge(df.ivar.snv, df.ivar.snv.max)
df.ivar.snv <- df.ivar.snv[with(df.ivar.snv, order(POS)),]

merged.col <- merge(df.ivar.snv[,c("POS", "FREQ", "mode")], df.varscan[,c("POS", "FREQ", "mode")], by="POS", suffixes=c(".ivar", ".varscan"))
merged.col[,"FREQ.ivar"] <- round(merged.col[,"FREQ.ivar"], 4)
merged.col[,"FREQ.varscan"] <- round(merged.col[,"FREQ.varscan"], 4)

apply(merged.col, 1, function(x){
    all.equal(as.numeric(x["FREQ.ivar"]), as.numeric(x["FREQ.ivar"]))
})

cor(merged.col$FREQ.ivar, merged.col$FREQ.varscan, method="spearman")
## ggplot(merged.col, aes(x=FREQ.ivar, y=FREQ.varscan)) + geom_jitter(width = 0.01)

## Indels
df.ivar.indel <- df.ivar.snv <- df.ivar[grep("^-|^\\+", df.ivar$ALT, value = FALSE),]

merged <- rbind(df.ivar.snv[,c("POS", "FREQ", "mode")], df.varscan[,c("POS", "FREQ", "mode")])

## ggplot(merged, aes(x=POS, y=FREQ, color=mode)) + geom_jitter(width = 100, alpha = 0.75)

## Get only maximum per position
df.ivar.indel.max <- aggregate(FREQ ~ POS, df.ivar.indel, max)
df.ivar.indel <- merge(df.ivar.indel, df.ivar.indel.max)
df.ivar.indel <- df.ivar.indel[with(df.ivar.indel, order(POS)),]

merged.col <- merge(df.ivar.indel[,c("POS", "FREQ", "mode", "REF", "ALT")], df.varscan[,c("POS", "FREQ", "mode", "REF", "ALT")], by="POS", suffixes=c(".ivar", ".varscan"))
merged.col[,"FREQ.ivar"] <- round(merged.col[,"FREQ.ivar"], 4)
merged.col[,"FREQ.varscan"] <- round(merged.col[,"FREQ.varscan"], 4)

apply(merged.col, 1, function(x){
    all.equal(as.numeric(x["FREQ.ivar"]), as.numeric(x["FREQ.ivar"]))
})

cor(merged.col$FREQ.ivar, merged.col$FREQ.varscan, method="spearman")
