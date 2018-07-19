library("seqinr")
library(reshape2)
library(ggplot2)

no.diff.list <- c("?", "-", "N", "n")

get.diff <- function(p){
    algn <- read.fasta(file = p)
    n <- names(algn)
    seq1 <- as.character(get(n[1], algn))
    seq2 <- as.character(get(n[2], algn))
    seq3 <- as.character(get(n[3], algn))
    l1 <- sapply(c(1:length(seq1)), function(x){
        if(seq1[x] != seq2[x]){
            if(!(seq1[x] %in% no.diff.list && seq2[x] %in% no.diff.list))
                c(seq1[x], seq2[x])
        }
    })
    l2 <- sapply(c(1:length(seq1)), function(x){
        if(seq2[x] != seq3[x]){
            if(!(seq2[x] %in% no.diff.list && seq3[x] %in% no.diff.list))
                c(seq2[x], seq3[x])
        }
    })
    l3 <- sapply(c(1:length(seq1)), function(x){
        if(seq3[x] != seq1[x]){
            if(!(seq3[x] %in% no.diff.list && seq1[x] %in% no.diff.list))
                c(seq3[x], seq1[x])
        }
    })
    l <- c(length(Filter(Negate(is.null), l1)), length(Filter(Negate(is.null), l2)), length(Filter(Negate(is.null), l3)))
    return(l)
}

files <- list.files(path="../data/alignments/", pattern="*.fasta", full.names=T, recursive=FALSE)

diff <- lapply(files, function(x){
    get.diff(x)
})

diff <- data.frame(diff)

colnames(diff) <- sapply(files, function(x){
    y <- strsplit(strsplit(x, "/")[[1]][5], "\\.")[[1]][1]
    paste(paste(unlist(strsplit(y, "_t")), collapse = " - Thres: "), "%")
});

rownames(diff) <- c("Geneious vs iVar", "iVar vs Reference", "Geneious vs Reference")

pdf("../plots/consensus.pdf")
ggplot(melt(t(diff)), aes(y=Var1, x=Var2)) + geom_tile(aes(fill = value), color="#000000") + geom_text(aes(y=Var1, x = Var2, label = value))  + scale_fill_gradient(low="#FFFFFF", high = "#FF0000") + theme_classic() + theme(axis.line <- element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Differences in Aligment") + ylab("Consensus Sequences") + ggtitle("Comparison of consensus calling using Geneious and iVar")
dev.off()
