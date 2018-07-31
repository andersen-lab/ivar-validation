library(ggplot2)
library(reshape2)
library(patchwork)

compare.depth <- function(prefix){
    bed <- read.table("../data/test_primer.bed", sep="\t")
    title <- strsplit(prefix,"\\.")[[1]]
    colnames(bed) <- c("Regon", "Start", "End", "Name", "Qual", "Dir")
    bed["Start"] = bed["Start"] + 1
    bed["End"] = bed["End"]
    primer1 <- bed[bed$End < 1000,]
    primer2 <- bed[bed$End > 6500,]
    col.names <- c("Region", "POS", "Depth")
    depth <- read.table(paste("../data/",prefix,".select_region.depth", sep=""), sep="\t")
    colnames(depth) <- col.names
    depth.cutadapt <- read.table(paste("../data/",prefix,".select_region.cutadapt.depth", sep=""), sep="\t")
    colnames(depth.cutadapt) <- col.names
    depth.ivar <- read.table(paste("../data/",prefix,".select_region.ivar.depth", sep=""), sep="\t")
    colnames(depth.ivar) <- col.names
    depth.trimmomatic <- read.table(paste("../data/",prefix,".select_region.trimmomatic.depth", sep=""), sep="\t")
    colnames(depth.trimmomatic) <- col.names
    depth.merged <- merge(depth.cutadapt, depth.ivar, by="POS",suffixes =  c("_cutadapt", "_ivar"))
    depth.merged <- merge(depth.merged, depth, by="POS")
    depth.merged <- merge(depth.trimmomatic, depth.merged, by="POS", suffixes = c("_trimmomatic", "_untrimmed"))
    depth.m <- melt(depth.merged, id.vars <- c("POS", "Region_untrimmed", "Region_ivar", "Region_cutadapt", "Region_trimmomatic"))
    max.y1 <- max(depth.m[depth.m$POS >= primer1[,"Start"][1] & depth.m$POS <= primer1[,"End"][2],"value"])
    diff.c <- max(depth.m[depth.m$POS >= primer1[,"Start"][1] & depth.m$POS <= primer1[,"End"][1] & depth.m$variable=="Depth_cutadapt","value"])
    diff.i <- max(depth.m[depth.m$POS >= primer1[,"Start"][1] & depth.m$POS <= primer1[,"End"][1] & depth.m$variable=="Depth_ivar","value"])
    diff.t <- max(depth.m[depth.m$POS >= primer1[,"Start"][1] & depth.m$POS <= primer1[,"End"][1] & depth.m$variable=="Depth","value"])
    print(paste("Primer1%: ", (diff.c - diff.i)/diff.t))
        diff.c <- max(depth.m[depth.m$POS >= primer2[,"Start"][1] & depth.m$POS <= primer2[,"End"][1] & depth.m$variable=="Depth_cutadapt","value"])
    diff.i <- max(depth.m[depth.m$POS >= primer2[,"Start"][1] & depth.m$POS <= primer2[,"End"][1] & depth.m$variable=="Depth_ivar","value"])
    diff.t <- max(depth.m[depth.m$POS >= primer2[,"Start"][1] & depth.m$POS <= primer2[,"End"][1] & depth.m$variable=="Depth","value"])
    print(paste("Primer2%: ", (diff.c - diff.i)/diff.t))
    p1 <- ggplot(depth.m[depth.m$POS >= primer1[,"Start"][1] & depth.m$POS <= primer1[,"End"][2],], aes(x=POS, y=value, color=variable)) + geom_line(size = 0.75, alpha = 0.75) + geom_rect(data=primer1, mapping=aes(xmin=Start, xmax=End, ymin=0, ymax=max.y1), fill="#2F4F4F", alpha=0.2, inherit.aes=FALSE) + theme_classic(base_size=18) + ggtitle(paste(title,"Amplicon 1")) + scale_y_continuous(expand = c(0, 0)) + guides(fill=guide_legend(title="")) + labs(x="Position", y="Depth") + scale_colour_brewer(palette = "Set1", name="", labels=c("Trimmed by trimmomatic","Trimmed by cutadapt", "Trimmed by iVar", "Untrimmed Alignment"))
    max.y2 <- max(depth.m[depth.m$POS >= primer2[,"Start"][1] & depth.m$POS <= primer2[,"End"][2],"value"])
    p2 <- ggplot(depth.m[depth.m$POS >= primer2[,"Start"][1] & depth.m$POS <= primer2[,"End"][2],], aes(x=POS, y=value, color=variable)) + geom_line(size = 0.75, alpha = 0.75) + geom_rect(data=primer2, mapping=aes(xmin=Start, xmax=End, ymin=0, ymax=max.y2), fill="#2F4F4F", alpha=0.2, inherit.aes=FALSE) + theme_classic(base_size = 18) + ggtitle(paste(title,"Amplicon 2")) + scale_y_continuous(expand = c(0, 0)) + guides(fill=guide_legend(title="")) + labs(x="Position", y="Depth") + scale_colour_brewer(palette = "Set1", name="", labels=c("Trimmed by trimmomatic","Trimmed by cutadapt", "Trimmed by iVar", "Untrimmed Alignment"))
    plts <- list()
    plts[[1]] <- p1
    plts[[2]] <- p2
    return(plts);
}

plts <- compare.depth("ZI-merge-26_a.sorted")
p <- plts[[1]] + plts[[2]]
plts <- compare.depth("ZI-merge-27_a.sorted")
pdf("../plots/depth.pdf", width = 16, height = 15)
p + plts[[1]] + plts[[2]] + plot_layout(ncol=2)
dev.off()
