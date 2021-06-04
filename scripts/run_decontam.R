#!/usr/bin/env Rscript

library("optparse")
library("decontam")
library("reshape2")
library("stats")
#library(ggplot2); packageVersion("ggplot2")

parser <- OptionParser()
parser <- add_option(parser, c("-o", "--resout"), default="decontam_out.tsv", type="character", help="File to output results")
parser <- add_option(parser, c("-d", "--modout"), default="decontam_mod.tsv", type="character", help="File to output models")
parser <- add_option(parser, c("-i", "--counts"), default="", type="character", help="Input count table")
parser <- add_option(parser, c("-c", "--concentrations"), default="", type="character", help="Input table with DNA concentration")
parser <- add_option(parser, c("-n", "--controls"), default="", type="character", help="Input list with control sample ids")
parser <- add_option(parser, c("-m", "--method"), default="frequency", type="character", help="Method to use: frequecy, prevalence, combined")
parser <- add_option(parser, c("-t", "--threshold"), default=0.1, type="double", help="Threshold")
parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]")
args <- parse_args(parser)

generate_plot_frequency_values <- function(seqtab, taxa, conc, neg=NULL, normalize=FALSE, showModels=TRUE, log=TRUE, facet=TRUE){
  # Code copied and adapted from https://github.com/benjjneb/decontam/blob/master/R/plotting.R
  # v1.1.2 6a242fc7fc452a971b7b60b6757ea81a86ade7b5

  #print(seqtab)

  if(any(rowSums(seqtab) == 0)) { # Catch and remove zero-count samples
    zero.count <- rowSums(seqtab) == 0
    seqtab <- seqtab[!zero.count,]
    conc <- conc[!zero.count]
    if(!is.null(neg)) neg <- neg[!zero.count]
    warning("Removed ", sum(zero.count), " samples with zero total counts (or frequency).")
  }

  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  if(!(is.numeric(conc) && all(conc>0))) stop("conc must be positive numeric.")
  if(is.null(neg)) neg <- rep(FALSE, length(conc)) # Don't ignore any samples
  if(is.character(taxa)) {
    seqtab <- seqtab[,colnames(seqtab) %in% taxa,drop=FALSE]
  } else {
    stop("taxa must be a vector of taxa names.")
  }
  ntax.plot <- ncol(seqtab)
  if(ntax.plot == 0) stop("None of the provided taxa were present in seqtab.")
  # Prepare plotting data.frame
  plotdf <- cbind(data.frame(seqtab, check.names=FALSE), DNA_conc=conc, Type=ifelse(neg, "Negative", "Sample"))
  plot_melt <- melt(plotdf, measure.vars=1:ntax.plot, variable.name="taxa", value.name="taxon_abundance")
  taxon_levels <- taxa

  plot_melt$taxa <- factor(plot_melt$taxa, levels = taxon_levels)
  if(showModels) {
    mod_melts <- split(plot_melt, plot_melt$taxa)
    logc <- log(seq(min(plotdf$DNA_conc), max(plotdf$DNA_conc), length.out=1000))
    for(tax in names(mod_melts)) {
      newdata <- data.frame(logc=logc, taxa=tax, DNA_conc=exp(logc))
      freq <- mod_melts[[tax]]$taxon_abundance
      conc <- mod_melts[[tax]]$DNA_conc
      df <- data.frame(logc=log(conc), logf=log(freq))
      df <- df[!neg | is.na(neg),]
      df <- df[freq>0,]
      if(sum(freq>0)>1) {
        lm1 <- lm(logf~offset(-1*logc), data=df)
        lm0 <- lm(logf~1, data=df)
        newdata$contam <- exp(predict(lm1, newdata=newdata))
        newdata$non.contam <- exp(predict(lm0, newdata=newdata))
      } else {
        newdata$contam <- NA
        newdata$non.contam <- NA
      }
      mod_melts[[tax]] <- newdata
    }
    mod_melt <- do.call(rbind, mod_melts)
  }

  # p1 <- ggplot(data=plot_melt, aes_string("DNA_conc", "taxon_abundance")) + xlab("DNA Concentration")
  # p1 <- p1 + ylab(ifelse(normalize, "Frequency", "Relative Abundance"))
  # if(log) p1 <- p1 + scale_x_log10()
  # if(log) p1 <- p1 + scale_y_log10()
  # if(nlevels(factor(neg))>1) p1 <- p1 + aes_string(color="Type")
  # if(facet && ntax.plot > 1) p1 <- p1 + facet_wrap(~taxa)
  # if(showModels) p1 <- p1 + geom_line(data=mod_melt, aes_string(y="contam"), color="red", linetype="solid")
  # if(showModels) p1 <- p1 + geom_line(data=mod_melt, aes_string(y="non.contam"), color="black", linetype="dashed")
  # p1 + geom_point()
  # ggsave("test.png")

  # Get first and last points of the models
  idx <- sort(c(seq(1, length(mod_melt$taxa), 1000), seq(1000, length(mod_melt$taxa), 1000)))
  return(mod_melt[idx,c("contam","non.contam")])
}

# Load count table
count_table <- read.table(file=args$counts, sep='\t', header=TRUE, check.names=FALSE)
rows_table <- count_table[,1]
count_matrix <- data.matrix(data.frame(count_table[,-1], row.names = rows_table, check.names=FALSE))

# Load concentration table
if(!args$concentrations==""){
	concentrations <- read.table(file=args$concentrations, sep='\t', header=FALSE, check.names=FALSE)
	concentrations_list <- concentrations[ , "V2"]
}

# Load list of controls
if(!args$controls==""){
	controls <- read.table(file=args$controls, sep='\t', header=FALSE, check.names=FALSE)
	controls_index <- rows_table %in% controls[ , "V1"]
}

# Run DECONTAM
if (args$method=="frequency"){
	decontam_out <- isContaminant(count_matrix, normalize=FALSE, conc=concentrations_list, method="frequency", threshold=args$threshold)
}else if (args$method=="prevalence") {
	decontam_out <- isContaminant(count_matrix, normalize=FALSE, neg=controls_index, method="prevalence", threshold=args$threshold)
}else if (args$method=="combined") {
	decontam_out <- isContaminant(count_matrix, normalize=FALSE, neg=controls_index, conc=concentrations_list, method="combined", threshold=args$threshold)
}

write.table(decontam_out, file=args$resout, sep="\t", quote=FALSE)

models <- generate_plot_frequency_values(count_matrix, colnames(count_table[-1]), normalize=FALSE, conc=concentrations_list)

write.table(models, file=args$modout, sep="\t", quote=FALSE)
