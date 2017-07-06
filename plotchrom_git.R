# Function for plotting moving average of a single chromosome
# for copynumber(CGHcall) object or a data frame
# User can specify the chromosome to be plotted or pass "all" for all autosomes
# mawindow is the number of bins the moving average window takes the average of
# colourlist is a list of plotting colours (not required)
# outlier is to specify below which value outliers have to be corrected (not required)
# TODO: something to show which regions are NA

plotchrom <- function(obj, chromnr, mawindow, chrlist, outlier) {
  
  ###################
  ## Libraries
  ###################
  library(RColorBrewer)
  
  ###################
  ## Check input, set variables
  ###################
  if (class(obj)=="cghCall") {
    # Get copynumber dataframe from CGHcall object
    obj <- copynumber(obj)
  }
  
  # if (missing(colourlist)) {
  #   colourlist <- brewer.pal(12, "Paired")
  # }
  
  if (missing(outlier)) {
    outlier <- -5
  }
  
  pdfwidth <- 25.6
  pdfheight <- 12.8
  filename <- paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,".pdf")
  colourlist <- brewer.pal(12, "Paired")
  counter <- 1
  mar.default <- c(5,4,4,2) + 0.1
  
  #################
  ## Other scripts
  #################
  source("/Volumes/Oscar_2TB/OLD/PDX_Merged_data/Plot_chromosome_lines.R")
  source("~/Documents/Scripts")
  
  # Remove sex chromosomes

  # Fix outliers: set all of the values below the outlier to NA to keep shape but make nicer plots
  d_outliers <- obj
  d_outliers[d_outliers < outlier] <- NA
  
  # Moving average for all chromosomes
  if (chromnr=="all") {
    if (missing(chrlist)) {
      stop("Please provide a list of chromosome numbers.")
    }
    for (x in 1:ncol(d_outliers)) {
      cat("Calculating ma for", colnames(d_outliers[x]),"\n")
      assign(paste0("ma_", mawindow,"_", x), movingAverage(test[,x], n=mawindow, centered = T))
      #assign(paste0("ma_", mawindow,"_", x), filter(d_outliers[,x], rep(1/mawindow, mawindow), sides=2))
    }
    
    # Median centering
    for (x in 1:ncol(d_outliers)) {
      cat("Median centering...\n")
      assign(paste0("med.cent_ma_", mawindow, "_", x), 
             get(paste0("ma_", mawindow,"_", x)) - median(get(paste0("ma_", mawindow,"_", x))))
      #assign(paste0("med.cent.ma", mawindow,"_", x), sweep(ma.df, 1, median(ma.df[[1]])))
    }
  
    # If filename doesn't exist create file
    if (!file.exists(filename)) {
      cat("Creating file...\n")
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_all_ma_", mawindow,".pdf"), width=pdfwidth, height=pdfheight)
    }
    # If filename already exists create file with (1) appended
    else if (!file.exists(paste0(deparse(substitute(obj)),"_chrom_all_ma_", mawindow,"(1).pdf"))) {
      cat("Creating file...\n")
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_all_ma_", mawindow,"(1).pdf"), width=pdfwidth, height=pdfheight)
    }
    # If 2 files from same command already exist put next number in parentheses
    else {
      cat("Creating file...\n")
      existing <- list.files(pattern = paste0(substr(filename, 1, nchar(filename)-4),"\\([0-9]+\\)"))
      maxnum <- max(sapply(existing, function(x) {substr(x, nchar(x) -5, nchar(x)-5)}))
      counter <- as.numeric(maxnum) + 1
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_all_ma_", mawindow,"(",counter,").pdf"), width=pdfwidth, height=pdfheight)
    }
    
    # Plot moving averages of all chromosomes
    par(mar = mar.default + c(0, 4, 0, 0))
    plot(get(paste0("med.cent_ma", mawindow,"_", 1)), type="l", lwd=0.5, ylim=c(-2,1.5), col=colourlist[1], 
         xaxt="n", xlab=paste("Chromosomes"), ylab="log2 read depth")
    # if (length(colourlist) < ncol(d_outliers)) {
    #   stop("Not enough colours in the palette. Please run again with a list of as many colours as samples.")
    # }
    for (x in 2:ncol(d_outliers)) {
      cat("Plotting line ", x,"\n")
      if (x < 2*length(colourlist)) {
        y <- 1
      }
      else {
        y <- 0
      }
      lines(get(paste0("med.cent_ma", mawindow,"_", x)), lwd=0.5, col = colourlist[(x%%length(colourlist)) + 1], 
            lty = (x %/% (length(colourlist) + y)) + 1)
    }
    Chr.names(chrlist)
    legend("bottomleft", legend = colnames(d_outliers), pch=19, col=colourlist[c(1,3:ncol(obj))])
    cat("Closing file...\n")
    dev.off()
  }
  
  # For individual chromosomes
  else {
    # Retrieve names from dataframe
    dotrows <- rownames(d_outliers)
    # List of chromosomes without features
    dotsrows <- as.numeric(gsub(":\\d*-\\d*", "", dotrows))
    # Select from d all rows with specified chromosome then calculate the moving average over all samples
    selection <- d_outliers[dotsrows==chromnr,]
    for (x in 1:ncol(selection)) {
      cat("Calculating ma for ", colnames(d_outliers[x]),"\n")
      assign(paste0("ma_", mawindow,"_", x), filter(selection[,x], rep(1/mawindow, mawindow), sides=2))
    }
    
    # Plotting prep
    filename <- paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,".pdf")
    if (!file.exists(filename)) {
      cat("Creating file...\n")
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,".pdf"), width=pdfwidth, height=pdfheight)
    }
    # If filename already exists create file with (1) appended
    else if (!file.exists(paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,"(1).pdf"))) {
      cat("Creating file...\n")
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,"(1).pdf"), width=pdfwidth, height=pdfheight)
    }
    # If 2 files from same command already exist put next number in parentheses
    else {
      cat("Creating file...\n")
      existing <- list.files(pattern = paste0(substr(filename, 1, nchar(filename)-4),"\\([0-9]+\\)"))
      maxnum <- max(sapply(existing, function(x) {substr(x, nchar(x) -5, nchar(x)-5)}))
      counter <- as.numeric(maxnum) + 1
      pdf(file=paste0(deparse(substitute(obj)),"_chrom_", chromnr,"_ma_", mawindow,"(",counter,").pdf"), width=pdfwidth, height=pdfheight)
    }
    
    # Plotting
    par(mar = mar.default + c(0, 4, 0, 0))
    plot(get(paste0("ma_", mawindow,"_", 1)), type="l", lwd=0.5, ylim=c(-2,1.5), xaxt="n", xlab=paste("Chromosome", chromnr), ylab="log2 read depth")
    # if (length(colourlist) < ncol(selection)) {
    #   stop("Not enough colours in the palette. Please run again with a list of as many colours as samples.")
    # }
    for (x in 2:ncol(selection)) {
      cat("Plotting line ", x,"\n")
      if (x < 2*length(colourlist)) {
        y <- 1
      }
      else {
        y <- 0
      }
      lines(get(paste0("ma_", mawindow,"_", x)), lwd=0.5, col = colourlist[(x%%length(colourlist)) + 1], 
            lty = (x %/% (length(colourlist) + y)) + 1)
    }
    Chr.names(chromosomes(d_outliers))
    legend("bottomleft", legend = colnames(d_outliers), pch=19, col=colourlist[c(1,3:ncol(obj))])
    cat("Closing file...\n")
    dev.off()
  }
}
