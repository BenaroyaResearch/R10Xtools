#' Estimate saturation of UMIs or genes based on rarefaction of reads
#' 
#' Estimate the saturation of UMI or gene detection based on rarefaction of the mapped read
#' counts from a 10X RNAseq sample. This function takes the read counts for each sample and
#' sequentially rarefies them at different levels to determine how thoroughly UMIs or genes
#' are being sampled. Optional settings include the number of intermediate points to sample
#' (default=6), the number of times to sample at each depth (default=5), and the minimum number
#' of counts for a gene to be counted as "detected" (default=1).
#' @param molecule_info the molecule_info data frame of a 10X gene-barcode matrix. Should include columns for barcode, gene, UMI, and reads.
#' @param max_reads the maximum number of reads to sample at. By default, this value is the maximum of total read counts across all barcodes, genes, and umis.
#' @param method character, either "division" or "sampling". Method "sampling" is slower but more realistic, and yields smoother curves. Method "division" is faster but more coarse and less realistic. See Details for more complete description.
#' @param depths the read depths to sample at. Either a vector of read depths at which to sample, or a single integer value giving the number of evenly-spaced depths at which to sample. 0 is always included as an additional depth for plotting facility.
#' @param nreps the number of samples to take for each library at each depth. With well-sampled libraries, 1 should be sufficient. With poorly-sampled libraries, sampling variance may be substantial, requiring higher values. Ignored if \code{method} is set to "division".
#' @param min_counts the minimum number of counts for a UMI/gene to be counted as detected. UMIs/genes with sample counts >= this value are considered detected. Defaults to 1. Set to NULL to use min_cpm.
#' @param min_cpm the minimum counts per million for a UMI/gene to be counted as detected. UMIs/genes with sample counts >= this value are considered detected. Either this or min_count should be specified, but not both; including both yields an error. Defaults to NULL.
#' @param verbose logical, whether to output the status of the estimation.
#' @import countSubsetNorm
#' @export
#' @details The \code{method} parameter determines the approach used to estimate the number of UMIs or genes detected at different read depths. Method "division" simply divides the counts for each UMI/gene by a series of scaling factors, then counts the genes whose adjusted counts exceed the detection threshold. Method "sampling" generates a number of sets (\code{nreps}) of simulated counts for each library at each sequencing depth, by probabilistically simulating counts using observed proportions. It then counts the number of genes that meet the detection threshold in each simulation, and takes the arithmetic mean of the values for each library at each depth.
#' @return A data frame containing \code{nreps} rows for each depth, with one row for each sample at each depth. Columns include "sample" (the name of the sample identifier), "depth" (the depth value for that iteration), and "sat" (the number of genes or UMIs detected at that depth for that sample).  For method "sampling", it includes an additional column with the variance of genes detected across all replicates of each sample at each depth.
#' @usage \code{
#' estimate_10X_saturation(
#'   counts, max_reads=Inf,
#'   method="sampling",
#'   depths=6, nreps=5,
#'   min_counts=1, min_cpm=NULL,
#'   verbose=FALSE)}
estimate_10X_saturation <-
  function(molecule_info, max_reads=Inf,
           method="sampling",
           depths=6, nreps=5,
           min_counts=1, min_cpm=NULL,
           verbose=FALSE) {
    if (sum(!is.null(min_counts), !is.null(min_cpm)) != 1)
      stop("One of min_counts or min_cpm must be specified, but not both.")
    method <- match.arg(method, choices=c("division", "sampling"))
    
    # add a check here to make sure the input molecule_info object is valid for what I'm doing
    
    # add functionality later for dealing with multiple samples (an aggregated molecule_info file)
    # for now, it needs to be subset to a single sample before feeding into this function
    # this would require a bunch of tweaks
    # if I want to do that, check out my estimate_saturation function for help
    molecule_info <- as.data.table(molecule_info)
    molecule_info$reads <- as.numeric(molecule_info$reads) # easier to manage if I convert to numeric
    
    
    readsums <- sum(molecule_info$reads)
    max_reads <- min(readsums, max_reads)
    n_barcodes <- length(unique(molecule_info$barcode))
    n_umis <- nrow(molecule_info)
    umi_probs <- molecule_info$reads
    
    if (length(depths) == 1) depths <- round(seq(from=0, to=max_reads, length.out=depths+1))
    saturation <-
      data.frame(depth=depths, n_barcodes=rep(n_barcodes, length(depths)))
      # data.frame(sample=as.vector(sapply(colnames(counts), FUN=rep, times=depths+1)),
      #            depth=rep(depths, time=ncol(counts)))
    
    genes.estimates <- numeric()
    umis.estimates <- numeric()
    
    if (method=="sampling") {
      genes.var.estimates <- numeric()
      umis.var.estimates <- numeric()
    }
    
    # adjust to min_cpm if specified
    if (!is.null(min_cpm)) {
      min_counts.lib <- readsums[i] / 1000000
    } else {
      min_counts.lib <- min_counts
    }
    
    for (j in depths) {
      if (verbose) cat("\nWorking on depth", j, "out of", length(depths), "different depths.\n")
      if (method=="division") {
        stop("Function not currently working with method 'division'")
        molecule_info.tmp <-
          molecule_info[
            (molecule_info$reads / max_reads * j) >= min_counts.lib,]
      } else if (method=="sampling") {
        if (all_equal(c(j,0))) {
          umis.estimates <- c(umis.estimates, 0)
          umis.var.estimates <- c(umis.var.estimates, 0)
          genes.estimates <- c(genes.estimates, 0)
          genes.var.estimates <- c(genes.var.estimates, 0)
        } else {
          umis.est <- as.numeric(rep(NA, nreps))
          genes.est <- as.numeric(rep(NA, nreps))
          
          for (k in 1:nreps) {
            if (verbose) cat("Working on rep", k, "of", nreps, "\n")
            
            ## using data.table functions
            molecule_info.tmp <-  # generate a resampled version of molecule_info
              molecule_info[
                sample.int(n_umis, size=j, replace=TRUE, prob=umi_probs),]
            umis.est[k] <- uniqueN(molecule_info.tmp[,c("barcode", "gene", "umi")])
            genes.est[k] <- uniqueN(molecule_info.tmp[,c("barcode", "gene")])
            
            ## using dplyr functions
            # molecule_info.tmp <-
            #   dplyr::sample_n(
            #     molecule_info,
            #     size=j, replace=TRUE,
            #     weight=molecule_info$reads)
            # umis.est[k] <- 
            #   molecule_info.tmp %>%
            #   dplyr::count(barcode, gene, umi) %>%
            #   nrow()
            # genes.est[k] <- 
            #   molecule_info.tmp %>%
            #   dplyr::count(barcode, gene) %>%
            #   nrow()
            
            ## using min_counts.lib
            # umis.est[k] <- sum(molecule_info.tmp$n >= min_counts.lib)
            # genes.est[k] <- sum(molecule_info.tmp$n >= min_counts.lib)
          }
          umis.estimates <- c(umis.estimates, mean(umis.est))
          umis.var.estimates <- c(umis.var.estimates, var(umis.est))
          genes.estimates <- c(genes.estimates, mean(genes.est))
          genes.var.estimates <- c(genes.var.estimates, var(genes.est))
        }
      }
    }
    
    saturation$umis <- umis.estimates
    saturation$genes <- genes.estimates
    if (exists("umis.var.estimates")) {
      saturation$umis.var <- umis.var.estimates
      saturation$genes.var <- genes.var.estimates
    }
    return(saturation)
  }
