#' Distinct colors for large categorical datasets
#' 
#' @param n Number of colors desired
#' 
#' @return Color palette
#' 
#' @export
distinct_colors <- function(n) {
    x <- c("#353E7C","#007094","#009B95","#00BE7D","#96D84B","#FDE333","#040404","#3E134F","#851170","#C53270",
          "#F36E35","#F8B83C","#A23D47","#A96C00","#9A9800","#75C165","#50E2BB","#B0F4FA","#26185F","#005D9E",
          "#18BDB0","#9ADCBB","#D7F4CF","#DD008F","#DB6AC0","#D5A6DB","#F6F6FC","#CBA079","#928261","#605F4C",
          "#363932","#001889","#87008D","#DAFF47","#88002D","#FE945C","#FFE2C0","#004533","#006F69","#0091AD",
          "#EDD788","#AB4A3D","#73243C","#AC0535","#EB2C31","#EF4868","#F56553","#404E9C","#3D7CB8","#4BA5BF",
          "#55C7B1","#A3E292","#FAEF8B","#0A1230","#3C2871","#7A3392","#B7509A","#E68375","#F1C687","#985277",
          "#A37B49","#9BA453","#86CBA0","#7DEBEA","#C0FCFC","#2C2C7D","#396ABC","#5DC6DB","#AAE6EA","#DEFDFD",
          "#CA40B4","#CF82E3","#D1B6F3","#F7FFFF","#C7AEAE","#928F93","#626C7B","#3A465F","#25309C","#7B31A9",
          "#DEFD99","#7D245B","#F3A698","#FDF0F3","#215061","#367A94","#4B9CD2","#EBE4C4","#A15D70","#6B3868",
          "#9F2D65","#DB4D6A","#DF6598","#E77C8B")
    return(rep(x, floor(n/length(x))+1)[1:n])
}

#' Normalize values between a given range
#' 
#' @param x Values to normalize
#' @param a Min of range
#' @param b max of range
#' 
#' @return Normalized values
#' 
#' @export
normalize_range <- function(x, a=0, b=1) {
    (b-a)*( (x-min(x)) / (max(x)-min(x)) )+a
}

#' Colorize numerical values
#' 
#' @param x Values to normalize
#' @param resolution Limit resolution for small values
#' 
#' @return Colorized values
#' 
#' @export
colorize <- function(values, resolution=4) {
    multiplier <- 100*resolution
    colors <- rev(heat.colors(multiplier+1))
    colors[round(normalize_range(values)*multiplier, 0)+1]
}

#' Network propagation wrapper
#' 
#' @param ig An igraph object
#' @param seeds One or more seed nodes
#' @param restart Restart probability between 0-1
#' @param normalise Normalization method for the adjacency matrix
#' @param verbose Use false for quiet mode
#' @param ... Additional keyword arguments
#' 
#' @return Propagation values
#' 
#' @import dnet
#' @export
do_rwr <- function(ig, seeds, restart=0.5, normalise="laplacian", verbose=FALSE, ...) {
    mat <- matrix(0, nrow=vcount(ig), ncol=1)
    rownames(mat) <- V(ig)$name
    mat[,1] <- as.integer(rownames(mat) %in% seeds)
    rwr <- dnet::dRWR(g=ig, normalise=normalise, setSeeds=mat, restart=restart, verbose=verbose, ...)
    p <- rwr[,1]
    names(p) <- rownames(mat)
    return(p)
}

#' Extend a signature in the context of a network
#' 
#' @param ig An igraph object
#' @param signature A character vector of symbols
#' @param restart Restart probability between 0-1
#' 
#' @return An igraph and an extended signature
#' 
#' @import hypeR
#' @import igraph
#' @export
do_signature <- function(ig, signature, restart=0.5) {
    seeds <- intersect(signature, V(ig)$name)
    stopifnot(length(seeds) > 0)
    rwr <- do_rwr(ig, seeds=seeds, restart=restart)
    V(ig)$p <- rwr[match(V(ig)$name, names(rwr))]
    V(ig)$pcolor <- colorize(V(ig)$p)
    signature.ranked <- rwr[order(-rwr)]
    data = list(ig=ig, signature=signature.ranked)
    return(data)
}

#' Do enrichment for a single community
#' 
#' @param ig An igraph object
#' @param community A community label
#' @param genesets A named list of genesets
#' @param restart Restart probability between 0-1
#' @param fdr Passed to hypeR
#' @param pval Passed to hypeR
#' @param plotting Passed to hypeR
#' 
#' @return An igraph and hyp object
#' 
#' @import hypeR
#' @import igraph
#' @export
do_community <- function(ig, community, genesets, restart=0.5, fdr=1, pval=1, plotting=FALSE) {
    seeds <- V(ig)$name[V(ig)$community == community]
    rwr <- do_rwr(ig, seeds=seeds, restart=restart)
    V(ig)$p <- rwr[match(V(ig)$name, names(rwr))]
    V(ig)$pcolor <- colorize(V(ig)$p)
    signature.ranked <- names(rwr[order(-rwr)])
    hyp <- hypeR(signature.ranked, genesets, test="kstest", pval=pval, fdr=fdr, plotting=plotting)
    data <- list(ig=ig, hyp=hyp)
    return(data)
}

#' Do enrichment for all communities of the network
#' 
#' @param ig An igraph object
#' @param genesets A named list of genesets
#' @param restart Restart probability between 0-1
#' @param fdr Passed to hypeR
#' @param pval Passed to hypeR
#' @param top Limit number of genesets shown
#' 
#' @return An igraph object
#' 
#' @import igraph
#' @export
do_enrichment <- function(ig, genesets, restart=0.5, fdr=1, pval=1, top=1) {
    communities.unique <- unique(V(ig)$community)
    labels <- mapply(function(community) {
        data <- do_community(ig=ig, communit=community, genesets=genesets, restart=restart, fdr=fdr, pval=pval)
        hyp <- data$hyp
        pathways <- hyp$data$label[1:min(nrow(hyp$data), top)]
        pathways <- ifelse(is.na(pathways), "", pathways)
        pathways <- paste(pathways, collapse="\n")
        return(pathways)
    }, communities.unique, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    names(labels) <- paste("C", communities.unique, sep="@")
    V(ig)$annotation = ""
    for (label in names(labels)) {
        community <- strsplit(label, "@")[[1]][[2]]
        annotation <- labels[[label]]
        centrality <- igraph::degree(ig)[V(ig)$community == community]
        central <- names(sort(centrality, decreasing=TRUE))[1]
        V(ig)$annotation[V(ig)$name == central] <- annotation
    }
    return(ig)
}
