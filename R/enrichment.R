#' Distinct colors for large categorical datasets
#' 
#' Colors derived from the R package colorspace
#' colorspace: A Toolbox for Manipulating and Assessing Colors and Palettes
#' https://cran.r-project.org/web/packages/colorspace/index.html
#' License: BSD_3_clause + https://cran.r-project.org/web/packages/colorspace/LICENSE
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
           "#DEFD99","#7D245B","#F3A698","#FDF0F3","#215061","#367A94","#4B9CD2","#EBE4C4","#A15D70","#6B3868")
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
#' @importFrom dnet dRWR
#' @importFrom igraph V vcount vertex_attr_names
#' @export
do_rwr <- function(ig, seeds, restart=0.5, normalise="laplacian", verbose=FALSE, ...) {
    stopifnot("name" %in% igraph::vertex_attr_names(ig))
    mat <- matrix(0, nrow=igraph::vcount(ig), ncol=1)
    rownames(mat) <- V(ig)$name
    mat[,1] <- as.integer(rownames(mat) %in% seeds)
    rwr <- dnet::dRWR(g=ig, normalise=normalise, setSeeds=mat, restart=restart, verbose=verbose, ...)
    p <- rwr[,1]
    names(p) <- rownames(mat)
    return(p)
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
#' @importFrom hypeR hypeR
#' @importFrom dplyr %>% arrange pull
#' @importFrom igraph V vertex_attr_names as_data_frame
#' @export
enrich_community <- function(ig, community, genesets, restart=0.5, fdr=1, pval=1, plotting=FALSE) {
    stopifnot("symbol" %in% igraph::vertex_attr_names(ig))
    stopifnot("community" %in% igraph::vertex_attr_names(ig))
    
    # Extract node identifiers
    seeds <- V(ig)$name[V(ig)$community == community]
    rwr <- do_rwr(ig, seeds=seeds, restart=restart)
    
    # Map node identifiers to propagation
    v.where <- match(V(ig)$name, names(rwr))
    
    # Copy over data
    V(ig)$p <- rwr[v.where]
    V(ig)$pcolor <- colorize(V(ig)$p)
    
    # Use symbols for the signature
    signature.ranked <-
        igraph::as_data_frame(ig, what="vertices") %>%
        dplyr::arrange(desc(p)) %>%
        dplyr::pull(symbol)
    
    # Ranked-based enrichment
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
#' 
#' @return A list of hyp objects
#' 
#' @importFrom igraph V vertex_attr_names
#' @export
enrich_communities <- function(ig, genesets, restart=0.5, fdr=1, pval=1) {
    stopifnot("symbol" %in% igraph::vertex_attr_names(ig))
    stopifnot("community" %in% igraph::vertex_attr_names(ig))

    # Enrichment for each community
    communities.unique <- unique(V(ig)$community)
    communities.hyp <- mapply(function(community) {
        data <- enrich_community(ig=ig, community=community, genesets=genesets, restart=restart, fdr=fdr, pval=pval)
        return(data$hyp)
    }, communities.unique, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    names(communities.hyp) <- paste("C", communities.unique, sep="@")
    return(communities.hyp)
}


#' Labelize community enrichment
#' 
#' @param communities.hyp A list of hyp objects
#' @param top Limit number of genesets shown
#' @param sep label delimeter
#' @param val A value to associate with genesets
#' 
#' @return A list of community enrichment labels
#' 
#' @export
labelize_communities <- function(communities.hyp, top=1, sep="\n", val=NULL) {
    communities.labels <- mapply(function(hyp) {
        if (nrow(hyp$data) == 0) return("")
        rows <- hyp$data[1:min(nrow(hyp$data), top),]
        labs <- rows$label
        if (!is.null(val)) {
            if (val %in% colnames(rows)) {
                labs <- paste(rows$label, " (", val, "=", rows[,val,drop=TRUE], ")", sep="")
            }
        }
        labs <- paste(labs, collapse=sep)
        return(labs)
    }, communities.hyp, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    return(communities.labels)
}

#' Static enrichment for all communities of the network
#' 
#' @param ig An igraph object
#' @param genesets A named list of genesets
#' @param restart Restart probability between 0-1
#' @param fdr Passed to hypeR
#' @param pval Passed to hypeR
#' @param top Limit number of genesets shown
#' @param val A value to associate with genesets
#' 
#' @return An igraph object
#' 
#' @importFrom igraph V
#' @export
enrich_communities_static <- function(ig, genesets, restart=0.5, fdr=1, pval=1, top=1, val=NULL) {
    communities.hyp <- enrich_communities(ig, genesets, restart=restart, fdr=fdr, pval=pval)
    communities.labels <- labelize_communities(communities.hyp, top=top, sep="\n", val=val)
    V(ig)$enrichment = ""
    for (lab in names(communities.labels)) {
        community <- strsplit(lab, "@")[[1]][[2]]
        enrichment <- communities.labels[[lab]]
        centrality <- igraph::degree(ig)[V(ig)$community == community]
        central <- names(sort(centrality, decreasing=TRUE))[1]
        V(ig)$enrichment[V(ig)$name == central] <- enrichment
    }
    return(ig)
}

#' Interactive enrichment for all communities of the network
#' 
#' @param ig An igraph object
#' @param genesets A named list of genesets
#' @param restart Restart probability between 0-1
#' @param fdr Passed to hypeR
#' @param pval Passed to hypeR
#' @param top Limit number of genesets shown
#' @param val A value to associate with genesets
#' 
#' @return An visnetwork object
#' 
#' @importFrom dplyr %>%
#' @importFrom visNetwork toVisNetworkData visNetwork visPhysics visEdges
#' @export
enrich_communities_interactive <- function(ig, genesets, restart=0.5, fdr=1, pval=1, top=1, val=NULL) {
    communities.hyp <- enrich_communities(ig, genesets, restart=restart, fdr=fdr, pval=pval)
    communities.labels <- labelize_communities(communities.hyp, top=top, sep="<br>", val=val)
    
    vd <- visNetwork::toVisNetworkData(ig)
    nodes <- vd$nodes
    edges <- vd$edges

    # Community labeling
    nodes$title <- sapply(nodes$community, function(x) {
        community <- paste("C", x, sep="@")
        enrichment <- communities.labels[[community]]
        return(paste(paste("<b>C", x, "</b>", sep=""), enrichment, sep="<br>"))
    })

    visNetwork(nodes, edges) %>%
    visPhysics(stabilization=FALSE) %>%
    visEdges(smooth=FALSE)
}

#' Meta interactive enrichment for all communities of the network
#' 
#' @param ig An igraph object
#' @param genesets A named list of genesets
#' @param restart Restart probability between 0-1
#' @param fdr Passed to hypeR
#' @param pval Passed to hypeR
#' @param top Limit number of genesets shown
#' @param mse Minimum shared edges threshold to visualize
#' @param val A value to associate with genesets
#' 
#' @return An visnetwork object
#' 
#' @importFrom dplyr %>% mutate select distinct
#' @importFrom tibble deframe
#' @importFrom visNetwork toVisNetworkData visNetwork visPhysics visEdges
#' @importFrom igraph V as_data_frame graph_from_adjacency_matrix
#' @export
enrich_communities_interactive_meta <- function(ig, genesets, restart=0.5, fdr=1, pval=1, top=1, mse=-inf, val=NULL) {
    communities.hyp <- enrich_communities(ig, genesets, restart=restart, fdr=fdr, pval=pval)
    communities.labels <- labelize_communities(communities.hyp, top=top, sep="\n", val=val)
    
    # Calculate edges between communities
    mat <- matrix(0, nrow=length(communities.labels), ncol=length(communities.labels))
    rownames(mat) <- colnames(mat) <- names(communities.labels)
    for (ir in rownames(mat)) {
        for (ic in colnames(mat)) {
            c1 <- as.numeric(strsplit(ir, "@")[[1]][[2]])
            c2 <- as.numeric(strsplit(ic, "@")[[1]][[2]])
            c1.v <- igraph::V(ig)[igraph::V(ig)$community == c1]
            c2.v <- igraph::V(ig)[igraph::V(ig)$community == c2]
            mat[ir,ic] <- length(E(ig)[to(c1.v) & from(c2.v)])
            mat[ic,ir] <- length(E(ig)[to(c2.v) & from(c1.v)])
        }
    }
    if (!base::isSymmetric(mat)) stop("Expected an undirected graph")
    diag(mat) <- 0

    # Remove weights below minimum shared edges threshold
    mat[mat < mse] <- 0

    ig.c <- igraph::graph_from_adjacency_matrix(mat, mode="undirected", weighted=TRUE, diag=FALSE)
    vd.c <- visNetwork::toVisNetworkData(ig.c)

    nodes <- vd.c$nodes
    edges <- vd.c$edges

    color.mapping <- igraph::as_data_frame(ig, what="vertices") %>%
        dplyr::mutate(community=paste("C", community, sep="@")) %>%
        dplyr::select(community, color) %>%
        dplyr::distinct(community, color) %>%
        tibble::deframe()
        
    nodes$color <- color.mapping[nodes$id]

    # Community labeling
    nodes$label <- sapply(nodes$id, function(x) {
        enrichment <- communities.labels[[x]]
        return(paste(gsub("@", "", x), enrichment, sep="\n"))
    })

    edges$color <- "grey"
    edges$width <- edges$weight

    visNetwork(nodes, edges) %>%
    visPhysics(stabilization=FALSE) %>%
    visEdges(smooth=FALSE, color="grey")
}

#' Extend a signature in the context of a network
#' 
#' @param ig An igraph object
#' @param signature A character vector of node identifiers
#' @param restart Restart probability between 0-1
#' 
#' @return An igrpah object
#' 
#' @export
extend_signature <- function(ig, signature, restart=0.5) {
    seeds <- V(ig)$name[match(signature, V(ig)$symbol)]
    stopifnot(length(seeds) > 0)
    rwr <- do_rwr(ig, seeds=seeds, restart=restart)
    V(ig)$p <- rwr[match(V(ig)$name, names(rwr))]
    V(ig)$pcolor <- colorize(V(ig)$p)
    return(ig)
}
