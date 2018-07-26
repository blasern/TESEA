#' Topological edge set enrichment analysis
#' 
#' @param EdgeCorScore A numeric vector. Each element is the differential correlation score of an edge.
#' @param pathwayEdge.db A character vector, the length of it is the number of pathways.
#' @param weighted.score.type A value. Edge enrichment correlation-based weighting: 0=no weight, 1=standard weigth, 2 = over-weigth. The default value is 1
#' @param pathway A character string of pathway database. Should be one of "kegg","reactome", "nci","huamncyc","biocarta","spike" and "panther". The default value is "kegg"
#' @param gs.size.threshold.min An integer. The minimum size (in edges) for pathways to be considered. The default value is 15. 
#' @param gs.size.threshold.max An integer. The maximum size (in edges) for pathways to be considered. The default value is 1000.
#' @param reshuffling.type A character string. The type of permutation reshuffling: "edge.labels" or "gene.labels". The default value is "edge.labels".
#' @param nperm An integer. The number of permutation reshuffling. The default value is 100.
#' @param p.val.threshold A value. The significance threshold of NOM p-value for pathways whose detail results of pathways to be presented. The  default value is -1, which means no threshold.
#' @param FDR.threshold A value. The significance threshold of FDR q-value for pathways whose detail results of pathways to be presented. The default value is 0.05.
#' @param topgs An integer. The number of top scoring gene sets used for detailed reports. The default value is 1.
#' @author Junwei Han, Xinrui Shi and Chunquan Li wrote the original in the ESEA package, small changes by Nello Blaser.
#' @importFrom stats rnorm sd
#' @export
#' @examples 
#' \dontrun{
#' # set random seed
#' set.seed(1)
#' # get data from ESEA
#' ESEA::initializeESEA()
#' edgesbackground <- ESEA::GetEdgesBackgrandData()
#' pathwayEdge.db <- ESEA::GetPathwayEdgeData()
#' dataset <- ESEA::GetExampleData("dataset")
#' class.labels <- ESEA::GetExampleData("class.labels")
#' controlcharacter <- ESEA::GetExampleData("controlcharactor") 
#' # calculate edge score (triangle version)
#' EdgeTriScore <- triangle_creation_score(dataset,
#'                                         class.labels, 
#'                                         controlcharacter, 
#'                                         edgesbackground)
#' # topological edge set enrichment analysis
#' Results_Tri <- TESEA.Main(
#'   EdgeTriScore, 
#'   pathwayEdge.db, weighted.score.type = 1,
#'   pathway = "kegg", gs.size.threshold.min = 15,
#'   gs.size.threshold.max = 1000, 
#'   reshuffling.type = "edge.labels", nperm = 1000,
#'   p.val.threshold= -1, FDR.threshold = 0.05, topgs =1)
#' # look at results
#' res_Tri <- Results_Tri[["summary"]]
#' res_Tri[res_Tri[, 'FDR q-val'] < 0.05, c("GS", "SIZE", "ES", "NES")]
#' ##                                       GS SIZE       ES       NES
#' ## 103 Maturity onset diabetes of the young   17  0.90011  2.373214
#' ## 15   Antigen processing and presentation   56 -0.54973 -1.931995
#' ## 126                p53 signaling pathway   70  0.46139  1.638150
#' ## 147               Rap1 signaling pathway  980 -0.27030 -1.296834
#' ## 148                Ras signaling pathway  976 -0.26085 -1.245997
#' }
TESEA.Main <-
  function (EdgeCorScore, pathwayEdge.db, weighted.score.type = 1, 
            pathway = "kegg", gs.size.threshold.min = 15, gs.size.threshold.max = 1000, 
            reshuffling.type = "edge.labels", nperm = 100, p.val.threshold = -1, 
            FDR.threshold = 0.05, topgs = 1) 
  {
    # initialize
    print("Running TESEA Analysis...")
    if (gs.size.threshold.min <= 1) gs.size.threshold.min <- 1
    edge.labels <- names(EdgeCorScore)
    pathway_list <- strsplit(pathwayEdge.db, "\\t")
    
    # update pathway list based on pathway and edge overlap
    pathway_list_merged <- lapply(pathway_list, function(pathway_str){
      if (pathway != "null" & pathway_str[2] != pathway){
        return(NULL)
      }
      edge.set.tags <- as.character(pathway_str[-(1:2)])
      existing.set <- is.element(edge.set.tags, edge.labels)
      return(edge.set.tags[existing.set])
    })
    pathway_size_G <- sapply(pathway_list_merged, length)
    pathway_subset <- ((pathway_size_G >= gs.size.threshold.min) &
                         (pathway_size_G <= gs.size.threshold.max)) 
    pathway_size_G <- pathway_size_G[pathway_subset]
    pathway_list_merged <- pathway_list_merged[pathway_subset]
    pathway_names <- sapply(pathway_list[pathway_subset], `[`, 1)
    pathway_descs <- sapply(pathway_list[pathway_subset], `[`, 2)
    max_size_G <- max(pathway_size_G[pathway_subset])
    n_pathways <- length(pathway_list_merged)
    Ng <- length(pathway_list_merged)
    if (length(pathway_list_merged) == 0) stop("Edge labels not matching or no pathways of specified size.")
    
    # sort edge scores 
    N <- length(EdgeCorScore)
    obs.s2n <- as.vector(EdgeCorScore)
    obs.edge.list2 <- order(obs.s2n, decreasing = T)
    obs.s2n <- obs.s2n[obs.edge.list2]
    obs.edge.labels <- edge.labels[obs.edge.list2]
    
    # Edge Enrichment Scores
    TESEA.results <- lapply(pathway_list_merged, function(edge_set){
      edge_set2 <- match(edge_set, edge.labels)
      TESEA.EnrichmentScore(edge.list = obs.edge.list2, 
                            edge.set = edge_set2, 
                            weighted.score.type = weighted.score.type, 
                            correl.vector = obs.s2n)
    })
    Obs.ES <- sapply(TESEA.results, "[[", "ES")
    Obs.arg.ES <- sapply(TESEA.results, "[[", "arg.ES")
    Obs.RES <- t(sapply(TESEA.results, "[[", "RES"))
    Obs.indicator <- t(sapply(TESEA.results, "[[", "indicator"))
    edge.frac <- ifelse(Obs.ES >= 0, Obs.arg.ES, N - Obs.arg.ES + 1)/N
    tag.frac <- mapply(function(indicator, ES, arg.ES){
      if (ES >= 0) return(sum(indicator[1:arg.ES]))
      else return(sum(indicator[arg.ES:N]))
    }, as.data.frame(t(Obs.indicator)), Obs.ES, Obs.arg.ES, 
    SIMPLIFY = TRUE, USE.NAMES = FALSE)/pathway_size_G
    signal.strength <- tag.frac * (1 - edge.frac) * (N/(N - pathway_size_G))
    
    # Permutation
    obs.phi <- matrix(rep(Obs.ES, nperm), ncol = nperm)
    phi <- t(sapply(pathway_list_merged, function(edge_set){
      edge_set2 <- match(edge_set, edge.labels)
      if (reshuffling.type == "gene.labels") {
        replicate(nperm, {
          s2n <- stats::rnorm(length(EdgeCorScore), mean = mean(EdgeCorScore), 
                              sd = stats::sd(EdgeCorScore))
          edge_list2 <- order(s2n, decreasing = T)
          s2n <- s2n[edge_list2]
          return(TESEA.EnrichmentScore2(edge.list = edge_list2, 
                                        edge.set = edge_set2, 
                                        weighted.score.type = weighted.score.type, 
                                        correl.vector = s2n))
        })
      }
      else if (reshuffling.type == "edge.labels") {
        replicate(nperm, {
          edge_list2 <- sample(1:N)
          TESEA.results <- TESEA.EnrichmentScore2(edge.list = edge_list2, 
                                                  edge.set = edge_set2, 
                                                  weighted.score.type = weighted.score.type, 
                                                  correl.vector = obs.s2n)
          stop_counter <- 1
          while (is.na(TESEA.results)){
            if (stop_counter > 100) stop("Enrichment score is NA")
            edge.list2 <- sample(1:N)
            TESEA.results <- TESEA.EnrichmentScore2(edge.list = edge_list2, 
                                                    edge.set = edge_set2, 
                                                    weighted.score.type = weighted.score.type, 
                                                    correl.vector = obs.s2n)
            stop_counter <- stop_counter + 1
          }
          return(TESEA.results)
        })
      }
    }))
    
    # p-values
    p.val <- rep(0, Ng)
    for (i in 1:Ng) {
      if (Obs.ES[i] >= 0) {
        p.val[i] <- sum(phi[i, ] >= Obs.ES[i])/length(phi[i, ])
      }
      else {
        p.val[i] <- sum(phi[i, ] <= Obs.ES[i])/length(phi[i, ])
      }
    }
    FDR.val <- p.adjust(p.val, method = "fdr")
    
    # Normalized scores
    Obs.ES.norm <- vector(length = Ng, mode = "numeric") 
    for (i in 1:Ng) {
      if (Obs.ES[i] >= 0) {
        pos.m <- mean(c(Obs.ES[i], phi[i, phi[i, ] >= 0]))
        Obs.ES.norm[i] <- Obs.ES[i]/pos.m
      }
      else {
        neg.m <- mean(abs(c(Obs.ES[i], phi[i, phi[i, ] < 0])))
        Obs.ES.norm[i] <- Obs.ES[i]/neg.m
      }
    }
    Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
    
    # Summarize results
    report <- cbind.data.frame(pathway_names, pathway_descs, pathway_size_G, Obs.ES, 
                               Obs.ES.norm, p.val, FDR.val, tag.frac, edge.frac, signal.strength)
    names(report) <- c("GS", "SOURCE", "SIZE", "ES", "NES", "NOM p-val", 
                       "FDR q-val", "Tag %", "Edge %", "Signal")
    result1 <- report[order(abs(Obs.ES.norm), decreasing = TRUE), ]
    
    # Summarize pathways
    output_pathway <- (p.val <= p.val.threshold) | (FDR.val <= FDR.threshold)
    output_pathway[c(Obs.ES.index[1:topgs], Obs.ES.index[(Ng - topgs + 1):Ng])] <- TRUE
    
    result2 <- lapply(which(output_pathway), function(i){
      kk <- 1
      edge.number <- vector(length = pathway_size_G[i], mode = "character")
      edge.names <- vector(length = pathway_size_G[i], mode = "character")
      edge.list.loc <- vector(length = pathway_size_G[i], mode = "numeric")
      core.enrichment <- vector(length = pathway_size_G[i], mode = "character")
      edge.s2n <- vector(length = pathway_size_G[i], mode = "numeric")
      edge.RES <- vector(length = pathway_size_G[i], mode = "numeric")
      rank.list <- seq(1, N)
      if (Obs.ES[i] >= 0) {
        set.k <- seq(1, N, 1)
        loc <- match(i, Obs.ES.index)
      }
      else {
        set.k <- seq(N, 1, -1)
        loc <- Ng - match(i, Obs.ES.index) + 1
      }
      for (k in set.k) {
        if (Obs.indicator[i, k] == 1) {
          edge.number[kk] <- kk
          edge.names[kk] <- obs.edge.labels[k]
          edge.list.loc[kk] <- k
          edge.s2n[kk] <- signif(obs.s2n[k], digits = 3)
          edge.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
          if (Obs.ES[i] >= 0) {
            core.enrichment[kk] <- ifelse(edge.list.loc[kk] <= 
                                            Obs.arg.ES[i], "YES", "NO")
          }
          else {
            core.enrichment[kk] <- ifelse(edge.list.loc[kk] > 
                                            Obs.arg.ES[i], "YES", "NO")
          }
          kk <- kk + 1
        }
      }
      edge.report <- data.frame(cbind(edge.number, edge.names, 
                                      edge.list.loc, edge.s2n, edge.RES, core.enrichment))
      names(edge.report) <- c("#", "EdgeID", "List Loc", 
                              "EdgeCorScore", "RES", "CORE_ENRICHMENT")
      edge.report
    })
    names(result2) <- pathway_names[output_pathway]
  
    # return result
    result <- list(summary = result1, pathways = result2)
    return(result)
  }

TESEA.EnrichmentScore <- function(edge.list, edge.set, weighted.score.type = 1, 
                                  correl.vector = NULL) {
  tag.indicator <- sign(match(edge.list, edge.set, nomatch = 0))
  no.tag.indicator <- 1 - tag.indicator
  N <- length(edge.list)
  Nh <- length(edge.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 
                                        1])
  norm.tag <- 1/sum.correl.tag
  if (norm.tag == Inf) {
    norm.tag <- 1
  }
  norm.no.tag <- 1/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - 
                  no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  }
  else {
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

TESEA.EnrichmentScore2 <- function(edge.list, edge.set, 
                                   weighted.score.type = 1, 
                                   correl.vector = NULL) {
  # set sizes
  N <- length(edge.list)
  Nh <- length(edge.set)
  Nm <- N - Nh
  
  # initialize
  loc.vector <- vector(length = N, mode = "numeric")
  loc.vector[edge.list] <- seq(1, N)
  tag.loc.vector <- sort(loc.vector[edge.set], decreasing = FALSE)
  
  # weighting
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  }
  else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
  }
  else {
    tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
  }
  tag.correl.vector <- abs(tag.correl.vector)
  
  # calculate enrichment score
  norm.tag <- 1/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1/Nm
  tag.diff.vector <- c(tag.loc.vector[1], 
                       tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)]) - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- ifelse(max.ES > -min.ES, max.ES, min.ES)
  return(ES)
}

#' Plot the pathway-result network diagram
#' 
#' @param graph A dataframe of pathway result obtained from the TESEA.main function.
#' @param margin A numeric. The value is usually between -0.5 and 0.5, which is able to zoom in or out a pathway graph. The default is 0.
#' @param vertex.label.cex 	A numeric vector of node label size.
#' @param vertex.label.font A numeric vector of label font.
#' @param vertex.size A numeric vector of Node size. See plot.igraph
#' @param vertex.size2 A numeric vector of Node size.
#' @param vertex.shape A vector of node shape. The default is graphics_type.
#' @param layout A matrix of x-y coordinates with two dims. Determine the placement of the nodes for drawing a graph.The default is layout.random.
#' @param vertex.label.color A vector of node label colors. The default is black.
#' @param vertex.color A vector of node colors. The default is the KEGG node color.
#' @param vertex.frame.color A vector of node frame color. The default is dimgray.
#' @param edge.color A vector of edge color. The default is dimgray.
#' @param axes A logical. whether to plot axes. The default is FALSE.
#' @param xlab A character string. The label of the horizontal axis. The default is the empty string.
#' @param ylab A character string. The label of the vertical axis. The default is the empty string.
#' @param sub A character string of subtitle.
#' @param main A character string of main title.
#' @param ... The arguments passed to or from methods. See plot.igraph and see plot.
#' @author Junwei Han, Xinrui Shi and Chunquan Li wrote the original in the ESEA package, small changes by Nello Blaser.
#' @importFrom graphics par
#' @export
PlotTESEAPathwayGraph <- 
  function (graph, margin = 0, vertex.label.cex = 0.6, vertex.label.font = 1, 
            vertex.size = 8, vertex.size2 = 6, vertex.shape = "rectangle", 
            layout = igraph::layout.random, vertex.label.color = "black", edge.color = "dimgray", 
            vertex.color = "#C1FFC1", vertex.frame.color = "dimgray", 
            axes = FALSE, xlab = "", ylab = "", sub = NULL, main = NULL, 
            ...) 
  {
    if (!requireNamespace("igraph", quietly = TRUE)){
      stop('Package "igraph" needed for "PlotTESEAPathwayGraph" to work. Please install it.')
    }
    
    Split <- unlist(strsplit(as.character(graph[, 2]), split = "|", 
                             fixed = T))
    graph <- igraph::graph.data.frame(
      as.data.frame(cbind(first = Split[(1:(length(Split)/2)) * 2 - 1], 
                          sencond = Split[(1:(length(Split)/2)) * 2], 
                          graph[, c(4, 6)]), stringsAsFactors = FALSE), directed = F)
    if (class(graph) != "igraph") 
      stop("the graph should be a igraph graph.")
    if (igraph::vcount(graph) == 0) {
      print("the graph is an empty graph.")
    }
    else {
      igraph::E(graph)$color <- ifelse(igraph::E(graph)$CORE_ENRICHMENT == 
                                         "YES", "red", "dimgray")
      if (length(vertex.shape) == 0) 
        vertex.shape <- NULL
      if (length(vertex.color) == 0) 
        vertex.color <- NULL
      if (length(layout) == 0) 
        layout <- NULL
      if ((axes == FALSE) && xlab == "" && ylab == "" && is.null(sub) && 
          is.null(main)) {
        old.mai <- graphics::par(mai = c(0.01, 0.25, 0.01, 0.3))
        on.exit(graphics::par(mai = old.mai), add = TRUE)
      }
      igraph::plot.igraph(graph, margin = margin, vertex.label.cex = vertex.label.cex, 
                          vertex.label.font = vertex.label.font, vertex.size = vertex.size, 
                          vertex.size2 = vertex.size2, vertex.shape = vertex.shape, 
                          layout = layout, vertex.label.color = vertex.label.color, 
                          edge.color = igraph::E(graph)$color, vertex.color = vertex.color, 
                          vertex.frame.color = vertex.frame.color, 
                          axes = axes, xlab = xlab, ylab = ylab, sub = sub, main = main, 
                          ...)
    }
  }
