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
#' @export
TESEA.Main <-
  function (EdgeCorScore, pathwayEdge.db, weighted.score.type = 1, 
            pathway = "kegg", gs.size.threshold.min = 15, gs.size.threshold.max = 1000, 
            reshuffling.type = "edge.labels", nperm = 100, p.val.threshold = -1, 
            FDR.threshold = 0.05, topgs = 1) 
  {
    print("Running TESEA Analysis...")
    edge.labels <- names(EdgeCorScore)
    temp <- pathwayEdge.db
    max.Ng <- length(temp)
    temp.size.G <- vector(length = max.Ng, mode = "numeric")
    for (i in 1:max.Ng) {
      temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\\t"))) - 
        2
    }
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, 
                 ncol = max.size.G)
    temp.names <- vector(length = max.Ng, mode = "character")
    temp.desc <- vector(length = max.Ng, mode = "character")
    gs.count <- 1
    for (i in 1:max.Ng) {
      edge.set.size <- length(unlist(strsplit(temp[[i]], "\\t"))) - 
        2
      gs.line <- noquote(unlist(strsplit(temp[[i]], "\\t")))
      if (pathway == "null") {
        edge.set.name <- gs.line[1]
        edge.set.desc <- gs.line[2]
        edge.set.tags <- vector(length = edge.set.size, mode = "character")
        for (j in 1:edge.set.size) {
          edge.set.tags[j] <- gs.line[j + 2]
        }
        existing.set <- is.element(edge.set.tags, edge.labels)
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < gs.size.threshold.min) || (set.size > 
                                                   gs.size.threshold.max)) 
          next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(edge.set.tags[existing.set], 
                            rep("null", max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- edge.set.name
        temp.desc[gs.count] <- edge.set.desc
        gs.count <- gs.count + 1
      }
      else if (gs.line[2] == pathway) {
        edge.set.name <- gs.line[1]
        edge.set.desc <- gs.line[2]
        edge.set.tags <- vector(length = edge.set.size, mode = "character")
        for (j in 1:edge.set.size) {
          edge.set.tags[j] <- gs.line[j + 2]
        }
        existing.set <- is.element(edge.set.tags, edge.labels)
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < gs.size.threshold.min) || (set.size > 
                                                   gs.size.threshold.max)) 
          next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(edge.set.tags[existing.set], 
                            rep("null", max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- edge.set.name
        temp.desc[gs.count] <- edge.set.desc
        gs.count <- gs.count + 1
      }
    }
    Ng <- gs.count - 1
    gs.names <- vector(length = Ng, mode = "character")
    gs.desc <- vector(length = Ng, mode = "character")
    size.G <- vector(length = Ng, mode = "numeric")
    gs.names <- temp.names[1:Ng]
    gs.desc <- temp.desc[1:Ng]
    size.G <- temp.size.G[1:Ng]
    rm(temp)
    gc()
    N <- length(EdgeCorScore)
    Obs.indicator <- matrix(nrow = Ng, ncol = N)
    Obs.RES <- matrix(nrow = Ng, ncol = N)
    Obs.ES <- vector(length = Ng, mode = "numeric")
    Obs.arg.ES <- vector(length = Ng, mode = "numeric")
    Obs.ES.norm <- vector(length = Ng, mode = "numeric")
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
    TESEA.EnrichmentScore2 <- function(edge.list, edge.set, weighted.score.type = 1, 
                                      correl.vector = NULL) {
      N <- length(edge.list)
      Nh <- length(edge.set)
      Nm <- N - Nh
      loc.vector <- vector(length = N, mode = "numeric")
      peak.res.vector <- vector(length = Nh, mode = "numeric")
      valley.res.vector <- vector(length = Nh, mode = "numeric")
      tag.correl.vector <- vector(length = Nh, mode = "numeric")
      tag.diff.vector <- vector(length = Nh, mode = "numeric")
      tag.loc.vector <- vector(length = Nh, mode = "numeric")
      loc.vector[edge.list] <- seq(1, N)
      tag.loc.vector <- loc.vector[edge.set]
      tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
      if (weighted.score.type == 0) {
        tag.correl.vector <- rep(1, Nh)
      }
      else if (weighted.score.type == 1) {
        tag.correl.vector <- correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
      }
      else if (weighted.score.type == 2) {
        tag.correl.vector <- correl.vector[tag.loc.vector] * 
          correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
      }
      else {
        tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
        tag.correl.vector <- abs(tag.correl.vector)
      }
      norm.tag <- 1/sum(tag.correl.vector)
      tag.correl.vector <- tag.correl.vector * norm.tag
      norm.no.tag <- 1/Nm
      tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
      tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 
                                                                          1)] - 1
      tag.diff.vector <- tag.diff.vector * norm.no.tag
      peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
      valley.res.vector <- peak.res.vector - tag.correl.vector
      max.ES <- max(peak.res.vector)
      min.ES <- min(valley.res.vector)
      ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), 
                   digits = 5)
      return(list(ES = ES))
    }
    obs.s2n <- vector(length = N, mode = "numeric")
    signal.strength <- vector(length = Ng, mode = "numeric")
    tag.frac <- vector(length = Ng, mode = "numeric")
    edge.frac <- vector(length = Ng, mode = "numeric")
    obs.s2n <- as.vector(EdgeCorScore)
    obs.index <- order(obs.s2n, decreasing = T)
    obs.s2n <- sort(obs.s2n, decreasing = T)
    obs.edge.labels <- edge.labels[obs.index]
    obs.edge.list2 <- obs.index
    for (i in 1:Ng) {
      edge.set <- gs[i, gs[i, ] != "null"]
      edge.set2 <- vector(length = length(edge.set), mode = "numeric")
      edge.set2 <- match(edge.set, edge.labels)
      TESEA.results <- TESEA.EnrichmentScore(edge.list = obs.edge.list2, 
                                           edge.set = edge.set2, weighted.score.type = weighted.score.type, 
                                           correl.vector = obs.s2n)
      Obs.ES[i] <- TESEA.results$ES
      Obs.arg.ES[i] <- TESEA.results$arg.ES
      Obs.RES[i, ] <- TESEA.results$RES
      Obs.indicator[i, ] <- TESEA.results$indicator
      if (Obs.ES[i] >= 0) {
        tag.frac[i] <- sum(Obs.indicator[i, 1:Obs.arg.ES[i]])/size.G[i]
        edge.frac[i] <- Obs.arg.ES[i]/N
      }
      else {
        tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
        edge.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
      }
      signal.strength[i] <- tag.frac[i] * (1 - edge.frac[i]) * 
        (N/(N - size.G[i]))
    }
    phi <- matrix(nrow = Ng, ncol = nperm)
    obs.phi <- matrix(nrow = Ng, ncol = nperm)
    for (i in 1:nperm) {
      obs.phi[, i] <- Obs.ES
    }
    if (reshuffling.type == "gene.labels") {
      for (i in 1:Ng) {
        edge.set <- gs[i, gs[i, ] != "null"]
        edge.set2 <- vector(length = length(edge.set), mode = "numeric")
        edge.set2 <- match(edge.set, edge.labels)
        for (r in 1:nperm) {
          s2n <- rnorm(length(EdgeCorScore), mean = mean(EdgeCorScore), 
                       sd = sd(EdgeCorScore))
          index <- order(s2n, decreasing = T)
          s2n <- sort(s2n, decreasing = T)
          edge.list2 <- index
          TESEA.results <- TESEA.EnrichmentScore2(edge.list = edge.list2, 
                                                edge.set = edge.set2, weighted.score.type = weighted.score.type, 
                                                correl.vector = s2n)
          phi[i, r] <- TESEA.results$ES
        }
        gc()
      }
    }
    else if (reshuffling.type == "edge.labels") {
      for (i in 1:Ng) {
        edge.set <- gs[i, gs[i, ] != "null"]
        edge.set2 <- vector(length = length(edge.set), mode = "numeric")
        edge.set2 <- match(edge.set, edge.labels)
        for (r in 1:nperm) {
          edge.list2 <- sample(1:N)
          TESEA.results <- TESEA.EnrichmentScore2(edge.list = edge.list2, 
                                                edge.set = edge.set2, weighted.score.type = weighted.score.type, 
                                                correl.vector = obs.s2n)
          stop_counter <- 1
          while (is.na(TESEA.results$ES)){
            if (stop_counter > 100) stop("Enrichment score is NA")
            edge.list2 <- sample(1:N)
            TESEA.results <- TESEA.EnrichmentScore2(edge.list = edge.list2, 
                                                  edge.set = edge.set2, weighted.score.type = weighted.score.type, 
                                                  correl.vector = obs.s2n)
            stop_counter <- stop_counter + 1
          }
          phi[i, r] <- TESEA.results$ES
        }
        gc()
      }
    }
    p.vals <- matrix(0, nrow = Ng, ncol = 2)
    for (i in 1:Ng) {
      if (Obs.ES[i] >= 0) {
        p.vals[i, 1] <- sum(phi[i, ] >= Obs.ES[i])/length(phi[i, 
                                                              ])
        p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
      }
      else {
        p.vals[i, 1] <- sum(phi[i, ] <= Obs.ES[i])/length(phi[i, 
                                                              ])
        p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
      }
    }
    p.vals[, 2] <- p.adjust(p.vals[, 1], method = "fdr")
    for (i in 1:Ng) {
      pos.phi <- Obs.ES[i]
      neg.phi <- Obs.ES[i]
      for (j in 1:nperm) {
        if (phi[i, j] >= 0) {
          pos.phi <- c(pos.phi, phi[i, j])
        }
        else {
          neg.phi <- c(neg.phi, phi[i, j])
        }
      }
      pos.m <- mean(pos.phi)
      neg.m <- mean(abs(neg.phi))
      if (Obs.ES[i] >= 0) {
        Obs.ES.norm[i] <- Obs.ES[i]/pos.m
      }
      else {
        Obs.ES.norm[i] <- Obs.ES[i]/neg.m
      }
    }
    Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
    p.val <- vector(length = Ng, mode = "numeric")
    FDR.val <- vector(length = Ng, mode = "numeric")
    Obs.ES <- signif(Obs.ES, digits = 5)
    Obs.ES.norm <- signif(Obs.ES.norm, digits = 5)
    p.val <- signif(p.vals[, 1], digits = 4)
    FDR.val <- signif(p.vals[, 2], digits = 4)
    signal.strength <- signif(signal.strength, digits = 3)
    tag.frac <- signif(tag.frac, digits = 3)
    edge.frac <- signif(edge.frac, digits = 3)
    report <- data.frame(cbind(gs.names, gs.desc, size.G, Obs.ES, 
                               Obs.ES.norm, p.val, FDR.val, tag.frac, edge.frac, signal.strength))
    names(report) <- c("GS", "SOURCE", "SIZE", "ES", "NES", "NOM p-val", 
                       "FDR q-val", "Tag \\\\%", "Edge \\\\%", "Signal")
    report2 <- report
    report.index2 <- order(Obs.ES.norm, decreasing = T)
    for (i in 1:Ng) {
      report2[i, ] <- report[report.index2[i], ]
    }
    report3 <- report
    report.index3 <- order(Obs.ES.norm, decreasing = F)
    for (i in 1:Ng) {
      report3[i, ] <- report[report.index3[i], ]
    }
    phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
    phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
    report.phen1 <- report2[1:phen1.rows, ]
    report.phen2 <- report3[1:phen2.rows, ]
    result1 <- rbind(cbind(report.phen1, Character = "+"), 
                     cbind(report.phen2, Character = "-"))
    
    n <- 1
    result2 <- list()
    for (i in 1:Ng) {
      if ((p.vals[i, 1] <= p.val.threshold) || (p.vals[i, 2] <= 
                                                FDR.threshold) || (is.element(i, c(Obs.ES.index[1:topgs], 
                                                                                   Obs.ES.index[(Ng - topgs + 1):Ng])))) {
        kk <- 1
        edge.number <- vector(length = size.G[i], mode = "character")
        edge.names <- vector(length = size.G[i], mode = "character")
        edge.list.loc <- vector(length = size.G[i], mode = "numeric")
        core.enrichment <- vector(length = size.G[i], mode = "character")
        edge.s2n <- vector(length = size.G[i], mode = "numeric")
        edge.RES <- vector(length = size.G[i], mode = "numeric")
        rank.list <- seq(1, N)
        if (Obs.ES[i] >= 0) {
          set.k <- seq(1, N, 1)
          phen.tag <- "Gain-of-correlation"
          loc <- match(i, Obs.ES.index)
        }
        else {
          set.k <- seq(N, 1, -1)
          phen.tag <- "Loss-of-correlation"
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
        filename <- paste(gs.names[i], ".", phen.tag, sep = "", 
                          collapse = "")
        result2[[n]] <- edge.report
        names(result2)[n] <- filename
        n <- n + 1
      }
    }
    result <- list(SUMMARY.RESULTS = result1, `Pathway results` = result2)
    return(result)
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
        old.mai <- par(mai = c(0.01, 0.25, 0.01, 0.3))
        on.exit(par(mai = old.mai), add = TRUE)
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