#' Triangle creation score
#' 
#' Calculate the difference between the proportions of triangles an edge 
#' creates in the exposed group versus the control group. This is an 
#' alternative to \code{\link[ESEA]{calEdgeCorScore}} from the \code{ESEA} package 
#' to use in \code{\link{TESEA.Main}}. 
#' 
#' @param dataset A marix of gene expression data whose row names are genes 
#' symbols and whose column names are samples. 
#' @param class.labels A vector of binary labels. 
#' @param controlcharacter A character string of control sample label.
#' @param edgesbackground A marix which deposits the data of background set of 
#' edges.
#' @param return_p Should p value attributes be calculated. 
#' @param metric The distance metric to use. Defaults to "angular". 
#' 
#' @importFrom rdist cdist
#' @importFrom stats fisher.test na.omit p.adjust
#' @import dplyr
#' @export
triangle_creation_score <- function(dataset, class.labels, controlcharacter, edgesbackground, 
                                    return_p=TRUE, metric="angular"){
  combined_triangles <- triangle_creation_data(dataset=dataset, class.labels=class.labels, 
                                               controlcharacter=controlcharacter, edgesbackground=edgesbackground, 
                                               return_p=return_p, metric=metric)
  EdgeScore <- combined_triangles[, "score"]
  names(EdgeScore) <- combined_triangles[, "EdgeID"]
  if (return_p){
    attr(EdgeScore, "p_adj") <- unlist(combined_triangles[, "p_adj"], use.names = FALSE)
  }
  return(EdgeScore)
}

triangle_creation_data <- function(dataset, class.labels, controlcharacter, edgesbackground, 
                                   return_p=TRUE, metric="angular"){
  # get controls
  contr <- class.labels == controlcharacter
  # calculate cosine distance
  dists <- network_cosine_distance(dataset = dataset, edgesbackground = edgesbackground, 
                                   groups = list(control_dist = contr, exposed_dist = !contr), 
                                   metric = metric)
  # number of triangles
  control_triangles <- n_triangles_created(dists[, c("gene1", "gene2")], dists[, "control_dist"])
  exposed_triangles <- n_triangles_created(dists[, c("gene1", "gene2")], dists[, "exposed_dist"])
  # combine resutls
  combined_triangles <- dists %>% 
    dplyr::select(-dplyr::contains("location")) %>% 
    dplyr::right_join(control_triangles, by = c("gene1", "gene2")) %>%
    dplyr::full_join(exposed_triangles, by = c("gene1", "gene2"), 
                     suffix = c("_control", "_exposed"))
  combined_triangles <- combined_triangles %>%
    dplyr::mutate_(score = 'created_exposed/total_exposed - created_control/total_control', 
                   score = 'ifelse(is.na(score), 0, score)', 
                   EdgeID = 'paste(gene1, gene2, sep = "|")')
  # add test
  if (return_p){
    test_res <- triangle_test(combined_triangles)
    combined_triangles <- dplyr::bind_cols(combined_triangles, test_res)
  }
  return(combined_triangles)
}

triangle_test <- function(combined_triangles){
  combined_triangles %>% 
    dplyr::rowwise() %>% 
    dplyr::do({
      created <- rbind(.data[["created_control"]], .data[["created_exposed"]])
      total <- rbind(.data[["total_control"]], .data[["total_exposed"]])
      cont_table <- cbind(created, total - created)
      ftest <- stats::fisher.test(cont_table)
      data.frame(OR = ftest$estimate, p = ftest$p.value, 
                 lower = ftest$conf.int[1], upper = ftest$conf.int[2])
    }) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate_(p_adj = 'stats::p.adjust(p, method = "holm")')
}

network_cosine_distance <- function(dataset, edgesbackground, groups, metric="angular"){
  # generate two matched datasets
  location <- get_location(dataset, edgesbackground)
  dataset.1 <- dplyr::as_tibble(dataset[location[, "location_gene1"], ])
  dataset.2 <- dplyr::as_tibble(dataset[location[, "location_gene2"], ])
  # calculate cosine distance
  dist_fct <- function(row){
    rdist::cdist(matrix(row[seq(length(row)/2)], nrow = 1), 
                 matrix(row[seq((length(row)/2)+1, length(row))], nrow = 1), 
                 metric = metric)
  }
  dists <- lapply(groups, function(group){
    apply(dplyr::bind_cols(dataset.1[, group], dataset.2[, group]), 1, dist_fct)
  })
  return(cbind(location, as.data.frame(dists)))
}

get_location <- function(dataset, edgesbackground){
  location <- data.frame(
    location_gene1 = match(edgesbackground[, 1], rownames(dataset)), 
    location_gene2 = match(edgesbackground[, 2], rownames(dataset)))
  location <- stats::na.omit(location)
  if (nrow(location) == 0){
    stop("Not enough overlap between genes in dataset and background network")
  }
  # gene names
  location[, "gene1"] <- rownames(dataset)[location[, "location_gene1"]]
  location[, "gene2"] <- rownames(dataset)[location[, "location_gene2"]]
  return(location)
}

n_triangles_created <- function(edges, dists){
  # sort edges and dists
  dist_order <- order(dists)
  sorted_dists <- dists[dist_order]
  sorted_edges <- edges[dist_order, ]
  # neighbors
  neighbors <- list()
  triangles_created <- rep(0, nrow(sorted_edges))
  total_triangles <- rep(0, nrow(sorted_edges))
  for (i in seq(nrow(sorted_edges))){
    # add neighbors
    neighbors[[sorted_edges[i, 1]]] <- c(neighbors[[sorted_edges[i, 1]]], sorted_edges[i, 2])
    neighbors[[sorted_edges[i, 2]]] <- c(neighbors[[sorted_edges[i, 2]]], sorted_edges[i, 1])
    # check for triangles
    common_neighbors <- intersect(neighbors[[sorted_edges[i, 1]]], neighbors[[sorted_edges[i, 2]]])
    triangles_created[i] <- length(common_neighbors)
  }
  # total number of neighbors
  for (i in seq(nrow(sorted_edges))){
    # check for triangles
    common_neighbors <- intersect(neighbors[[sorted_edges[i, 1]]], neighbors[[sorted_edges[i, 2]]])
    total_triangles[i] <- length(common_neighbors)
  }
  data.frame(gene1 = sorted_edges[, 1], gene2 = sorted_edges[, 2], 
             created = triangles_created, total = total_triangles, 
             stringsAsFactors = FALSE)
}
