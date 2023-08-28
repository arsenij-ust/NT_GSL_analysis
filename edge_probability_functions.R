library(dplyr)

resample <- function(x, ...) x[sample.int(length(x), ...)]

#################################################
# assignSampleExpression function documentation #
#################################################

# Description:
#The assignSampleExpression function is designed to calculate and assign sample expression values to edges in an input graph using the specified sample expression data. This function can be used to augment an input graph with additional edge attributes representing aggregated sample expression values based on gene information associated with the edges.

# Arguments:
# igraph: The input igraph object representing the graph to which sample expression data should be assigned.
# sample_expr: A named numeric vector containing sample expression values for individual genes. The names of the vector should correspond to gene identifiers used in the gene_column.
# gene_column: The name of the edge attribute in the graph that contains gene identifiers. Defaults to "symbol".
# attr_name: The name of the edge attribute to which the aggregated sample expression values should be assigned. Defaults to "edge_sum".
# output_graph: A logical value indicating whether the function should return the graph with assigned edge attributes (TRUE) or just a vector of aggregated expression values (FALSE). Defaults to FALSE.

# Value: 
# The function returns either a modified igraph object with edge attributes assigned (if output_graph is TRUE), or a numeric vector containing the aggregated sample expression values for each edge (if output_graph is FALSE).

# Note:
# The function assumes that the input graph has edge attributes containing gene identifiers (default column name: "symbol"), and the sample_expr vector contains numeric expression values for those genes.
# The function uses the igraph package for working with graphs and edge attributes.

assignSampleExpression <- function(igraph, sample_expr, gene_column="symbol", attr_name = "edge_sum", output_graph=FALSE){
  
  if(output_graph){
    for(i in seq(1:length(igraph::E(igraph)))){
      genes<- unique(unlist(igraph::edge_attr(igraph, gene_column, i)))
      sum_edge_weihgt <- sum(sample_expr[genes], na.rm=TRUE)
      igraph::edge_attr(graph = igraph, name = attr_name, index = i) <- sum_edge_weihgt
    }
    return(igraph)
  } else {
    edge_vec <- c()
    for(i in seq(1:length(igraph::E(igraph)))){
      genes<- unique(unlist(igraph::edge_attr(igraph, gene_column, i)))
      sum_edge_weihgt <- sum(sample_expr[genes], na.rm=TRUE)
      edge_vec <- c(edge_vec, sum_edge_weihgt)
    }
    return(edge_vec)
  }
}

####################################################
# compute_reaction_activity function documentation #
####################################################

# Description:
# The compute_reaction_activity function is designed to compute and assign reaction activity values to edges in an input graph based on the provided gene expression counts data. This function facilitates the computation of reaction activities and the incorporation of these values as edge attributes within the graph.

# Attributes:
# igraph: The input igraph object representing the graph to which reaction activity values should be assigned.
# counts: A numeric matrix or data frame containing gene expression counts data. Each column represents a sample, and each row corresponds to a gene.
# gene_column: The name of the edge attribute in the graph that contains gene identifiers. Defaults to "symbol".
# attr_name: The name of the edge attribute to which the computed reaction activity values should be assigned. Defaults to "edge_sum".
# output_graph: A logical value indicating whether the function should return the modified graph with assigned edge attributes (TRUE) or a data frame of computed reaction activities (FALSE). Defaults to FALSE.

# Value:
# The function returns either a modified igraph object with assigned reaction activity edge attributes (if output_graph is TRUE), or a data frame containing computed reaction activity values for each sample and reaction (if output_graph is FALSE).

compute_reaction_activity <- function(igraph, counts, gene_column="symbol", attr_name = "edge_sum", output_graph=FALSE){
  
  res <- apply(counts, 2, function(x){
    assignSampleExpression(igraph, x, gene_column=gene_column, attr_name = attr_name, output_graph=output_graph)
  })
  
  if(!output_graph){
    rownames(res) <- unlist(igraph::edge_attr(igraph, "miriam.kegg.reaction"))
    res <- res[-which(duplicated(rownames(res))),]
  }
  
  return(res)
}

###########
# Note for GSL pathways:
# Reactions with more than one gene where the gene reaction rule is AND:
# R01281
###########

#######################################################
# get_edge_weight_n_steps_away function documentation #
#######################################################

# Description:
# The get_edge_weight_n_steps_away function is designed to find the weight of an edge that is a certain number of steps away from a specified node in a graph. This function is utilized within the context of the compute_transition_probability function to calculate transition probabilities along paths.

# Arguments:
# g: An adjacency matrix representing the graph. The matrix should contain edge weights as its values.
# start_col: The column index of the starting node.
# start_row: The row index of the starting node.

# Value:
# If the edge is not reachable due to no incoming edges, the function returns NULL.
# If the edge is reachable, the function returns the weight of the edge that is a certain number of steps away.

get_edge_weight_n_steps_away <- function(g, start_col, start_row) {
  current_col <- start_col
  current_row <- start_row
  
  prev_cols <- current_row
  prev_rows <- which(g[, prev_cols] != 0)
  
  if (length(prev_rows) == 0) {
    # No incoming edge, n steps away not reachable
    return(NULL)
  } else {
    # Find the maximum weight among incoming edges
    max_weight <- max(g[prev_rows, prev_cols])
    # Find the indices of the incoming edge(s) with the maximum weight
    max_indices <- which(g[prev_rows, prev_cols] == max_weight)
    # If there are multiple edges with the same maximum weight, select one randomly
    max_index <- resample(max_indices, 1)
    
    # Follow the selected incoming edge to the previous node
    current_col <- prev_cols[1]
    current_row <- prev_rows[max_index]
    
    if (g[current_row, current_col] == 1) {
      return(get_edge_weight_n_steps_away(g, current_col, current_row))
    } else {
      return(g[current_row, current_col])
    }
  }
  return(NULL)
}

###################################
# assignTP function documentation #
###################################

# Description:
# The assignTP function is designed to assign transition probabilities to edges in an input graph based on the specified transition probability data. This function aids in enhancing a graph's edges with additional attributes representing transition probabilities associated with specific reactions or entities. This function is utilized within the context of the compute_transition_probability function.

# Arguments:
# igraph: The input igraph object representing the graph to which transition probabilities should be assigned.
# transition_probability: A named numeric vector containing transition probabilities for individual reactions or entities. The names of the vector should correspond to the identifiers used in the specified column.
# column: The name of the edge attribute in the graph that contains identifiers for reactions or entities. Defaults to "miriam.kegg.reaction".
# attr_name: The name of the edge attribute to which the transition probabilities should be assigned. Defaults to "tp".


# Value:
# The function returns the input igraph object with the transition probabilities assigned as edge attributes.

# Notes:
# The function assumes that the input graph has edge attributes containing identifiers for reactions or entities (default column name: "miriam.kegg.reaction"), and the transition_probability vector contains numeric transition probability values for those identifiers.
# The function uses the igraph package for working with graphs and edge attributes.

assignTP <- function(igraph, transition_probability, column="miriam.kegg.reaction", attr_name = "tp"){
  
  for(i in seq(1:length(igraph::E(igraph)))){
    reaction<- unique(unlist(igraph::edge_attr(igraph, column, i)))
    igraph::edge_attr(graph = igraph, name = attr_name, index = i) <- transition_probability[reaction]
  }
  return(igraph)
}

#########################################################
# compute_transition_probability function documentation #
#########################################################

# Description:
# The compute_transition_probability function is designed to compute transition probabilities for edges in a list of igraph objects based on adjacency and edge attributes. This function enables the computation of transition probabilities and the option to adjust transition probability matrices for paths originating from a target node or the recursive adjustment method described in the manuscript.

# Arguments:
# igraph_list: A list of igraph objects for which transition probabilities are to be computed.
# attr_name: The name of the edge attribute in the graph that contains transition probabilities.
# attr_rownames: The name of the edge attribute that provides row names for the transition probability matrix.
# target_node: The target node from which paths are considered. 
# pass_through: A logical value indicating whether to compute transition probabilities with an additional step (recursive adjustment). Defaults to FALSE.

# Value
# The function returns a matrix of computed transition probabilities. Rows correspond to edge attributes specified by attr_rownames, and columns correspond to the igraph objects in igraph_list (samples).

compute_transition_probablity <- function(igraph_list, attr_name, attr_rownames, target_node=NULL, pass_through=FALSE){
  
  if(pass_through) target_node <- NULL
  
  # compute trasition probabilites with step length 1
  res <- sapply(igraph_list, function(g){
    adm <- as_adjacency_matrix(g, attr=attr_name, sparse = FALSE)
    adm[is.na(adm)] <- 0
    adm_norm <- t(apply(adm, 1, function(x) x/sum(x)))
    adm_norm[is.nan(adm_norm)] <- 0
    
    edge_df <- igraph::as_data_frame(g)
    edge_df$weight <- NA 
    for(r in 1:nrow(edge_df)){
      from_node <- unlist(edge_df[r,c("from")])
      to_node <- unlist(edge_df[r,c("to")])
      edge_df[r, "weight"] <- adm_norm[from_node, to_node]
    }
    return(edge_df$weight)
  })
  
  rownames(res) <- unlist(igraph::edge_attr(igraph_list[[1]], attr_rownames))
  res <- na.omit(res)
  res <- res[-which(duplicated(rownames(res))),]
  # compute paths from target 
  if(!is.null(target_node)){
    path_list <- igraph::all_simple_paths(igraph_list[[1]], from = target_node)
    prob_list <- list()
    last_reaction <- c()
    for(i in path_list){
      path_name <- paste0(names(i), collapse = "_")

      VP = names(i)
      EP = rep(VP, each=2)[-1]
      EP = EP[-length(EP)]

      eid <- get.edge.ids(mgraph, EP)
      r <- unlist(edge_attr(mgraph, attr_rownames, eid))
      last_reaction <- c(last_reaction, tail(r, n=1))
      if(length(r)>1){
        res2 <- apply(res[r,], 2, prod)
      } else {
        res2 <- res[r,]
      }
      prob_list[[path_name]] <- res2
    }

    res3 <- bind_rows(prob_list) %>% as.data.frame()
    # set rownames
    rownames(res3) <- last_reaction
    
    # fill not existing edges
    missing_edges <- setdiff(unlist(igraph::edge_attr(igraph_list[[1]], attr_rownames)), rownames(res3))
    newdf <- matrix(0, nrow=length(missing_edges), ncol = ncol(res3)) %>% as.data.frame()
    rownames(newdf) <- missing_edges
    colnames(newdf) <- colnames(res3)
    res3 <- rbind(res3, newdf)
    
    # sort by rownames 
    res <- res3[rownames(res),]
  }
  
  if(pass_through){
    g_df <- igraph::as_data_frame(igraph_list[[1]])
    r_tp_1_list <- apply(res, 2, function(x){
      graph <- assignTP(igraph_list[[1]], x)
      adm <- as_adjacency_matrix(graph, attr="tp", sparse = FALSE)
      
      t <- sapply(1:length(x), function(i){
        if(x[i] == 1){
          from_node <- unlist(g_df[which(g_df[[attr_rownames]] == names(x)[i]),"from"])
          to_node <- unlist(g_df[which(g_df[[attr_rownames]] == names(x)[i]),"to"])
          if(length(from_node)==1){
            new_weight <- tryCatch(
              get_edge_weight_n_steps_away(adm, to_node, from_node),
              error = function(e) {return(1)})
            if(!is.null(new_weight)){
              return(new_weight)
            } else {
              return(1)
            }
          } else {
            return(1)
          }
        } else {
          return(x[i])
        }
      })
    })
    rownames(r_tp_1_list) <- rownames(res)
    res <- r_tp_1_list
  }
  return(as.matrix(res))
}
