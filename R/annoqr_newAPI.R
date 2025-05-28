library(httr)
library(jsonlite)

# Define the Base URL
Annotations_URL <- "https://annoq.org/api-v2/graphql"

#' create_annotations_query_string
#'
#' Convert a list of fields into a GraphQL query format.
#'
#' @param annotations A character vector of fields.
#' @return A string formatted for a GraphQL query.
#' @export
create_annotations_query_string <- function(annotations) {
  paste(annotations, collapse = "\n")
}

#' perform_graphql_query
#'
#' Perform a GraphQL query.
#'
#' @param query A character string representing the GraphQL query.
#' @return A character string with the query result.
#' @export
perform_graphql_query <- function(query) {
  response <- POST(
    Annotations_URL, 
    content_type_json(), 
    body = list(query = query),
    encode = "json"
  )
  stop_for_status(response)
  content(response, "text", encoding = "UTF-8")
}

#' regionQuery
#'
#' Get SNPs by region.
#'
#' @param chr A character string for the chromosome.
#' @param start An integer for the start position.
#' @param end An integer for the end position.
#' @param annotations_to_retrieve A character vector of fields to retrieve.
#' @return A data frame with the SNPs.
#' @export
regionQ <- function(chr, start, end, annotations_to_retrieve) {
  annotations_query_string <- create_annotations_query_string(annotations_to_retrieve)
  query <- sprintf('
  query {
    get_SNPs_by_chromosome(chr: "%s", start: %d, end: %d, query_type_option: SNPS, page_args: {size: 10000}) {
      snps {
        %s
      }
    }
  }', chr, start, end, annotations_query_string)
  
  response_content <- perform_graphql_query(query)
  data <- fromJSON(response_content, flatten = TRUE)
  data$data$get_SNPs_by_chromosome$snps
}

#' rsidQuery
#'
#' Get SNP by rsID.
#'
#' @param rsID A character string for the rsID.
#' @param annotations_to_retrieve A character vector of fields to retrieve.
#' @return A data frame with the SNPs.
#' @export
rsidQ <- function(rsID, annotations_to_retrieve) {
  annotations_query_string <- create_annotations_query_string(annotations_to_retrieve)
  query <- sprintf('
  query {
    get_SNPs_by_RsID(rsID: "%s", query_type_option: SNPS, filter_args: {exists: ["rs_dbSNP151"]}) {
      snps {
        %s
      }
    }
  }', rsID, annotations_query_string)
  
  response_content <- perform_graphql_query(query)
  data <- fromJSON(response_content, flatten = TRUE)
  data$data$get_SNPs_by_RsID$snps
}

#' rsidsQuery
#'
#' Get SNPs by multiple rsIDs.
#'
#' @param rsIDs A character vector of rsIDs.
#' @param annotations_to_retrieve A character vector of fields to retrieve.
#' @return A data frame with the SNPs.
#' @export
rsidsQ <- function(rsIDs, annotations_to_retrieve) {
  annotations_query_string <- create_annotations_query_string(annotations_to_retrieve)
  rsIDs_string <- paste(sprintf('"%s"', rsIDs), collapse = ", ")
  query <- sprintf('
  query {
    get_SNPs_by_RsIDs(rsIDs: [%s], query_type_option: SNPS, filter_args: {exists: ["rs_dbSNP151"]}, page_args: {size: 10000}) {
      snps {
        %s
      }
    }
  }', rsIDs_string, annotations_query_string)
  
  response_content <- perform_graphql_query(query)
  data <- fromJSON(response_content, flatten = TRUE)
  data$data$get_SNPs_by_RsIDs$snps
}

# Example usage:
# annotations_to_retrieve <- c("chr", "pos", "ref", "alt")
# snps <- regionQuery("18", 1, 50000, annotations_to_retrieve)
# head(snps)

# snp <- rsidQuery("rs559687999", annotations_to_retrieve)
# print(snp)

# snps <- rsidsQuery(c("rs115366554", "rs189126619"), annotations_to_retrieve)
# print(snps)
