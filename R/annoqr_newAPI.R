library(httr)
library(jsonlite)

# Define the Base URL
Annotations_URL <- "http://annoq.org/api-v2/graphql"

# Convert a list of fields into a GraphQL query format
create_annotations_query_string <- function(annotations) {
  paste(annotations, collapse = "\n")
}

# Perform a GraphQL query
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

# Get SNPs by region
regionQuery <- function(chr, start, end, annotations_to_retrieve) {
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

# Get SNP by rsID
rsidQuery <- function(rsID, annotations_to_retrieve) {
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

# Get SNPs by multiple rsIDs
rsidsQuery <- function(rsIDs, annotations_to_retrieve) {
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