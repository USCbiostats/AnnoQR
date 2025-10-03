# Annoq API R Client
#
# An R package for accessing SNP data from Annoq.org

# Load required libraries
if (!requireNamespace("httr", quietly = TRUE)) {
  stop("Package 'httr' is required but not installed.")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required but not installed.")
}

# Base URL for the Annoq API
BASE_URL <- "https://api-v2-dev.annoq.org"


# Process the fields parameter to handle the three possible input types:
# 1. JSON string: {"_source":["Basic Info","chr","pos","ref","alt","rs_dbSNP151"]}
# 2. File path: path to a file containing the JSON config
# 3. List of attributes: ["Basic Info", "chr", "pos", "ref", "alt", "rs_dbSNP151"]
#
# Returns the JSON string representation or NULL if fields is NULL.
.process_fields_param <- function(fields) {
  if (is.null(fields)) {
    return(NULL)
  }
  
  if (is.character(fields)) {
    # Check if it's a file path by attempting to read it
    if (startsWith(fields, "{") && endsWith(fields, "}")) {
      # It's a JSON string
      return(fields)
    } else {
      # It might be a file path, try to read it
      if (file.exists(fields)) {
        content <- readChar(fields, file.info(fields)$size)
        return(content)
      } else {
        # If it's not a valid file path, treat it as a JSON string (though invalid)
        stop(paste("Fields parameter appears to be a file path but file not found:", fields))
      }
    }
  } else if (is.vector(fields) && !is.null(names(fields))) {
    # This is not a list, so it's probably a named vector, not what we want
    stop(paste("Fields parameter must be a character string (JSON or file path), vector of attributes, or NULL. Got:", class(fields)))
  } else if (is.vector(fields) || is.list(fields)) {
    # Convert vector/list to the required JSON format
    json_fields <- jsonlite::toJSON(list("_source" = fields), auto_unbox = TRUE)
    return(jsonlite::toJSON(jsonlite::fromJSON(json_fields), auto_unbox = TRUE))
  } else {
    stop(paste("Fields parameter must be a character string (JSON or file path), vector of attributes, or NULL. Got:", class(fields)))
  }
}


#' Retrieve available list of SNP attributes.
#'
#' @return A list containing the available SNP attributes.
#'
#' @export
snpAttributesQuery <- function() {
  url <- paste0(BASE_URL, "/fastapi/snpAttributes")
  
  response <- httr::POST(url)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"results" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$results)
}


#' Search for SNPs by chromosome id and position range.
#'
#' @param chromosome_identifier Chromosome id to search (e.g., "1", "2", "X")
#' @param start_position Start position region of search (default: 1)
#' @param end_position End position region of search (default: 100000)
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return A list containing the SNP information.
#'
#' @export
regionQuery <- function(chromosome_identifier,
                            start_position = NULL,
                            end_position = NULL,
                            fields = NULL,
                            filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/snp/chr")
  
  params <- list("chromosome_identifier" = chromosome_identifier)
  
  if (!is.null(start_position)) {
    params[["start_position"]] <- as.character(start_position)
  }
  if (!is.null(end_position)) {
    params[["end_position"]] <- as.character(end_position)
  }
  
  processed_fields <- .process_fields_param(fields)
  if (!is.null(processed_fields)) {
    params[["fields"]] <- processed_fields
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  # Note: pagination parameters are ignored as they don't function
  # But they are still required by the API (dummy values used)
  params[["pagination_from"]] <- "0"
  params[["pagination_size"]] <- "10"
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}


#' Search for specified list of RSIDs.
#'
#' @param rsid_list List of RSIDs to search, can be comma-separated string or vector of strings
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return A list containing the SNP information.
#'
#' @export
rsidsQuery <- function(rsid_list = NULL,
                                  fields = NULL,
                                  filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/snp/rsidList")
  
  params <- list()
  
  if (!is.null(rsid_list)) {
    if (is.vector(rsid_list)) {
      params[["rsid_list"]] <- paste(rsid_list, collapse = ",")
    } else {
      params[["rsid_list"]] <- rsid_list
    }
  }
  
  processed_fields <- .process_fields_param(fields)
  if (!is.null(processed_fields)) {
    params[["fields"]] <- processed_fields
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  # Note: pagination parameters are ignored as they don't function
  # But they are still required by the API (dummy values used)
  params[["pagination_from"]] <- "0"
  params[["pagination_size"]] <- "100"
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}


#' Search for specified gene product; this can be a gene id, gene symbol or UniProt id.
#'
#' @param gene Gene product to search
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return A list containing the SNP information.
#'
#' @export
geneQuery <- function(gene = NULL,
                                     fields = NULL,
                                     filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/snp/gene_product")
  
  params <- list()
  
  if (!is.null(gene)) {
    params[["gene"]] <- gene
  }
  
  processed_fields <- .process_fields_param(fields)
  if (!is.null(processed_fields)) {
    params[["fields"]] <- processed_fields
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  # Note: pagination parameters are ignored as they don't function
  # But they are still required by the API (dummy values used)
  params[["pagination_from"]] <- "0"
  params[["pagination_size"]] <- "100"
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}


#' Count SNPs based on specified chromosome, start position, end position and filter arguments.
#'
#' @param chromosome_identifier The chromosome number (or 'X' for the X-chromosome)
#' @param start_position Start position region of search (default: 1)
#' @param end_position End position region of search (default: 100000)
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return The count of SNPs matching the criteria.
#'
#' @export
countRegionQuery <- function(chromosome_identifier,
                              start_position = NULL,
                              end_position = NULL,
                              filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/count/chr")
  
  params <- list("chromosome_identifier" = chromosome_identifier)
  
  if (!is.null(start_position)) {
    params[["start_position"]] <- as.character(start_position)
  }
  if (!is.null(end_position)) {
    params[["end_position"]] <- as.character(end_position)
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}


#' Count the number of SNPs defined in the system that have matching RSIDs from the specified list.
#'
#' @param rsid_list List of RSIDs to search, can be comma-separated string or vector of strings
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return The count of SNPs matching the criteria.
#'
#' @export
countRsidsQuery <- function(rsid_list = NULL,
                                    filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/count/rsidList")
  
  params <- list()
  
  if (!is.null(rsid_list)) {
    if (is.vector(rsid_list)) {
      params[["rsid_list"]] <- paste(rsid_list, collapse = ",")
    } else {
      params[["rsid_list"]] <- rsid_list
    }
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}


#' Count the number of SNPs defined in the system that have been associated for the specified gene product.
#'
#' @param gene Gene product to search (gene id, gene symbol or UniProt id)
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#'
#' @return The count of SNPs matching the criteria.
#'
#' @export
countGeneQuery <- function(gene = NULL, filter_fields = NULL) {
  url <- paste0(BASE_URL, "/fastapi/count/gene_product")
  
  params <- list()
  
  if (!is.null(gene)) {
    params[["gene"]] <- gene
  }
  
  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }
  
  response <- httr::POST(url, query = params)
  httr::stop_for_status(response)
  
  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }
  
  return(response_content$details)
}
