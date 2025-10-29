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

  if (is.character(fields) && length(fields) == 1) {
    # Trim whitespace
    fields <- trimws(fields)
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


#' Helper function to download all SNPs using the download API endpoint.
#'
#' @param url The API endpoint URL for downloading SNPs.
#' @param params The parameters to be sent with the request.
#'
#' @return A list containing the SNP information.
#' @keywords internal
#' @noRd 
.download_all_snps <- function(url, params) {
  params[["format"]] <- "ndjson"

  response <- httr::GET(url, query = params)
  httr::stop_for_status(response)

  content_text <- httr::content(response, "text", encoding = "UTF-8")
  lines <- strsplit(content_text, "\n")[[1]]

  snp_list <- list()
  for (line in lines) {
    if (nchar(trimws(line)) > 0) {
      snp_record <- jsonlite::fromJSON(line)
      snp_list <- c(snp_list, list(snp_record))
    }
  }

  return(snp_list)
}


#' Retrieve available list of SNP attributes.
#'
#' @return A list containing the available SNP attributes.
#'
#' @examples
#' # Retrieve all available SNP attributes
#' attributes <- snpAttributesQuery()
#' print(attributes)
#'
#' @export
snpAttributesQuery <- function() {
  url <- paste0(BASE_URL, "/snpAttributes")

  response <- httr::GET(url)
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
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes. Number of fields is limited to 20.
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#' @param pagination_from Pagination start index (default: 0)
#' @param pagination_size Pagination page size (default: 1000)
#' @param fetch_all If TRUE, retrieves all matching SNPs by downloading all pages (default: FALSE)
#'
#' @return A list containing the SNP information.
#'
#' @details
#' If fetch_all is TRUE, pagination_from and pagination_size are ignored.
#' The function will return all matching SNPs in a single list.
#' It only supports up to 1,000,000 SNPs being fetched in total.
#'
#' If using pagination (fetch_all=FALSE), you cannot fetch more than the first 10,000 SNPs over all pages.
#' pagination_from + pagination_size must be <= 10,000.
#'
#' @examples
#' # Search for SNPs on chromosome 1 between positions 10000 and 20000
#' snps <- regionQuery(
#'   chromosome_identifier = "1",
#'   start_position = 10000,
#'   end_position = 20000
#' )
#' print(snps)
#'
#' # Search with specific fields returned
#' snps <- regionQuery(
#'   chromosome_identifier = "X",
#'   start_position = 100000,
#'   end_position = 200000,
#'   fields = c("Basic Info", "chr", "pos", "ref", "alt")
#' )
#' print(snps)
#'
#' # Fetch all results
#' snps <- regionQuery(
#'   chromosome_identifier = "1",
#'   start_position = 10000,
#'   end_position = 20000,
#'   fetch_all = TRUE
#' )
#' print(snps)
#'
#' @export
regionQuery <- function(chromosome_identifier,
                        start_position = 1,
                        end_position = 100000,
                        fields = NULL,
                        filter_fields = NULL,
                        pagination_from = 0,
                        pagination_size = 1000,
                        fetch_all = FALSE) {
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

  if (fetch_all) {
    # Use the download api to fetch all results
    url <- paste0(BASE_URL, "/snp/chr/download")
    return(.download_all_snps(url, params))
  }

  if (pagination_from < 0 || pagination_size <= 0) {
    stop("pagination_from must be >= 0 and pagination_size must be > 0.")
  }

  if (pagination_from + pagination_size > 10000) {
    stop("When fetch_all is FALSE, pagination_from + pagination_size must be <= 10,000.")
  }

  url <- paste0(BASE_URL, "/snp/chr")

  params[["pagination_from"]] <- as.character(pagination_from)
  params[["pagination_size"]] <- as.character(pagination_size)

  response <- httr::GET(url, query = params)
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
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes. Number of fields is limited to 20.
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#' @param pagination_from Pagination start index (default: 0)
#' @param pagination_size Pagination page size (default: 1000)
#' @param fetch_all If TRUE, retrieves all matching SNPs by downloading all pages (default: FALSE)
#'
#' @return A list containing the SNP information.
#'
#' @details
#' If fetch_all is TRUE, pagination_from and pagination_size are ignored.
#' The function will return all matching SNPs in a single list.
#' It only supports up to 1,000,000 SNPs being fetched in total.
#'
#' If using pagination (fetch_all=FALSE), you cannot fetch more than the first 10,000 SNPs over all pages.
#' pagination_from + pagination_size must be <= 10,000.
#'
#' @examples
#' # Search for specific RSIDs
#' rsid_results <- rsidsQuery(rsid_list = c("rs123456", "rs789012"))
#' print(rsid_results)
#'
#' # Search with specific fields returned
#' rsid_results <- rsidsQuery(
#'   rsid_list = "rs123456,rs789012",
#'   fields = c("Basic Info", "chr", "pos")
#' )
#' print(rsid_results)
#'
#' # Fetch all results
#' rsid_results <- rsidsQuery(
#'   rsid_list = c("rs123456", "rs789012"),
#'   fetch_all = TRUE
#' )
#' print(rsid_results)
#'
#' @export
rsidsQuery <- function(rsid_list,
                       fields = NULL,
                       filter_fields = NULL,
                       pagination_from = 0,
                       pagination_size = 1000,
                       fetch_all = FALSE) {
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

  if (fetch_all) {
    # Use the download api to fetch all results
    url <- paste0(BASE_URL, "/snp/rsidList/download")
    return(.download_all_snps(url, params))
  }

  if (pagination_from < 0 || pagination_size <= 0) {
    stop("pagination_from must be >= 0 and pagination_size must be > 0.")
  }

  if (pagination_from + pagination_size > 10000) {
    stop("When fetch_all is FALSE, pagination_from + pagination_size must be <= 10,000.")
  }

  url <- paste0(BASE_URL, "/snp/rsidList")

  params[["pagination_from"]] <- as.character(pagination_from)
  params[["pagination_size"]] <- as.character(pagination_size)

  response <- httr::GET(url, query = params)
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
#' @param fields Fields to return, can be JSON string, file path, or vector of attributes. Number of fields is limited to 20.
#' @param filter_fields SNP attribute labels that should not be empty for the record to be retrieved
#' @param pagination_from Pagination start index (default: 0)
#' @param pagination_size Pagination page size (default: 1000)
#' @param fetch_all If TRUE, retrieves all matching SNPs by downloading all pages (default: FALSE)
#'
#' @return A list containing the SNP information.
#'
#' @details
#' If fetch_all is TRUE, pagination_from and pagination_size are ignored.
#' The function will return all matching SNPs in a single list.
#' It only supports up to 1,000,000 SNPs being fetched in total.
#'
#' If using pagination (fetch_all=FALSE), you cannot fetch more than the first 10,000 SNPs over all pages.
#' pagination_from + pagination_size must be <= 10,000.
#'
#' @examples
#' # Search for SNPs associated with a specific gene
#' gene_results <- geneQuery(gene = "BRCA1")
#' print(gene_results)
#'
#' # Search with specific fields returned
#' gene_results <- geneQuery(
#'   gene = "TP53",
#'   fields = c("Basic Info", "chr", "pos", "ref", "alt")
#' )
#' print(gene_results)
#'
#' # Fetch all results
#' gene_results <- geneQuery(gene = "BRCA1", fetch_all = TRUE)
#' print(gene_results)
#'
#' @export
geneQuery <- function(gene,
                      fields = NULL,
                      filter_fields = NULL,
                      pagination_from = 0,
                      pagination_size = 1000,
                      fetch_all = FALSE) {
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

  if (fetch_all) {
    # Use the download api to fetch all results
    url <- paste0(BASE_URL, "/snp/gene_product/download")
    return(.download_all_snps(url, params))
  }

  if (pagination_from < 0 || pagination_size <= 0) {
    stop("pagination_from must be >= 0 and pagination_size must be > 0.")
  }

  if (pagination_from + pagination_size > 10000) {
    stop("When fetch_all is FALSE, pagination_from + pagination_size must be <= 10,000.")
  }

  url <- paste0(BASE_URL, "/snp/gene_product")

  params[["pagination_from"]] <- as.character(pagination_from)
  params[["pagination_size"]] <- as.character(pagination_size)

  response <- httr::GET(url, query = params)
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
#' @examples
#' # Count SNPs on chromosome 1 between positions 10000 and 20000
#' count <- countRegionQuery(
#'   chromosome_identifier = "1",
#'   start_position = 10000,
#'   end_position = 20000
#' )
#' print(paste("Number of SNPs:", count))
#'
#' # Count with filter fields
#' count <- countRegionQuery(
#'   chromosome_identifier = "X",
#'   start_position = 100000,
#'   end_position = 200000,
#'   filter_fields = c("pos", "chr")
#' )
#' print(paste("Number of SNPs:", count))
#'
#' @export
countRegionQuery <- function(chromosome_identifier,
                             start_position = 1,
                             end_position = 100000,
                             filter_fields = NULL) {
  url <- paste0(BASE_URL, "/count/chr")

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

  response <- httr::GET(url, query = params)
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
#' @examples
#' # Count SNPs for specific RSIDs
#' count <- countRsidsQuery(rsid_list = c("rs123456", "rs789012"))
#' print(paste("Number of SNPs:", count))
#'
#' # Count with filter fields
#' count <- countRsidsQuery(
#'   rsid_list = "rs123456,rs789012",
#'   filter_fields = c("pos")
#' )
#' print(paste("Number of SNPs:", count))
#'
#' @export
countRsidsQuery <- function(rsid_list,
                            filter_fields = NULL) {
  url <- paste0(BASE_URL, "/count/rsidList")

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

  response <- httr::GET(url, query = params)
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
#' @examples
#' # Count SNPs associated with a specific gene
#' count <- countGeneQuery(gene = "BRCA1")
#' print(paste("Number of SNPs:", count))
#'
#' # Count with filter fields
#' count <- countGeneQuery(
#'   gene = "BRCA1",
#'   filter_fields = c("pos", "chr")
#' )
#' print(paste("Number of SNPs:", count))
#'
#' @export
countGeneQuery <- function(gene, filter_fields = NULL) {
  url <- paste0(BASE_URL, "/count/gene_product")

  params <- list()

  if (!is.null(gene)) {
    params[["gene"]] <- gene
  }

  if (!is.null(filter_fields)) {
    params[["filter_fields"]] <- paste(filter_fields, collapse = ",")
  }

  response <- httr::GET(url, query = params)
  httr::stop_for_status(response)

  response_content <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

  if (!"details" %in% names(response_content)) {
    stop(paste("Unexpected response from server:", jsonlite::toJSON(response_content)))
  }

  return(response_content$details)
}
