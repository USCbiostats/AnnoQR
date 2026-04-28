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
BASE_URL <- "https://api-v2.annoq.org"
SNPWAY_BASE_URL_DEFAULT <- "https://enrichment-dev.annoq.org"


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

  resp <- httr::POST(url, query = params)
  httr::stop_for_status(resp)

  txt <- httr::content(resp, as = "text", encoding = "UTF-8")

  con <- textConnection(txt)
  on.exit(close(con), add = TRUE)

  # Parse NDJSON (one JSON object per line) into a data.frame
  df <- jsonlite::stream_in(con, verbose = FALSE)

  # Ensure it's a base data.frame (not tibble), just in case
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  df
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
    params[["start_position"]] <- sprintf("%d", start_position)
  }
  if (!is.null(end_position)) {
    params[["end_position"]] <- sprintf("%d", end_position)
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

  params[["pagination_from"]] <- sprintf("%d", pagination_from)
  params[["pagination_size"]] <- sprintf("%d", pagination_size)

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

  params[["pagination_from"]] <- sprintf("%d", pagination_from)
  params[["pagination_size"]] <- sprintf("%d", pagination_size)

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

  params[["pagination_from"]] <- sprintf("%d", pagination_from)
  params[["pagination_size"]] <- sprintf("%d", pagination_size)

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
    params[["start_position"]] <- sprintf("%d", start_position)
  }
  if (!is.null(end_position)) {
    params[["end_position"]] <- sprintf("%d", end_position)
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


.resolve_snpway_base_url <- function(base_url = NULL) {
  if (!is.null(base_url) && nzchar(trimws(base_url))) {
    return(sub("/+$", "", trimws(base_url)))
  }

  env_value <- trimws(Sys.getenv("ANNOQ_SNPWAY_BASE_URL", unset = ""))
  if (nzchar(env_value)) {
    return(sub("/+$", "", env_value))
  }

  return(sub("/+$", "", SNPWAY_BASE_URL_DEFAULT))
}


.normalize_chr_label <- function(raw_chr) {
  normalized <- trimws(as.character(raw_chr))
  if (!nzchar(normalized)) {
    return("")
  }
  return(gsub("^chr", "", normalized, ignore.case = TRUE))
}


.unique_preserve_order <- function(values) {
  values[!duplicated(values)]
}


.parse_vcf_to_ids <- function(vcf_text) {
  ids <- c()
  lines <- unlist(strsplit(as.character(vcf_text), "\n", fixed = TRUE))

  for (line in lines) {
    trimmed <- trimws(line)
    if (!nzchar(trimmed) || startsWith(trimmed, "#")) {
      next
    }

    fields <- strsplit(trimmed, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) {
      fields <- strsplit(trimmed, "\\s+")[[1]]
    }

    if (length(fields) < 2) {
      next
    }

    chrom <- .normalize_chr_label(fields[[1]])
    pos <- trimws(as.character(fields[[2]]))

    if (!nzchar(chrom) || !grepl("^[0-9]+$", pos)) {
      next
    }

    ids <- c(ids, paste0(chrom, ":", pos))
  }

  .unique_preserve_order(ids)
}


.normalize_rsid_list <- function(rsid_list) {
  if (is.character(rsid_list) && length(rsid_list) == 1) {
    parsed <- unlist(strsplit(rsid_list, "[[:space:],]+"))
    parsed <- trimws(parsed)
    parsed <- parsed[nzchar(parsed)]
    return(.unique_preserve_order(parsed))
  }

  if (is.vector(rsid_list) || is.list(rsid_list)) {
    parsed <- trimws(as.character(unlist(rsid_list)))
    parsed <- parsed[nzchar(parsed)]
    return(.unique_preserve_order(parsed))
  }

  character(0)
}


.build_snpway_payload <- function(vcf_text = NULL,
                                  chrom_pos_ids = NULL,
                                  chromosome_identifier = NULL,
                                  start_position = NULL,
                                  end_position = NULL,
                                  rsid_list = NULL) {
  has_vcf <- !is.null(vcf_text) && nzchar(trimws(as.character(vcf_text)))
  has_chrom_pos <- !is.null(chrom_pos_ids) && length(chrom_pos_ids) > 0
  has_region <- !is.null(chromosome_identifier) || !is.null(start_position) || !is.null(end_position)
  has_rsid <- !is.null(rsid_list)

  modes_used <- sum(c(has_vcf, has_chrom_pos, has_region, has_rsid))

  if (modes_used != 1) {
    stop(
      "Provide exactly one input mode: vcf_text, chrom_pos_ids, chromosome range, or rsid_list.",
      call. = FALSE
    )
  }

  if (has_vcf) {
    ids <- .parse_vcf_to_ids(vcf_text)
    if (length(ids) == 0) {
      stop("No valid CHROM/POS entries were found in vcf_text.", call. = FALSE)
    }

    return(list(
      input_type = "ids",
      idsQuery = list(ids = unname(as.list(ids)))
    ))
  }

  if (has_chrom_pos) {
    ids <- trimws(as.character(unlist(chrom_pos_ids)))
    ids <- ids[nzchar(ids)]
    ids <- .unique_preserve_order(ids)

    if (length(ids) == 0) {
      stop("chrom_pos_ids must contain at least one non-empty ID.", call. = FALSE)
    }

    return(list(
      input_type = "ids",
      idsQuery = list(ids = unname(as.list(ids)))
    ))
  }

  if (has_region) {
    if (is.null(chromosome_identifier) || is.null(start_position) || is.null(end_position)) {
      stop(
        "chromosome_identifier, start_position, and end_position are required for chromosome range input mode.",
        call. = FALSE
      )
    }

    normalized_chr <- .normalize_chr_label(chromosome_identifier)
    if (!nzchar(normalized_chr)) {
      stop("chromosome_identifier cannot be empty.", call. = FALSE)
    }

    return(list(
      input_type = "chromosome",
      chrQuery = list(
        chr = normalized_chr,
        start = as.integer(start_position),
        end = as.integer(end_position)
      )
    ))
  }

  normalized_rsids <- .normalize_rsid_list(rsid_list)
  if (length(normalized_rsids) == 0) {
    stop("rsid_list must contain at least one rsID.", call. = FALSE)
  }

  list(
    input_type = "rsIdList",
    rsIdListQuery = list(rsIdList = unname(as.list(normalized_rsids)))
  )
}


.parse_snpway_error <- function(response) {
  text_body <- httr::content(response, as = "text", encoding = "UTF-8")

  if (nzchar(text_body)) {
    parsed <- tryCatch(
      jsonlite::fromJSON(text_body),
      error = function(e) NULL
    )

    if (is.list(parsed) && !is.null(parsed$detail)) {
      return(as.character(parsed$detail))
    }

    return(text_body)
  }

  paste("Request failed with status", httr::status_code(response))
}


.post_snpway_request <- function(endpoint,
                                 payload,
                                 base_url = NULL,
                                 timeout_seconds = 120) {
  url <- paste0(.resolve_snpway_base_url(base_url), endpoint)

  response <- httr::POST(
    url,
    body = payload,
    encode = "json",
    httr::timeout(timeout_seconds)
  )

  if (httr::status_code(response) >= 400) {
    stop(.parse_snpway_error(response), call. = FALSE)
  }

  response_text <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed <- jsonlite::fromJSON(response_text, simplifyVector = FALSE)

  if (!is.list(parsed)) {
    stop("SNPWay API returned an unexpected response shape.", call. = FALSE)
  }

  parsed
}

.get_relevant_columns <- function(annotation_dataset = NULL) {
  base_columns <- c("rsId", "PANTHER_ID", "mappedGenes")

  all_columns <- c(
    base_columns,
    "PANTHER_family",
    "PANTHER_Subfamily",
    "PANTHER_Pathway",
    "Protein_Class",
    "Reactome_Pathway",
    "GO_database_MF_complete",
    "GO_database_BP_complete",
    "GO_database_CC_complete",
    "PANTHER_GO_slim_Molecular_Function",
    "PANTHER_GO_slim_Biological_Process",
    "PANTHER_GO_slim_Cellular_Component"
  )

  if (is.null(annotation_dataset)) {
    return(all_columns)
  }

  dataset_column_map <- list(
    "GO:0008150" = c(base_columns, "GO_database_BP_complete"),
    "GO:0003674" = c(base_columns, "GO_database_MF_complete"),
    "GO:0005575" = c(base_columns, "GO_database_CC_complete"),
    "ANNOT_TYPE_ID_PANTHER_PATHWAY" = c(base_columns, "PANTHER_Pathway"),
    "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF" = c(base_columns, "PANTHER_GO_slim_Molecular_Function"),
    "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" = c(base_columns, "PANTHER_GO_slim_Biological_Process"),
    "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC" = c(base_columns, "PANTHER_GO_slim_Cellular_Component"),
    "ANNOT_TYPE_ID_PANTHER_PC" = c(base_columns, "Protein_Class"),
    "ANNOT_TYPE_ID_REACTOME_PATHWAY" = c(base_columns, "Reactome_Pathway")
  )

  selected <- dataset_column_map[[annotation_dataset]]
  if (is.null(selected)) {
    return(all_columns)
  }

  selected
}

.create_results_table_data <- function(response_data,
                                       panther_ids_to_include = NULL,
                                       annotation_dataset = NULL) {
  rs_id_genes_map <- response_data$rsId_genes_map
  if (is.null(rs_id_genes_map) || !is.list(rs_id_genes_map)) {
    rs_id_genes_map <- list()
  }

  panther_gene_info <- response_data$panther_gene_info
  if (is.null(panther_gene_info) || !is.list(panther_gene_info)) {
    panther_gene_info <- list()
  }

  gene_panther_mapping <- response_data$gene_panther_mapping
  if (is.null(gene_panther_mapping) || !is.list(gene_panther_mapping)) {
    gene_panther_mapping <- list()
  }

  panther_ids_set <- if (is.null(panther_ids_to_include)) NULL else unique(as.character(panther_ids_to_include))
  table_data <- list()
  processed_pairs <- new.env(parent = emptyenv())

  for (rs_id in names(rs_id_genes_map)) {
    genes_for_rs_id <- rs_id_genes_map[[rs_id]]
    if (is.null(genes_for_rs_id)) {
      next
    }

    for (gene in genes_for_rs_id) {
      panther_ids_for_gene <- gene_panther_mapping[[gene]]
      if (is.null(panther_ids_for_gene)) {
        next
      }

      for (panther_id in panther_ids_for_gene) {
        panther_id <- trimws(as.character(panther_id))
        if (!nzchar(panther_id)) {
          next
        }

        if (!is.null(panther_ids_set) && !(panther_id %in% panther_ids_set)) {
          next
        }

        pair_key <- paste(rs_id, panther_id, sep = "-")
        if (exists(pair_key, envir = processed_pairs, inherits = FALSE)) {
          next
        }
        assign(pair_key, TRUE, envir = processed_pairs)

        mapped_genes <- genes_for_rs_id[vapply(
          genes_for_rs_id,
          function(mapped_gene) {
            mapped_gene_panther_ids <- gene_panther_mapping[[mapped_gene]]
            !is.null(mapped_gene_panther_ids) && panther_id %in% as.character(mapped_gene_panther_ids)
          },
          logical(1)
        )]

        gene_info <- panther_gene_info[[panther_id]]
        if (is.null(gene_info) || !is.list(gene_info)) {
          next
        }

        table_data[[length(table_data) + 1]] <- list(
          rsId = rs_id,
          PANTHER_ID = panther_id,
          mappedGenes = as.list(mapped_genes),
          PANTHER_family = if (!is.null(gene_info$PANTHER_family)) gene_info$PANTHER_family else "",
          PANTHER_Subfamily = if (!is.null(gene_info$PANTHER_Subfamily)) gene_info$PANTHER_Subfamily else "",
          PANTHER_Pathway = if (!is.null(gene_info$PANTHER_Pathway)) gene_info$PANTHER_Pathway else "",
          Protein_Class = if (!is.null(gene_info$Protein_Class)) gene_info$Protein_Class else "",
          Reactome_Pathway = if (!is.null(gene_info$Reactome_Pathway)) gene_info$Reactome_Pathway else "",
          GO_database_MF_complete = if (!is.null(gene_info$GO_database_MF_complete)) gene_info$GO_database_MF_complete else "",
          GO_database_BP_complete = if (!is.null(gene_info$GO_database_BP_complete)) gene_info$GO_database_BP_complete else "",
          GO_database_CC_complete = if (!is.null(gene_info$GO_database_CC_complete)) gene_info$GO_database_CC_complete else "",
          PANTHER_GO_slim_Molecular_Function = if (!is.null(gene_info$PANTHER_GO_slim_Molecular_Function)) gene_info$PANTHER_GO_slim_Molecular_Function else "",
          PANTHER_GO_slim_Biological_Process = if (!is.null(gene_info$PANTHER_GO_slim_Biological_Process)) gene_info$PANTHER_GO_slim_Biological_Process else "",
          PANTHER_GO_slim_Cellular_Component = if (!is.null(gene_info$PANTHER_GO_slim_Cellular_Component)) gene_info$PANTHER_GO_slim_Cellular_Component else ""
        )
      }
    }
  }

  if (is.null(annotation_dataset)) {
    return(table_data)
  }

  selected_columns <- .get_relevant_columns(annotation_dataset)
  lapply(table_data, function(row) row[selected_columns])
}

.get_significant_results <- function(rows, correction) {
  if (is.null(rows) || !is.list(rows)) {
    return(list())
  }

  if (tolower(as.character(correction)) == "fdr") {
    return(Filter(function(row) {
      !is.null(row$fdr) && as.numeric(row$fdr) < 0.05
    }, rows))
  }

  Filter(function(row) {
    !is.null(row$pValue) && as.numeric(row$pValue) < 0.05
  }, rows)
}

.get_significant_genes <- function(rows) {
  if (is.null(rows) || !is.list(rows)) {
    return(character())
  }

  unique_rows_by_term <- list()
  for (row in rows) {
    key <- as.character(if (!is.null(row$termId) && nzchar(as.character(row$termId))) row$termId else row$process)
    if (!nzchar(key) || !is.null(unique_rows_by_term[[key]])) {
      next
    }
    unique_rows_by_term[[key]] <- row
  }

  genes <- character()
  for (row in unique_rows_by_term) {
    mapped_ids <- row$mapped_ids
    if (is.null(mapped_ids)) {
      next
    }
    for (gene in mapped_ids) {
      normalized_gene <- trimws(as.character(gene))
      if (nzchar(normalized_gene) && !(normalized_gene %in% genes)) {
        genes <- c(genes, normalized_gene)
      }
    }
  }

  genes
}

.build_panther_id_filter_from_genes <- function(gene_panther_mapping, genes_to_include) {
  if (is.null(genes_to_include) || length(genes_to_include) == 0) {
    return(character())
  }

  genes_to_include <- unique(trimws(as.character(genes_to_include)))
  genes_to_include <- genes_to_include[nzchar(genes_to_include)]
  if (length(genes_to_include) == 0) {
    return(character())
  }

  panther_ids <- character()
  for (gene in genes_to_include) {
    panther_ids_for_gene <- gene_panther_mapping[[gene]]
    if (is.null(panther_ids_for_gene)) {
      next
    }

    for (panther_id in panther_ids_for_gene) {
      normalized_panther_id <- trimws(as.character(panther_id))
      if (nzchar(normalized_panther_id) && !(normalized_panther_id %in% panther_ids)) {
        panther_ids <- c(panther_ids, normalized_panther_id)
      }
    }
  }

  panther_ids
}


#' Run SNPWay gene mapping workflow.
#'
#' This function calls the SNPWay backend endpoint that maps SNP input to genes
#' and returns PANTHER mapping metadata.
#'
#' Provide exactly one input mode:
#' - `vcf_text`
#' - `chrom_pos_ids`
#' - chromosome range (`chromosome_identifier`, `start_position`, `end_position`)
#' - `rsid_list`
#'
#' @param vcf_text Raw VCF text. CHROM and POS columns are used.
#' @param chrom_pos_ids Character vector of IDs in `chr:pos` format.
#' @param chromosome_identifier Chromosome label for region mode.
#' @param start_position Region start coordinate.
#' @param end_position Region end coordinate.
#' @param rsid_list RSID list as character vector or comma/space/newline string.
#' @param base_url Optional SNPWay base URL override.
#' @param timeout_seconds HTTP timeout in seconds.
#'
#' @return A list containing SNP-to-gene mappings and PANTHER metadata.
#'
#' @export
snpwayGeneMappingsQuery <- function(vcf_text = NULL,
                                    chrom_pos_ids = NULL,
                                    chromosome_identifier = NULL,
                                    start_position = NULL,
                                    end_position = NULL,
                                    rsid_list = NULL,
                                    base_url = NULL,
                                    timeout_seconds = 120) {
  payload <- .build_snpway_payload(
    vcf_text = vcf_text,
    chrom_pos_ids = chrom_pos_ids,
    chromosome_identifier = chromosome_identifier,
    start_position = start_position,
    end_position = end_position,
    rsid_list = rsid_list
  )

  .post_snpway_request(
    endpoint = "/workflow/gene_mappings",
    payload = payload,
    base_url = base_url,
    timeout_seconds = timeout_seconds
  )
}


#' Run full SNPWay overrepresentation workflow.
#'
#' This function calls the SNPWay backend endpoint that runs mapping,
#' overrepresentation, and CSV-ready download row generation.
#'
#' @param annot_data_set Annotation dataset ID (default: GO:0008150).
#' @param correction Multiple-testing correction (`FDR`, `BONFERRONI`, `NONE`).
#' @param enrichment_test_type Enrichment test type (`FISHER`, `BINOMIAL`).
#' @param vcf_text Raw VCF text. CHROM and POS columns are used.
#' @param chrom_pos_ids Character vector of IDs in `chr:pos` format.
#' @param chromosome_identifier Chromosome label for region mode.
#' @param start_position Region start coordinate.
#' @param end_position Region end coordinate.
#' @param rsid_list RSID list as character vector or comma/space/newline string.
#' @param base_url Optional SNPWay base URL override.
#' @param timeout_seconds HTTP timeout in seconds.
#'
#' @return A list with all overrepresentation results, significant results,
#' and CSV-ready mapping rows.
#'
#' @export
snpwayOverrepresentationWorkflowQuery <- function(annot_data_set = "GO:0008150",
                                                  correction = "FDR",
                                                  enrichment_test_type = "FISHER",
                                                  vcf_text = NULL,
                                                  chrom_pos_ids = NULL,
                                                  chromosome_identifier = NULL,
                                                  start_position = NULL,
                                                  end_position = NULL,
                                                  rsid_list = NULL,
                                                  base_url = NULL,
                                                  timeout_seconds = 300) {
  payload <- .build_snpway_payload(
    vcf_text = vcf_text,
    chrom_pos_ids = chrom_pos_ids,
    chromosome_identifier = chromosome_identifier,
    start_position = start_position,
    end_position = end_position,
    rsid_list = rsid_list
  )

  payload$annotDataSet <- annot_data_set
  payload$correction <- correction
  payload$enrichmentTestType <- enrichment_test_type

  response_data <- .post_snpway_request(
    endpoint = "/workflow/overrepresentation",
    payload = payload,
    base_url = base_url,
    timeout_seconds = timeout_seconds
  )

  overrepresentation_results <- response_data$overrepresentation_results
  if (is.null(overrepresentation_results) || !is.list(overrepresentation_results)) {
    overrepresentation_results <- list()
  }

  significant_results <- .get_significant_results(overrepresentation_results, correction)
  significant_genes <- .get_significant_genes(significant_results)
  significant_panther_ids <- .build_panther_id_filter_from_genes(
    response_data$gene_panther_mapping,
    significant_genes
  )

  response_data$overrepresentation_all_results <- overrepresentation_results
  response_data$overrepresentation_significant_results <- significant_results
  response_data$csv_all_mappings <- .create_results_table_data(
    response_data,
    annotation_dataset = annot_data_set
  )
  response_data$csv_all_mappings_all_columns <- .create_results_table_data(
    response_data,
    annotation_dataset = NULL
  )
  response_data$csv_significant_mappings <- .create_results_table_data(
    response_data,
    panther_ids_to_include = significant_panther_ids,
    annotation_dataset = annot_data_set
  )
  response_data$csv_significant_mappings_all_columns <- .create_results_table_data(
    response_data,
    panther_ids_to_include = significant_panther_ids,
    annotation_dataset = NULL
  )

  response_data
}
