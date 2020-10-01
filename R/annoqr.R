#' @import jsonlite httr
#'
NULL

host = 'http://annoq.org:3404'

#' init_query_js_body
#'
#' This function init a empty query body is R list format
#'
#' @usage init_query_js_body()
#' @return A list
#'
#' @export
init_query_js_body <- function(){
  q = list()
  q[['query']] = list()
  q[['query']][['bool']] = list()
  q[['query']][['bool']][['filter']] =list()
  q
}

#' add_source
#'
#' This function adds a source field to a query body in R obj
#'
#' @usage add_source(q, source)
#' @param q query body
#' @param source a vector contain source field names
#' @return A query body
#'
#' @export
add_source <- function(q, source){
  q[['_source']] = as.list(source)
  q
}

#' exists_filter
#'
#' This function returns filter depends on the existance of a certain field
#'
#' @usage exists_filter(key)
#' @param key A string, name for that field.
#' @return A filter
#'
#' @export
exists_filter <- function(key){
  q = list()
  q[['exists']] = list()
  q[['exists']][['field']] = key
  q
}

#' term_filter
#'
#' This function returns a filter to choose a certain field is exactly some value
#'
#' @usage term_filter(key, value)
#' @param key A string, name for that field.
#' @param value A string, value for that field.
#' @return A filter
#'
#' @export
term_filter <- function(key, value){
  q = list()
  q[['term']] = list()
  q[['term']][[key]] = tolower(value)
  q
}

#' range_filter
#'
#' This function return a filter to choose a certain numeric field within a range.
#'
#' @usage range_filter(key, gt, lt)
#' @param key A string, name for that field.
#' @param gt A number, optional, open lower boundary for that field
#' @param lt A number, optional, open upper boundary for that field
#' @return A filter
#'
#' @export
range_filter <- function(key, gt=NULL, lt=NULL){
  q = list()
  if (is.null(gt)&is.null(lt)) {
    return
  }
  q[['range']] = list()
  q[['range']][[key]] = list()
  if (!is.null(gt)) {
    q[['range']][[key]][['gt']] = gt
  }
  if (!is.null(gt)) {
    q[['range']][[key]][['lt']] = lt
  }
  q
}

#' read_config
#'
#' Reading source configure from a file
#'
#' @usage read_config(file_name)
#' @param file_name A string, full path to that configure file.
#' @return A vector contain source fiedl names
#'
#' @export
read_config <- function(file_name){
  unlist(read_json(file_name)[['_source']])
}

#' add_query_filter
#'
#' This function binds a filter to a query body.
#'
#' @usage add_query_filter(q, filter)
#' @param q A query object.
#' @param filter A filter object.
#' @return A query body
#'
#' @export
add_query_filter <- function(q, filter) {
  old = q[['query']][['bool']][['filter']]
  l = length(old)
  if ( l==0 ){
    q[['query']][['bool']][['filter']] = list(filter)
  }
  else {
    q[['query']][['bool']][['filter']][[l + 1]] = filter
  }
  q
}

#' query_obj_to_json
#'
#' This function convert a query obj to a json string.
#'
#' @usage query_obj_to_json(q)
#' @param q A qeury obj
#' @return A json string.
#'
#' @export
query_obj_to_json <- function(q){
  jsonlite::toJSON(q, auto_unbox = T)
}

#' perform_search
#'
#' execute http request with a query json string
#'
#' @usage perform_search(q)
#' @param q A qeury string
#' @return A list contain variants
#'
#' @export
perform_search <- function(q) {
  r <- POST(paste0(host, "/vs-index/_search"), content_type_json(), body = query_obj_to_json(q))
  stop_for_status(r)
  content(r, "parsed", "application/json")
}

#' regionQuery
#' @usage regionQuery(contig, start, end)
#' @description This function use a genome coordinate to query variants within that region
#' @param contig contig string
#' @param start  sgtart position in int
#' @param end end position in int
#' @param configFile  optional parameter that provide fields to return, FALSE by default
#' @return A list contain variants
#'
#' @export
regionQuery <- function(contig, start, end, configFile = FALSE) {
  body = init_query_js_body()
  chr = term_filter('chr', contig)
  pos = range_filter('pos', gt=start, lt=end)
  if (!isFALSE(configFile)) {
    body = add_source(body, read_config(configFile))
  }
  body = add_query_filter(body, chr)
  body = add_query_filter(body, pos)
  perform_search(body)
}

#' rsidQuery
#' @usage rsidQuery(rsid)
#' @description This function use a rsid to query a variant
#' @param rsid rsid a string
#' @param configFile  optional parameter that provide fields to return, FALSE by default
#' @return One variant if found
#' @export
rsidQuery <- function(rsid){
  body = init_query_js_body()
  rs_filter = term_filter('rs_dbSNP151', rsid)
  body = add_query_filter(body, rs_filter)
  perform_search(body)
}


#' keywordsQuery
#' @usage keywordsQuery(keywords)
#' @description Perform full text search in our variant annotation dataset
#' @param keywords keywords a string
#' @param configFile  optional parameter that provide fields to return, FALSE by default
#' @return A list contain variants
#' @export
keywordsQuery <- function(keywords){
  body = '{
  "query": {
  "multi_match": {"query":"Signaling by GPCR"}
  }}'
  body = parse_json(body)
  body[['query']][['multi_match']][['query']] = keywords
  r <- POST(paste0(host, "/vs-index/_search"), content_type_json(), body = toJSON(body, auto_unbox = T))
  stop_for_status(r)
  content(r, "parsed", "application/json")
}


