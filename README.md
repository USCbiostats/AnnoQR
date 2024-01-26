# AnnoQR: R client for AnnoQ Variant Query

## Introduction

This is an R client for performing queries with AnnoQ API.

## Installation

Install from github. Make sure you have installed `devtools`.

`install.packages("devtools")`

Then

`library(devtools)`

`install_github("USCbiostats/AnnoQR")`

## Function list

-   add_query_filter
-   add_source
-   exists_filter
-   init_query_js_body
-   keywordsQuery
-   perform_search
-   query_obj_to_json
-   range_filter
-   read_config
-   regionQuery
-   rsidQuery
-   term_filter

## Examples

Query Variants with `ANNOVAR_ensembl_Effect` Annotation

``` r
library(AnnoQR)
q = init_query_js_body()
ex = exists_filter("ANNOVAR_ensembl_Effect")
q = add_query_filter(q, ex)
variants = perform_search(q)
variants
```

Only retrieve `ANNOVAR_ensembl_Effect` column

``` r
q = add_source(q, c("ANNOVAR_ensembl_Effect"))
variants = perform_search(q)
variants
```

Query variants field `SnpEff_ensembl_Effect` marked as `intergenic_region`

``` r
q = init_query_js_body()
term_f = term_filter('SnpEff_ensembl_Effect' , 'intergenic_region')
q = add_query_filter(q, term_f)
variants = perform_search(q)
variants
```

Query variants field `SnpEff_ensembl_Effect` marked as \`intergenic_region with in chromosome 20

``` r
q = init_query_js_body()
term_f1 = term_filter('SnpEff_ensembl_Effect' , 'intergenic_region')
term_f2 = term_filter('chr' , '20')
q = add_query_filter(q, term_f1)
q = add_query_filter(q, term_f2)
#q = add_source(q, c('SnpEff_ensembl_Effect'))
variants = perform_search(q)
variants
```

Query variants with 1000 genome allel count `1000Gp3_AC` larger than 5

``` r
q = init_query_js_body()
range_f = range_filter(key='1000Gp3_AC' , gt=5)
q = add_query_filter(q, range_f)
variants = perform_search(q)
variants
```

Chromosome range query

``` r
variants = regionQuery(contig = '20', start=31710367, end=31820367)
variants
```

rsID query

``` r
variant = rsidQuery('rs193031179')
variant
```

keywordsQuery

``` r
keywordsQuery('protein_coding')
```

------------------------------------------------------------------------

## **Guidance on Using Our Elasticsearch-based API**

Our API leverages the powerful features of Elasticsearch, but it's important to be aware of certain behaviors related to query results:

### Default Behavior with `perform_search(q)`

-   **Limited Results**: Utilizing `perform_search(q)` by default yields only the first 10 matches. This is a standard constraint imposed by Elasticsearch on query results.
-   **Ideal for Quick Queries**: This function is optimal for concise queries where a limited dataset suffices or where only a preview of results is needed.

For example, the following R code snippet demonstrates how to use `perform_search(q)` effectively:

``` r
q = init_query_js_body()
ex = exists_filter("ANNOVAR_ensembl_Effect")
q = add_query_filter(q, ex)
variants = perform_search(q)
hits = variants$hits$hits
length(hits$`_index`)
```

Running this snippet typically results in:

```         
10
```

This indicates the successful retrieval of the first 10 matches, aligning with Elasticsearch's default result limit.

### Retrieving All Matches with `perform_search_with_count(q)`

-   **For Comprehensive Results**: To fetch all corresponding matches for a query, use `perform_search_with_count(q)`. This is particularly useful for exhaustive queries where the entire dataset is necessary.
-   **Handling Large Result Sets**: Be cautious with queries matching a large number of documents (over 10,000 hits), as this may lead to an HTTP 400 error. This occurs due to Elasticsearch's cap on the number of results returned in a single query.

Consider this code snippet:

``` r
q = init_query_js_body()
ex = exists_filter("ANNOVAR_ensembl_Effect")
q = add_query_filter(q, ex)
variants = perform_search_with_count(q)
```

This may result in an error if the result set is too large:

```         
Error in perform_search_with_count(q) : Bad Request (HTTP 400)
```

The error is attributed to the attempt of `perform_search_with_count(q)` to retrieve all matches, surpassing Elasticsearch's maximum limit.

### Diagnosing Large Queries with `perform_search_find_count(q)`

To ascertain the size of your query's result set, you can use:

``` r
q = init_query_js_body()
ex = exists_filter("ANNOVAR_ensembl_Effect")
q = add_query_filter(q, ex)
variants = perform_search_find_count(q)
```

This will generate:

```         
Debug: Query JSON:
 {"query":{"bool":{"filter":[{"exists":{"field":"ANNOVAR_ensembl_Effect"}}]}},"size":40405505} 
Response [http://annoq.org/api/annoq-test/_search]
  Date: 
  Status: 400
  Content-Type: application/json;charset=utf-8
  Size: 1.48 kB
```

Here, the "size" parameter exceeds 40 million, explaining the error. Such a large result set is beyond the permissible range for a single Elasticsearch query.

### Managing Large Datasets

-   **Total Match Assessment**: If your query is likely to yield an extensive number of matches, start with `perform_search_find_count(q)` to determine the total count.
-   **Pagination Development (In Progress)**: Recognizing the necessity to manage queries returning millions of results, we are developing a pagination function. This will facilitate accessing large datasets in smaller, sequential segments.
-   **Your Feedback Matters**: Your patience and input are invaluable as we strive to enhance our API. We are dedicated to continuously improving our services to better suit your needs.

