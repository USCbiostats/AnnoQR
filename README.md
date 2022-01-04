# AnnoQR: R client for AnnoQ Variant Query

## Introdcution

This is an R client for performing queries with AnnoQ API.

## Installation

Install from github.
Make sure you have installed `devtools`. 

`install.packages("devtools")`

Then 

`library(devtools)`

` install_github("USCbiostats/AnnoQR")`

## Function list

* add\_query\_filter
* add\_source
* exists\_filter
* init\_query\_js\_body
* keywordsQuery
* perform\_search
* query\_obj\_to\_json
* range\_filter
* read\_config
* regionQuery
* rsidQuery
* term\_filter

## Examples

Query Variants with `ANNOVAR_ensembl_Effect ` Annotation

```R
library(AnnoQR)
q = init_query_js_body()
ex = exists_filter("ANNOVAR_ensembl_Effect")
q = add_query_filter(q, ex)
variants = perform_search(q)
variants
```




Only retrieve `ANNOVAR_ensembl_Effect ` column

```R
q = add_source(q, c("ANNOVAR_ensembl_Effect"))
variants = perform_search(q)
variants
```



Query variants field `SnpEff_ensembl_Effect` marked as `intergenic_region`

```R
q = init_query_js_body()
term_f = term_filter('SnpEff_ensembl_Effect' , 'intergenic_region')
q = add_query_filter(q, term_f)
variants = perform_search(q)
variants
```


Query variants field `SnpEff_ensembl_Effect` marked as `intergenic_region with in chromosome 20

```R
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

```R
q = init_query_js_body()
range_f = range_filter(key='1000Gp3_AC' , gt=5)
q = add_query_filter(q, range_f)
variants = perform_search(q)
variants
```

Chromosome range query

```R
variants = regionQuery(contig = '20', start=31710367, end=31820367)
variants
```

rsID query

```R
variant = rsidQuery('rs193031179')
variant
```

keywordsQuery

```R
keywordsQuery('protein_coding')
```
