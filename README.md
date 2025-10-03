# AnnoQ GraphQL API Query Functions

This package provides functions to query the AnnoQ GraphQL API for SNP annotations.

## Installation

You can install the package directly from GitHub using the `devtools` package:

```R
install.packages("devtools")
devtools::install_github("USCbiostats/AnnoQR")
```

## Usage

Load the package and use the provided functions to perform queries.

```R
library(AnnoQR)

# Define the annotations to retrieve
annotations_to_retrieve <- c("chr", "pos", "ref", "alt")

# Perform a region query
snps <- regionQuery("18", 1, 50000, annotations_to_retrieve)
head(snps)

# Perform an rsID query
snp <- rsidQuery("rs559687999", annotations_to_retrieve)
print(snp)

# Perform a multiple rsIDs query
snps <- rsidsQuery(c("rs115366554", "rs189126619"), annotations_to_retrieve)
print(snps)
```
