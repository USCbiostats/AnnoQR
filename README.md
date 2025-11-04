# AnnoQR

An R package for programmatically accessing SNP data from the AnnoQ API.

## Installation

Install directly from GitHub using the `devtools` package:

```R
install.packages("devtools")
devtools::install_github("USCbiostats/AnnoQR")
```

## Requirements

- R 3.5 or higher
- Required packages: `httr`, `jsonlite` (automatically installed)

## Quick Start

```R
library(AnnoQR)

# Get available SNP attributes
attributes <- snpAttributesQuery()

# Search SNPs on chromosome 1
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 100000,
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151")
)
```

## Core Functions

The package provides 7 main functions organized into three categories:

### Attribute Discovery

- `snpAttributesQuery()` - List all available SNP attributes

### SNP Retrieval

- `regionQuery()` - Query by chromosome and position range
- `rsidsQuery()` - Query by RSID identifiers
- `geneQuery()` - Query by gene information

### SNP Counting

- `countRegionQuery()` - Count SNPs by chromosome
- `countRsidsQuery()` - Count SNPs by RSID list
- `countGeneQuery()` - Count SNPs by gene

---

## Detailed Usage

### 1. Getting SNP Attributes

Retrieve the list of all available SNP attributes that can be queried.

```R
library(AnnoQR)

# Get all available attributes
attributes <- snpAttributesQuery()

# attributes is a list of attribute metadata
for (i in seq_along(attributes)) {
  cat(sprintf("%s: %s\n", attributes[[i]]$label, attributes[[i]]$description))
}
```

### 2. Querying SNPs by Chromosome

Search for SNPs within a specific chromosome region.

#### Basic Usage

```R
# Query chromosome 1 from position 1 to 100,000 and get basic fields
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 100000,
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151")
)

# Query the X chromosome from position 1,000 to 50,000 and get default fields
snps <- regionQuery(
  chromosome_identifier = "X",
  start_position = 1000,
  end_position = 50000
)
```

#### Selecting Specific Fields

You can specify which fields to return in three different ways:

**As a vector of field names:**

```R
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 10000,
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151")
)
```

**As a string config exported from [AnnoQ](https://annoq.org):**

```R
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 10000,
  fields = '{"_source":["chr", "pos", "ref", "alt", "rs_dbSNP151"]}'
)
```

**From a JSON config exported from [AnnoQ](https://annoq.org):**

```R
# Export the config file: config.txt from AnnoQ
# {"_source":["chr", "pos", "ref", "alt", "rs_dbSNP151"]}

snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 10000,
  fields = "/path/to/config.txt"
)
```

**Note:** The maximum number of fields you can request is **20**. For more fields you can make multiple queries and combine the results.

#### Filtering by Non-Empty Fields

Return only SNPs where specific annotation fields have values:

```R
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 100000,
  filter_fields = c("ANNOVAR_ucsc_Transcript_ID", "VEP_ensembl_Gene_ID")
)
```

#### Pagination

By default, the API returns 1,000 results per page with a maximum of 10,000 results across all pages.

```R
# Get first 500 results
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 1000000,
  pagination_from = 0,
  pagination_size = 500
)

# Get next 500 results
snps_page2 <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 1000000,
  pagination_from = 500,
  pagination_size = 500
)

# Note: pagination_from + pagination_size must be <= 10,000
```

#### Fetching All Results

To retrieve all matching SNPs (up to 1,000,000), use `fetch_all = TRUE`:

```R
# This will download all matching SNPs
all_snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 100000,
  fetch_all = TRUE
)

# When fetch_all = TRUE, the pagination parameters are ignored
```

**Important:** When `fetch_all = TRUE`, the function downloads a lot of data in a different format and may take a long time for large result sets.

### 3. Querying SNPs by RSID

Search for SNPs using RSID identifiers.

#### Basic Usage

```R
# Using a comma-separated string
snps <- rsidsQuery(
  rsid_list = "rs1219648,rs2912774,rs2981582"
)

# Using a vector
snps <- rsidsQuery(
  rsid_list = c("rs1219648", "rs2912774", "rs2981582")
)
```

#### With Custom Fields

```R
snps <- rsidsQuery(
  rsid_list = c("rs1219648", "rs2912774", "rs2981582"),
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151")
)
```

#### With Filtering

```R
snps <- rsidsQuery(
  rsid_list = "rs1219648,rs2912774,rs2981582",
  filter_fields = c("VEP_ensembl_Gene_ID"),
  pagination_from = 0,
  pagination_size = 100
)
```

#### Fetching All Matching RSIDs

```R
# Get all SNPs for a large list of RSIDs
all_snps <- rsidsQuery(
  rsid_list = c("rs1219648", "rs2912774", "rs2981582", "rs123456", "rs789012"),
  fetch_all = TRUE
)
```

### 4. Querying SNPs by Gene Product

Search for SNPs associated with a gene using gene ID, gene symbol, or UniProt ID.

#### Basic Usage

```R
# Search by gene symbol
snps <- geneQuery(gene = "BRCA1")

# Search by gene ID or UniProt ID
snps <- geneQuery(gene = "ENSG00000012048")
```

#### With Custom Fields and Filtering

```R
snps <- geneQuery(
  gene = "TP53",
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151"),
  filter_fields = c("ANNOVAR_ucsc_Transcript_ID")
)
```

#### With Pagination

```R
# Get first 500 SNPs for a gene
snps <- geneQuery(
  gene = "APOE",
  pagination_from = 0,
  pagination_size = 500
)
```

#### Fetching All Gene-Associated SNPs

```R
# Get all SNPs associated with a gene
all_snps <- geneQuery(
  gene = "ZMYND11",
  fetch_all = TRUE
)
```

### 5. Counting SNPs

Count functions return the number of matching SNPs without retrieving the actual data.

#### Count by Chromosome

```R
# Count all SNPs in a region
count <- countRegionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 100000
)
cat(sprintf("Found %d SNPs\n", count))

# Count with filters
count <- countRegionQuery(
  chromosome_identifier = "X",
  start_position = 1000,
  end_position = 50000,
  filter_fields = c("VEP_ensembl_Gene_ID", "ANNOVAR_ucsc_Transcript_ID")
)
```

#### Count by RSID List

```R
# Count matching RSIDs
count <- countRsidsQuery(
  rsid_list = c("rs1219648", "rs2912774", "rs2981582")
)

# Count with filters
count <- countRsidsQuery(
  rsid_list = "rs1219648,rs2912774,rs2981582",
  filter_fields = c("ANNOVAR_ucsc_Transcript_ID")
)
```

#### Count by Gene Product

```R
# Count SNPs for a gene
count <- countGeneQuery(gene = "BRCA1")

# Count with filters
count <- countGeneQuery(
  gene = "TP53",
  filter_fields = c("VEP_ensembl_Gene_ID")
)
```

---

## Common Patterns

### Example 1: Progressive Filtering

```R
# First, count to see how many SNPs match
total <- countRegionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 1000000
)
cat(sprintf("Total SNPs: %d\n", total))

# Count with filters applied
filtered_count <- countRegionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 1000000,
  filter_fields = c("VEP_ensembl_Gene_ID")
)
cat(sprintf("Filtered SNPs: %d\n", filtered_count))

# Retrieve the filtered data
snps <- regionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 1000000,
  filter_fields = c("VEP_ensembl_Gene_ID"),
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151", "VEP_ensembl_Gene_ID")
)
```

### Example 2: Working with Large Datasets

```R
# For large regions, first check the count
count <- countRegionQuery(
  chromosome_identifier = "1",
  start_position = 1,
  end_position = 10000000
)

if (count > 1000000) {
  cat(sprintf("Warning: %d SNPs found. Consider narrowing your search.\n", count))
} else if (count > 10000) {
  # Use fetch_all for counts between 10K and 1M
  snps <- regionQuery(
    chromosome_identifier = "1",
    start_position = 1,
    end_position = 10000000,
    fetch_all = TRUE
  )
} else {
  # Use regular pagination for smaller datasets
  snps <- regionQuery(
    chromosome_identifier = "1",
    start_position = 1,
    end_position = 10000000,
    pagination_size = count  # Get all in one go
  )
}
```

### Example 3: Gene-Focused Analysis

```R
# Get all SNPs for multiple genes
genes <- c("BRCA1", "BRCA2", "TP53")
all_gene_snps <- list()

for (gene in genes) {
  count <- countGeneQuery(gene = gene)
  cat(sprintf("%s: %d SNPs\n", gene, count))
  
  all_gene_snps[[gene]] <- geneQuery(
    gene = gene,
    fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151"),
    fetch_all = TRUE
  )
}
```

### Example 4: Batch RSID Lookup

```R
# Read RSIDs from a file
rsids <- readLines("rsid_list.txt")
rsids <- rsids[nchar(rsids) > 0]  # Remove empty lines

# Check how many exist in the database
count <- countRsidsQuery(rsid_list = rsids)
cat(sprintf("%d out of %d RSIDs found\n", count, length(rsids)))

# Retrieve all matching SNPs
snps <- rsidsQuery(
  rsid_list = rsids,
  fields = c("chr", "pos", "ref", "alt", "rs_dbSNP151"),
  fetch_all = TRUE
)
```

---

## Important Limitations

### Pagination Constraints

- **Regular queries**: Maximum of 10,000 results across all pages (`pagination_from + pagination_size <= 10,000`)
- **Fetch all queries**: Maximum of 1,000,000 total results
- **Note**: For large datasets, the results may be too large and could lead to performance issues. It is recommended to narrow down the query if possible.

### Field Selection

- Maximum of **20 fields** can be requested per query
- Use the `snpAttributesQuery()` function to see all available fields

### Rate Limiting

- The API may implement rate limiting for excessive requests
- Use count functions before large retrievals to estimate data size

---

## Error Handling

All functions raise errors for common error cases:

```R
# Pagination error
tryCatch({
  snps <- regionQuery(
    chromosome_identifier = "1",
    start_position = 1,
    end_position = 100000,
    pagination_from = 9500,
    pagination_size = 1000  # This exceeds the 10,000 limit
  )
}, error = function(e) {
  cat(sprintf("Pagination error: %s\n", e$message))
})

# File error
tryCatch({
  snps <- regionQuery(
    chromosome_identifier = "1",
    start_position = 1,
    end_position = 100000,
    fields = "/nonexistent/file.json"
  )
}, error = function(e) {
  cat(sprintf("File error: %s\n", e$message))
})

# API error
tryCatch({
  snps <- regionQuery(
    chromosome_identifier = "invalid",
    start_position = 1,
    end_position = 100000
  )
}, error = function(e) {
  cat(sprintf("API error: %s\n", e$message))
})
```

---

## Contributing

Contributions are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on the [GitHub repository](https://github.com/USCbiostats/AnnoQR).

## License

This package is licensed under the MIT License.

## Support

For questions or issues related to AnnoQ itself, please visit the site [AnnoQ](https://annoq.org)
