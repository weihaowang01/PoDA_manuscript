library(PD16Sdata)
library(dbplyr)
Wallen=ps_list$Wallen1
# Wallen=ps_list$Wallen2 # Wallen2

otu= as(otu_table(Wallen),"matrix")
meta = sample_data(Wallen)
taxtable=tax_table(Wallen)
rownames(taxtable)=NULL
taxon_name = tax_table(Wallen)[,6]

oo = match( colnames(otu),rownames(taxon_name ))
colnames(otu)=taxon_name[oo]
otu=otu[,which(!is.na(colnames(otu)))]

rdind=which(meta$t==12)
otu=otu[rdind,]
meta=meta[rdind,]
otu1=as.data.frame(otu)

# Ensure that the column names are unique.
df=otu1
colnames(df) <- make.unique(colnames(df))

# Function to merge duplicate columns.
merge_duplicate_columns <- function(df) {
  # Escape special characters such as square brackets.
  escape_regex <- function(x) {
    gsub("([][\\^$.|?*+(){}])", "\\\\\\1", x)
  }
  
  # Extract unique column name prefixes while preserving square brackets
  unique_names <- unique(gsub("\\.\\d+$", "", names(df)))
  df_merged <- data.frame(matrix(nrow = nrow(df), ncol = 0))
  
  for (name in unique_names) {
    # Identify columns with the same prefix
    escaped_name <- escape_regex(name)
    columns_to_merge <- grep(paste0("^", escaped_name, "(\\.\\d+)?$"), names(df), value = TRUE)
    if (length(columns_to_merge) > 1) {
      # Merge those columns
      df_merged[[name]] <- rowSums(df[columns_to_merge], na.rm = TRUE)
    } else {
      # If there are no duplicates, retain the column as is
      df_merged[[name]] <- df[[columns_to_merge]]
    }
  }
  
  return(df_merged)
}

df_merged <- merge_duplicate_columns(df)

keepsam=which(!is.na(meta$PD))
df_merged=df_merged[keepsam,]
meta=meta[keepsam,]

datalist=list(otu=df_merged,meta=meta)
saveRDS(datalist,"~/PD/Wallen1.rds")
# saveRDS(datalist,"~/PD/Wallen2.rds")

