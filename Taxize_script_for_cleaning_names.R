## Taxize script for UNITE

setwd("~/Documents/FunDive/UNITE_taxize_cleanup")
library(tidyverse)

df <- read.table("Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax.fasta")

df_split <- as.data.frame(matrix(df$V1, 
                                 ncol = 2, 
                                 byrow = TRUE, 
                                 dimnames = list(NULL, c('Accession', 'Seq'))))

df_split_clean_names <- df_split %>% 
  separate(Accession, c("Accession", "Taxonomy"), sep = ";tax=k:") %>% 
  separate(Taxonomy, c("Kingdom", "Phylum" , "Class", "Order", "Family", "Genus", "Species"), 
           sep = ",[pcofgs]:")


# Split eukaryote, plant and fungal databases

df_fungi <- df_split_clean_names %>% 
  filter(Kingdom == "Fungi")

df_plant <- df_split_clean_names %>% 
  filter(Kingdom == "Viridiplantae")

df_other <- df_split_clean_names %>% 
  filter(!Kingdom %in% c("Fungi", "Viridiplantae"))

# Set up taxize on a subset of each database 

df_fungi_sub <- df_fungi

df_plant_sub <- df_plant[401:600, ]

tax <- df_fungi

# ==========================================
# 0. Libraries
# ==========================================
# library(dplyr)
# library(stringr)
# library(taxize)
# library(tibble)

# ==========================================
# 1. Preprocess species names
# ==========================================

# Convert Genus_species → Genus species
tax <- tax %>%
  mutate(Species = gsub("_", " ", Species))

# Split species vs genus-only (sp)
tax_species <- tax %>% 
  filter(!grepl(" sp$", Species))

tax_genus <- tax %>%
  filter(grepl(" sp$", Species)) %>%
  mutate(Species = gsub(" sp$", "", Species))

species_unique <- unique(tax_species$Species)
genus_unique   <- unique(tax_genus$Species)

# ===============================================================================
# 2. Index Fungorum verification, batched to avoid errors in large query vectors
# ===============================================================================

run_gna_batched <- function(names_vec, batch_size = 50, data_sources = 5) {
  split(names_vec, ceiling(seq_along(names_vec) / batch_size)) |>
    purrr::map_dfr(~ gna_verifier(names = .x, data_sources = data_sources))
}

ver_species <- function(species_vector) {
  res <- run_gna_batched(species_vector, data_sources = 5)
  
  res %>%
    dplyr::group_by(submittedName) %>%
    dplyr::slice(1) %>%   # keep best IF hit
    dplyr::ungroup() %>%
    dplyr::transmute(
      Species = submittedName,
      canonical_simple = currentCanonicalSimple,
      canonical_full = currentCanonicalFull,
      matched = matchedCanonicalFull,
      taxonomicStatus  = taxonomicStatus,
      name_changed     = is.na(currentCanonicalSimple) |
        submittedName != currentCanonicalSimple
    )
}

ver_genus <- function(genus_vector) {
  res <- run_gna_batched(genus_vector, data_sources = 5)
  
  res %>%
    dplyr::group_by(submittedName) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      Genus = submittedName,
      canonical_simple = currentCanonicalSimple,
      canonical_full = currentCanonicalFull,
      matched = matchedCanonicalFull,
      taxonomicStatus  = taxonomicStatus,
      name_changed     = is.na(currentCanonicalSimple) |
        submittedName != currentCanonicalSimple
    )
}


species_ver <- ver_species(species_unique)
table(species_ver$name_changed)
genus_ver   <- ver_genus(genus_unique)
table(genus_ver$name_changed)

# This concludes the check for synonyms against IndexFungorum

# ==========================================
# 4. Map results back to full taxonomy table
# ==========================================

tax_species_verified <- tax_species %>%
  left_join(species_ver, by = "Species")

tax_genus_verified <- tax_genus %>%
  left_join(genus_ver, by = c("Species" = "Genus"))

tax_verified <- bind_rows(
  tax_species_verified,
  tax_genus_verified
)

# This table now contains the original taxonomy still linked to the Accession number and DNA sequence

# ==========================================
# 5. Extract quality control candidates (manual review)
# ==========================================

qc_candidates <- tax_verified %>%
  filter(name_changed) %>%
  select(
    Accession,
    Seq,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    canonical_simple,
    canonical_full,
    taxonomicStatus
  ) %>%
  distinct()

# Extract the Accessions that had their name changed for review

# Many changes are new epithets int he same Genus, or an updated genus name which may be represented in UNITE in other accessions

# Link the changed names to an updated taxonomy:

update_taxonomy <- qc_candidates %>% select(Accession, Seq, Species, canonical_simple, canonical_full) %>% 
  separate(Species, c("Genus_UNITE", "Species_UNITE", sep = " ")) %>% 
  separate(canonical_simple, c("Genus_if", "Species_if"), sep =" ")

update_taxonomy$Same_genus <- ifelse(update_taxonomy$Genus_UNITE == update_taxonomy$Genus_if, "YES", "NO")

table(update_taxonomy$Same_genus)

# update the df_fungi object with the new epithets for these Accessions
same_genus <- update_taxonomy %>% 
  filter(Same_genus == "YES")

table(is.na(same_genus$Species_if))
# some of these have an NA species annotation from IF but a species annotation from UNITE
# Keep the species annotation from UNITE
for (i in 1:length(same_genus$Species_if)){
  if(is.na(same_genus$Species_if[i])){
    same_genus$Species_if[i] <- same_genus$Species_UNITE[i]
  }
  else{
    same_genus$Species_if[i] <- same_genus$Species_if[i] 
  }
}
table(is.na(same_genus$Species_if))

# Now set the updated species names in the original database object by Accession
# When the genus is the same and only the epithet changed, simply change the epithet

df_fungi_temp1 <- df_fungi %>% filter(!Accession %in% same_genus$Accession)
df_fungi_temp1$Changed <- "NO" # Add a column to be able to track which accessions changed

df_fungi_to_update <- df_fungi %>% filter(Accession %in% same_genus$Accession)
df_fungi_to_update$Species <- paste0(same_genus$Genus_if, "_", same_genus$Species_if)
df_fungi_to_update$Changed <- "YES"
df_fungi_updated1 <- bind_rows(df_fungi_temp1, df_fungi_to_update)

# Next step, check the Accessions where the genus changed

qc_genus_changed <- update_taxonomy %>% filter(Same_genus == "NO")

# Some of the genera may be represented still in UNITE. Rebuild the relevant taxonomy and update the original df

qc_genus_changed$Genus_in_UNITE <- "NO"

qc_genus_changed$Genus_in_UNITE[which(qc_genus_changed$Genus_if %in% df_fungi$Genus)] <- "YES"

table(qc_genus_changed$Genus_in_UNITE) 

# The majority of these genera are in UNITE. 

# Before updating, check for NAs in the species epithet

table(is.na(qc_genus_changed$Species_if))

# As these has had a genus change they cant simply inherit the epithet. Assign "sp"

qc_genus_changed$Species_if[is.na(qc_genus_changed$Species_if)] <- "sp"

parent_taxonomy <- df_fungi %>% filter(Genus %in% qc_genus_changed$Genus_if) %>% 
  group_by(Genus) %>% 
  select(-Accession, -Seq, -Species) %>% 
  slice(1)

genus_changed_updated_taxonomy <- qc_genus_changed %>% 
  filter(qc_genus_changed$Genus_in_UNITE == "YES") %>% 
  select(Accession, Seq, Genus_if, Species_if) %>% 
  rename(Species = Species_if,
         Genus = Genus_if) %>% 
  mutate(Species = paste0(Genus, "_", Species)) %>% 
  left_join(parent_taxonomy, by = "Genus", relationship = "many-to-many")

genus_changed_updated_taxonomy$Changed <- "YES"

# Add to the original database
df_fungi_temp2 <- df_fungi_updated1 %>% 
  filter(!Accession %in% genus_changed_updated_taxonomy$Accession)

df_fungi_updated2 <- bind_rows(df_fungi_temp2 ,genus_changed_updated_taxonomy)


# Now the new names changed from the Index Fungorum scans that were NOT already represented in UNITE needs to be classified

qc_new_name <- qc_genus_changed %>% filter(Genus_in_UNITE == "NO")


# ===============================
#   5. Classify by taxize
# ===============================

# Get the unique names
taxonomy_query <- unique(qc_new_name$canonical_full)

taxonomy_list <- classification(
  taxonomy_query,
  db = "gbif"
)


# ===============================
#   6. Combine the tables for QC
# ===============================


# define the canonical ranks we want
desired_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

gbif_taxonomy <- data.frame(matrix(ncol = length(desired_ranks), nrow = length(taxonomy_list)))
colnames(gbif_taxonomy) <- desired_ranks
taxonomy_list[1]

for (i in 1:length(taxonomy_list)){
  if(is.na(taxonomy_list[i])){
    gbif_taxonomy[i,1:6] = NA
    gbif_taxonomy[i, 7] = names(taxonomy_list[i])
  }
  else{
    taxonomy <- taxonomy_list[[i]] %>% 
      filter(rank %in% desired_ranks)
    taxonomy <- taxonomy$name
    gbif_taxonomy[i, 1:6] = taxonomy[1:6]
    gbif_taxonomy[i, 7] = names(taxonomy_list[i])
  }
}

gbif_taxonomy <- drop_na(gbif_taxonomy)

# merge the names extracted from gbif with the new names

qc_new_names_updated <- qc_new_name %>% 
  select(Accession, Seq, canonical_full) %>% 
  rename(species = canonical_full) %>% 
  left_join(gbif_taxonomy, by = "species") %>% 
  rename(Kingdom = kingdom,
         Phylum = phylum,
         Class = class,
         Order = order, 
         Family = family, 
         Genus = genus,
         Species = species)
qc_new_names_updated$Changed <- "YES"

# merge with the original data
df_fungi_temp3 <- df_fungi_updated2 %>% 
  filter(!Accession %in% qc_new_names_updated$Accession)

df_fungi_updated3 <- bind_rows(df_fungi_temp3 ,qc_new_names_updated)

length(which(df_fungi_updated3$Changed == "YES"))

df_fungi_updated <- df_fungi_updated3

changed_by_order <- data.frame(table(df_fungi_updated$Changed, df_fungi_updated$Order))

# Now all the accessions has been screened for outdated synonyms and updated where applicable 

# next up, check for proper dichotomous branching in the names

# Each identity may only have one parent. Each parent may have many children

# Check by rank:
names(df_fungi_updated)

## Phyla check
Phyla_deviant <- df_fungi_updated %>% 
  distinct(Phylum, Kingdom) %>%
  count(Phylum, name = "n_parents") %>%
  filter(n_parents > 1) # empty data frame. 


## Class check
Class_deviant <- df_fungi_updated %>% 
  distinct(Class, Phylum) %>%
  count(Class, name = "n_parents") %>%
  filter(n_parents > 1) # empty data frame


## Order check
Order_deviant <- df_fungi_updated %>% 
  distinct(Order, Class) %>%
  count(Order, name = "n_parents") %>%
  filter(n_parents > 1) # one order with more than one parent

Order_fix <- df_fungi_updated %>% 
  filter(Order %in% Order_deviant$Order) %>% 
  select(Accession, Kingdom, Phylum, Class, Order) 
# Paraglomerales has one entry with the Class "Glomeromycetes" and 175 entries with the Class "Paraglomeromycetes"
# Set the deviant entry to Paraglomeromycetes

Acc_fix <- Order_fix$Accession[Order_fix$Class == "Glomeromycetes"]

df_fungi_updated$Class[df_fungi_updated$Accession == Acc_fix] <- "Paraglomeromycetes"
# Rerun the dichotomy check
Order_deviant <- df_fungi_updated %>% 
  distinct(Order, Class) %>%
  count(Order, name = "n_parents") %>%
  filter(n_parents > 1) # empty data frame

## Family check
Family_deviant <- df_fungi_updated %>% 
  distinct(Family, Order) %>%
  count(Family, name = "n_parents") %>%
  filter(n_parents > 1) # Four families with more than two parents each
Family_deviant$Family

# Run down the families one by one:

# Australiascaceae
Family_fix <- df_fungi_updated %>% filter(Family == "Australiascaceae") %>% 
  select(Accession, Kingdom, Phylum, Class, Order, Family)
# Family Australiascaceae is placed in two orders; "Glomerellales" and "Chaetosphaeriales"
# Index Fungorum has the correct order as Chaetosphaeriales

df_fungi_updated$Order[df_fungi_updated$Family == "Australiascaceae"] <- "Chaetosphaeriales"

# Endochytriaceae
Family_fix <- df_fungi_updated %>% filter(Family == "Endochytriaceae") %>% 
  select(Accession, Kingdom, Phylum, Class, Order, Family)
# Family Endochytriaceae is placed in two orders; "Cladochytriales" and "Chytridiales"
# Index Fungorum has the correct order as Chytridiales

df_fungi_updated$Order[df_fungi_updated$Family == "Endochytriaceae"] <- "Chytridiales"

# Graphidaceae
Family_fix <- df_fungi_updated %>% filter(Family == "Graphidaceae") %>% 
  select(Accession, Kingdom, Phylum, Class, Order, Family)
# Family Graphidaceae is placed in two orders; "Ostropales" and "	Graphidales"
# Index Fungorum has the correct order as Ostropales while MycoBank used Graphidales. 
# Set to Ostropales as curation has been following Index Fungorum

df_fungi_updated$Order[df_fungi_updated$Family == "Graphidaceae"] <- "Ostropales"

# Synchytriaceae
Family_fix <- df_fungi_updated %>% filter(Family == "Synchytriaceae") %>% 
  select(Accession, Kingdom, Phylum, Class, Order, Family)
# Family Graphidaceae is placed in two orders; "Chytridiales" and "Synchytriales"
# Index Fungorum has the correct order as Chytridiales 

df_fungi_updated$Order[df_fungi_updated$Family == "Synchytriaceae"] <- "Chytridiales"

## Genus check
Genus_deviant <- df_fungi_updated %>% 
  distinct(Genus, Family) %>%
  count(Genus, name = "n_parents") %>%
  filter(n_parents > 1) # empty data frame

## Species check
Species_deviant <- df_fungi_updated %>% 
  distinct(Species, Genus) %>%
  count(Species, name = "n_parents") %>%
  filter(n_parents > 1) # Empty data frame


# There is some conflict in Candida and Entoloma sharing the name of a subsection leading to wrongly annotated Candida:
Candida_error <- df_fungi_updated %>% 
  filter(Changed == "YES") %>% 
  filter(Genus == "Entoloma") %>% 
  rename(Kingdom_n = Kingdom,
         Phylum_n = Phylum,
         Class_n = Class, 
         Order_n = Order, 
         Family_n = Family,
         Genus_n = Genus,
         Species_n = Species)

Candida_orig <- df_fungi %>% 
  filter(Accession %in% Candida_error$Accession) %>% 
  left_join(Candida_error, by = "Accession")

Candida_orig_Entoloma_error <- Candida_orig %>% 
  filter(Phylum == "Ascomycota")
# All of these are erraneous due to an issue where Index Fungorum routes  "Candida_sp" into "Entoloma sect. Candida"
# Remove from the changed list and keep original annotation

df_fungi_keep <- df_fungi %>% filter(Accession %in% c(Candida_orig_Entoloma_error$Accession, ">SH0880858.10FU")) %>% 
  add_column(Changed = "NO")

df_fungi_updated_temp <- df_fungi_updated %>% 
  filter(!Accession %in% c(Candida_orig_Entoloma_error$Accession, ">SH0880858.10FU")) %>% 
  bind_rows(df_fungi_keep) # Check that the number of entries are the same in the temp file before updating the original

df_fungi_updated <- df_fungi_updated_temp

# We have found that especially within Saccharomycetes there has been a split between Saccaromycetales and Serinales not caught in the current database
# Additionally, Index Fungorum has Candida within Serinales, currently in db set to Saccharomycetales_incertae_sedis

# Update Candida family:
df_fungi_updated$Family[df_fungi_updated$Genus == "Candida"] <- "Debaryomycetaceae" 

# Fix this as well to avoid inflated estimates of "Saccharomycetales" diversity or loss of "Serinales" diversity

df_saccharomycetales <- df_fungi_updated %>% filter(Order == "Saccharomycetales")

unique(df_saccharomycetales$Family)

# Move to orders:
Serinales <- c("Metschnikowiaceae", "Debaryomycetaceae")

# Running each family through index fungorum places only "Metschnikowiaceae" in Serinales while the remaining families stay within Saccharomycetales
df_fungi_updated$Order[df_fungi_updated$Family %in% Serinales] <- "Serinales"

# Check that this does not inrotoduce problems with Serinales:
df_serinales <- df_fungi_updated %>% filter(Order == "Serinales")
unique(df_serinales$Kingdom)
unique(df_serinales$Phylum)
unique(df_serinales$Class)
# no issues

# Rebuild the database fasta file keeping the SINTAX names clean:

head(df)

# First, tag the fungal accessions that have been changed
df_fungi_changed <- df_fungi_updated %>% filter(Changed == "YES") %>% 
  mutate(Accession = paste0(.$Accession, "_c")) %>% 
  select(-Changed)

# Make a table of old and new names for reference
df_fungi_old <- df_fungi %>% 
  filter(Accession %in% df_fungi_updated$Accession[df_fungi_updated$Changed == "YES"]) %>% 
  add_column(OLD_NAME = paste0(.$Kingdom, ",",
                               .$Phylum, ",",
                               .$Class, ",", 
                               .$Order, ",", 
                               .$Family, ",", 
                               .$Genus, ",",
                               .$Species)) %>% 
  select(Accession, OLD_NAME)
df_fungi_new <- df_fungi_updated %>% 
  filter(Changed == "YES") %>% 
  add_column(NEW_NAME = paste0(.$Kingdom, ",",
                               .$Phylum, ",",
                               .$Class, ",", 
                               .$Order, ",", 
                               .$Family, ",", 
                               .$Genus, ",",
                               .$Species)) %>% 
  select(Accession, NEW_NAME)

df_fungi_names_changed <- left_join(df_fungi_old, df_fungi_new, by = "Accession")

#write.xlsx(df_fungi_names_changed, "Name_changed_list.xlsx")

table(df_fungi_changed$Order)

df_fungi_no_change <- df_fungi_updated %>% filter(Changed == "NO") %>% 
  select(-Changed)

df_fungi_new <- bind_rows(df_fungi_changed, df_fungi_no_change)

new_df <- bind_rows(df_fungi_new, df_plant, df_other)
nrow(new_df) == nrow(df_split)

# Rebuild the SINTAX header
new_df$Header <- paste0(new_df$Accession, 
                        ";tax=k:", new_df$Kingdom,
                        ",p:", new_df$Phylum,
                        ",c:", new_df$Class,
                        ",o:", new_df$Order,
                        ",f:", new_df$Family,
                        ",g:", new_df$Genus,
                        ",s:", new_df$Species)

# Generate the FASTA formatted file

fasta_headers <- new_df %>% 
  drop_na() %>% 
  select(Header) %>% 
  rename(Seq = Header) %>% 
  add_column(X = seq(1:nrow(na.omit(new_df))))
fasta_seq <- new_df %>% 
  drop_na() %>% 
  select(Seq) %>% 
  add_column(X = seq(1:nrow(na.omit(new_df))))

fasta <- bind_rows(fasta_headers, fasta_seq) %>% 
  arrange(X) %>% 
  select(-X)

#write.table(fasta, "Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax_UPDATED_SYNONYMS_Feb_2026.fasta", quote = F, row.names = F, col.names = F)

table(df_fungi_updated$Changed)
