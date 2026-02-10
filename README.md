## READ ME FOR NAME SYNONYMY SCREENING OF UNITE ##


Files: 

Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax.fasta.zip - UNITE database combined with PLANiTS database. This has been cut with ITSx to only contain ITS2 barcodes. Additionally we have removed high ranking uncertain accessions with these identities from the original UNITE:
- Ascomycota_sp
- Basidiomycota_sp
- Fungi_phy_Incertae_sedis
- GS01_phy_Incertae_sedis
- Eukaryota_kgd_Incertae_sedis
- phy_Incertae_sedis_sp
- Viridiplantae (replaced with PLANiTS)


Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax_UPDATED_SYNONYMS_Feb_2026.fasta.zip - Post screening database. Can be replicated by running the database above through the R script 'Taxize_script_for_cleaning_names.R'. 
Method in short: species names from UNITE has been run through the package 'taxize' in R against the Index Fungorum database. Any names that are flagged as synonyms are given a full updated taxonomy string based on in one of two ways; 
1 - If the genera is represented in UNITE already, the updated names are given the same parent taxonomy. This is done for the majority of synonymous taxa.
2 - If the genus name is new to UNITE, the parent taxonomy is assigned through taxize::classification() using full names from 'gbif' 

Some groups within Saccharomycetales has been manually changed as for example Candida was classed as Saccharomycetales_Incertae_sedis while it should according to Index Fungorum be within Serinales -> Debaryomycetaceae -> Candida. Serinales and Saccharomycetales in general has been curated.

Additionally, a check for proper dichotomy in the taxonomy strings has been made and a few instances of Families and Genera with more than one parent has been fixed.

Changed accessions has had a "_c" pasted to the end of the UNITE accession

In total; 6 710 accessions has been updated, 130 278 is unchanged. 


Taxize_script_for_cleaning_names.R - R script that takes the unscreened fasta file as input and outputs the screened file. No external annotation beyond this has been made.


Name_changed_list.xlsx - List of all the OLD_NAME -> NEW_NAME changes from the database. Original accession numbers are included. 


UPDATE 10.Feb.2026:

Some name-changes flagged as erraneous and some non-dichotomour branching fixed due to same names for genera in plants/animals and fungi

Files have been updated accordingly
