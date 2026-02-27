## Name synonymy screening of UNITE ##


Files: 

### "Concat_UNITE_v2.fasta.zip"  
Curated ITS2 UNITE database with PLANiTS added as a replacement for UNITE Viridiplantae sequences. Ready for use with the SINTAX taxonomic assignment algorithm


#### Method description for curation: 

"Old_database/Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax.fasta.zip" 
UNITE database combined with PLANiTS database. This has been cut with ITSx to only contain ITS2 barcodes. Additionally we have removed high ranking uncertain accessions with these identities from the original UNITE:
- Ascomycota_sp
- Basidiomycota_sp
- Fungi_phy_Incertae_sedis
- GS01_phy_Incertae_sedis
- Eukaryota_kgd_Incertae_sedis
- phy_Incertae_sedis_sp
- Viridiplantae (replaced with PLANiTS)


"Old_database/Concat_UNITE_ITS2_PLANiTS_22.10.2025_for_sintax_UPDATED_SYNONYMS_Feb_2026.fasta.zip"
Post screening database. UPDATED 02. Feb. 2026.
Can be replicated by running the database above through the R script 'Taxize_script_for_cleaning_names.R'. 
Method in short: species names from UNITE has been run through the package 'taxize' in R against the Index Fungorum database. Any names that are flagged as synonyms are given a full updated taxonomy string based on in one of two ways; 
1 - If the genera is represented in UNITE already, the updated names are given the same parent taxonomy. This is done for the majority of synonymous taxa.
2 - If the genus name is new to UNITE, the parent taxonomy is assigned through taxize::classification() using full names from 'gbif' 

Some groups within Saccharomycetales has been manually changed as for example Candida was classed as Saccharomycetales_Incertae_sedis while it should according to Index Fungorum be within Serinales -> Debaryomycetaceae -> Candida. Serinales and Saccharomycetales in general has been curated.

Additionally, a check for proper dichotomy in the taxonomy strings has been made and a few instances of Families and Genera with more than one parent has been fixed.

Changed accessions has had a "_c" pasted to the end of the UNITE accession

In total; 6 710 accessions has been updated, 130 278 is unchanged. 


"Taxize_script_for_cleaning_names.R" 
R script that takes the unscreened fasta file as input and outputs the screened file. No external annotation beyond this has been made.


"Name_changed_list.xlsx" 
List of all the OLD_NAME -> NEW_NAME changes from the database. Original accession numbers are included. 


UPDATE 10.Feb.2026:

Some name-changes flagged as erroneous and some non-dichotomous branching fixed due to same names for genera in plants/animals and fungi

Files have been updated accordingly


UPDATE 26.Feb.2026

The original and screened database contains many reference sequences with different accessions that are very similar (99 to 100% sequence similarity) which have varying quality of the annotation. For example SH0126933.10FU with the name Tomentella_sp is 99.1% similar to SH0917554.10FU Tomentella_subtestacea. To avoid problems with many genera represented by large amounts of Genus_sp annotation with a 99% or higher sequence similarity to a better resolved accession we did a BLAST match of all sequences to each other in the database and extracted all the pairwise similarities of 99% or higher. Results are shown in Concat_UNITE_sim99.xlsx. Additionally, the accession number has been added to the end of the fasta headers so that the final annotations for each query can be tracked to an exact reference id the database. This also allows SINTAX to converge on a single unique name for each query regardless (instead of reporting a hit against one of several Tomentella_sp, with no way of tracking it back to which one). 

All the sequences that are suggested removed from a matching pair has been removed and the cleaned final database is named Concat_UNITE_v2.fasta. 
