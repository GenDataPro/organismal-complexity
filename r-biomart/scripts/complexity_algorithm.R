#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

# Load Biomart
library(biomaRt)
ensembl <- useMart("ensembl")
listDatasets(ensembl)
# Set the species vector 
species <- c("hsapiens", "mmusculus", "ggallus", "celegans", "dmelanogaster", "cintestinalis", 	
             "trubripes", "xtropicalis", "mmulatta")
# Make a connection to ensembl for all species
ensembl <- sapply(species, function(s) useMart("ensembl", paste0(s, "_gene_ensembl")))
names(ensembl) <- c(species)
# Get the human genes
hsapien_PC_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"), 
                          filters = "biotype", 
                          values = "protein_coding", 
                          mart = ensembl[["hsapiens"]])

hsapien_PC_genes <- hsapien_PC_genes[hsapien_PC_genes$chromosome_name %in% c(1:23, "X", "Y"), ]
ensembl_gene_ID <- hsapien_PC_genes$ensembl_gene_id

write.csv(hsapien_PC_genes, file = "Human_Protein_Coding_Genes_all.csv")

# Get the orthologues
# Exclude humans (these have already been retrieved) by usings species[-1]
all_homologues <- list()
all_homologues <- lapply(species[-1], function(s) getBM(attributes = c("ensembl_gene_id", 
                                                                       "external_gene_name", 
                                                                       paste0(s, c("_homolog_ensembl_gene",
                                                                                   "_homolog_associated_gene_name"))),
                                                        filters = "ensembl_gene_id",
                                                        values = c(ensembl_gene_ID),
                                                        mart = ensembl[["hsapiens"]]))
names(all_homologues) <- c(species[-1])

View(all_homologues$mmusculus)
# Fix blank homologues
homologue_fix <- list()
homologue_fix <- lapply(species[-1], function(s) all_homologues[[paste0(s)]][!(is.na(all_homologues[[paste0(s)]][paste0(s, "_homolog_ensembl_gene")])
                                                                                | all_homologues[[paste0(s)]][paste0(s, "_homolog_ensembl_gene")]==""), ])
all_homologues <- homologue_fix
names(all_homologues) <- c(species[-1])
# Remove Alternative 
# homologues with 1-many causes problems, have to select genes with highest target ID
# Extract Target IDs
target_id <- list()
target_id <- lapply(species[-1], function(s) getBM(attributes = c("ensembl_gene_id", 
                                                                  "external_gene_name", 
                                                                  "hsapiens_homolog_associated_gene_name", 
                                                                  "hsapiens_homolog_perc_id"), 
                                                   filters = "ensembl_gene_id", 
                                                   values = all_homologues[[paste0(s)]][paste0(s, "_homolog_ensembl_gene")], 
                                                   mart = ensembl[[paste0(s)]]))
names(target_id) <- c(species[-1])
# Select top target IDs for all species (excl hsapiens)
library(plyr)
top_target_id <- list()
top_target_id <- lapply(species[-1], function(s) ddply(target_id[[paste0(s)]],
                                                       .(hsapiens_homolog_associated_gene_name), 
                                                       function(x)x[which.max(x$hsapiens_homolog_perc_id), ]))
names(top_target_id) <- c(species[-1])
# put selected homologues into dataframe:
homologues <- list()
homologues <- lapply(species[-1], function (s) data.frame(Ensembl_Gene_ID = top_target_id[[paste0(s)]]$ensembl_gene_id,
                                                          Human_Gene_Name = top_target_id[[paste0(s)]]$hsapiens_homolog_associated_gene_name,
                                                          Gene_Name = top_target_id[[paste0(s)]]$external_gene_name))
names(homologues) <- c(species[-1])
all_homologues [["hsapiens"]] <- hsapien_PC_genes
homologues[["hsapiens"]] <- data.frame(Ensembl_Gene_ID = all_homologues[["hsapiens"]]$ensembl_gene_id,
                                       Human_Gene_Name = all_homologues[["hsapiens"]]$external_gene_name,
                                       Gene_Name = all_homologues[["hsapiens"]]$external_gene_name)
# Get the paralogues for all species using top target ID
paralogues <- list()
paralogues <- lapply(species, function(s) getBM(attributes = c("ensembl_gene_id",
                                                               paste0(s, "_paralog_ensembl_gene")), 
                                                filters = "ensembl_gene_id", 
                                                values = homologues[[paste0(s)]]["Ensembl_Gene_ID"] , mart = ensembl[[paste0(s)]]))
names(paralogues) <- c(species)
# Remove genes with no human gene name
library(dplyr) 
paralogue_fix <- list()
paralogue_fix <- lapply(species, function(s) paralogues[[paste0(s)]][!(paralogues[[paste0(s)]][paste0(s, "_paralog_ensembl_gene")]==""), ])
paralogues <- paralogue_fix
names(paralogues) <- c(species)
# Count paralogues
count_paralogues <- list()
count_paralogues <- lapply(species, function(s) as.data.frame(table(paralogues[[paste0(s)]]$ensembl_gene_id)))
names(count_paralogues) <- c(species)
paralogue_names <- c("Ensembl_Gene_ID", "Number_of_paralogues")
count_paralogues <- lapply(count_paralogues, setNames, paralogue_names)
# Get domains
domains <- list()
domains <- lapply(species, function(s) getBM(attributes = c("ensembl_gene_id",
                                                            "pfscan_start", "ensembl_transcript_id", "transcript_biotype"), 
                                             filters = "ensembl_gene_id", 
                                             values = homologues[[paste0(s)]]["Ensembl_Gene_ID"] , mart = ensembl[[paste0(s)]]))
# Subset protein coding transcripts 
names(domains) <- c(species)
domain_subset <- list()
domain_subset  <- lapply(species, function(s) subset(domains[[paste0(s)]], transcript_biotype=="protein_coding"))
names(domain_subset) <- c(species)
# Problem with biomart, some genes have blank domain entries 
# e.g. ENSG00000136044
# FIX: Remove domain blanks  
domain_fix <- list()
domain_fix <- lapply(species, function(s) domain_subset[[paste0(s)]][!(is.na(domain_subset[[paste0(s)]]$pfscan_start)
                                                                 | domain_subset[[paste0(s)]]$pfscan_start==""), ])
domains <- domain_fix
names(domains) <- c(species)
# Count the domains
count_domains <- list()
count_domains <- lapply(species, function(s) as.data.frame(table(domains[[paste0(s)]]$ensembl_gene_id)))
names(count_domains) <- c(species)
domain_names <- c("Ensembl_Gene_ID", "Number_of_Domains")
count_domains <- lapply(count_domains, setNames, domain_names)
# Get isoforms
isoforms <- list()
isoforms <- lapply(species, function(s) getBM(attributes = c("ensembl_gene_id", 
                                                             "ensembl_transcript_id",
                                                             "transcript_biotype", "chromosome_name"), 
                                              filters = "ensembl_gene_id", 
                                              values = homologues[[paste0(s)]]["Ensembl_Gene_ID"] , mart = ensembl[[paste0(s)]]))
names(isoforms) <- c(species)
# Subset protein coding genes and specific chromosome genes
isoform_subset <- list()
isoform_subset  <- lapply(species, function(s) subset(isoforms[[paste0(s)]], transcript_biotype=="protein_coding"))
names(isoform_subset) <- c(species)
# VARS2
# Count isoforms
count_isoforms <- list()
count_isoforms <- lapply(species, function(s) as.data.frame(table(isoforms[[paste0(s)]]$ensembl_gene_id)))
names(count_isoforms) <- c(species)
isoform_names <- c("Ensembl_Gene_ID", "Number_of_Isoforms")
count_isoforms <- lapply(count_isoforms, setNames, isoform_names)
## Make species specific master table 
master_table <- list()
master_table <- lapply(species, function (s) data.frame(Ensembl_Gene_ID = homologues[[paste0(s)]]$Ensembl_Gene_ID,
                                                        Human_Gene_Name = homologues[[paste0(s)]]$Human_Gene_Name,
                                                        Gene_Name = homologues[[paste0(s)]]$Gene_Name, 
                                                        Species = paste0(s)))
names(master_table) <- c(species)
master_table <- lapply(species, function(s) merge(count_paralogues[[paste0(s)]], master_table[[paste0(s)]], by = c("Ensembl_Gene_ID"), all.y = TRUE)[, union(names(master_table[[paste0(s)]]), names(count_paralogues[[paste0(s)]]))])
names(master_table) <- c(species)
master_table <- lapply(species, function(s) merge(count_domains[[paste0(s)]], master_table[[paste0(s)]], by = c("Ensembl_Gene_ID"), all.y = TRUE)[, union(names(master_table[[paste0(s)]]), names(count_domains[[paste0(s)]]))])
names(master_table) <- c(species)
master_table <- lapply(species, function(s) merge(count_isoforms[[paste0(s)]], master_table[[paste0(s)]], by = c("Ensembl_Gene_ID"), all.y = TRUE)[, union(names(master_table[[paste0(s)]]), names(count_isoforms[[paste0(s)]]))])
names(master_table) <- c(species)

# make all species table
master_table[["all_species"]] <- rbind(master_table[["hsapiens"]], master_table[["mmusculus"]], master_table[["ggallus"]], master_table[["celegans"]], master_table[["dmelanogaster"]], master_table[["cintestinalis"]], master_table[["trubripes"]], master_table[["xtropicalis"]], master_table[["mmulatta"]])
# Change NA to 0
master_table[["all_species"]][is.na(master_table[["all_species"]])] <- 0
# Add 1 to each paralogue so that the gene is counted as a paralogue
master_table[["all_species"]]$Number_of_paralogues <- master_table$all_species[, "Number_of_paralogues"] + 1
                       
# Work out complexity score
complexity_score <- master_table[["all_species"]] # Make complexity score dataframe 
complexity_score$log_paralogues <- log2(complexity_score[, 5]) # log paralogues
complexity_score$log_domains <- log2(complexity_score[, 6])    # log domains
complexity_score$log_isoforms <- log2(complexity_score[, 7])   # log domains
complexity_score[complexity_score=="-Inf"] <- 0 # make log(0) values 0
complexity_score$Functional_Diversity <- rowSums(complexity_score[, c(8:10)]) # Add log values to calculate DF
# Write complexity score 

# change table layout
   # Subset by species
complexity_species <- list()
complexity_species <- lapply(species, function(s) subset(complexity_score, Species==paste0(s)))
names(complexity_species) <- c(species)
  # Subset only Human gene name and DF
species_df <- list()
species_df <- lapply(species, function(s) complexity_species[[paste0(s)]][, c("Human_Gene_Name", "Functional_Diversity")])
names(species_df) <- c(species)
# Name columns
colnames(species_df[["hsapiens"]]) <- c("Human_Gene_Name", "hsapiens")
colnames(species_df[["mmusculus"]]) <- c("Human_Gene_Name", "mmusculus")
colnames(species_df[["ggallus"]]) <- c("Human_Gene_Name", "ggallus")
colnames(species_df[["celegans"]]) <- c("Human_Gene_Name", "celegans")
colnames(species_df[["dmelanogaster"]]) <- c("Human_Gene_Name", "dmelanogaster")
colnames(species_df[["cintestinalis"]]) <- c("Human_Gene_Name", "cintestinalis")
colnames(species_df[["trubripes"]]) <- c("Human_Gene_Name", "trubripes")
colnames(species_df[["xtropicalis"]]) <- c("Human_Gene_Name", "xtropicalis")
colnames(species_df[["mmulatta"]]) <- c("Human_Gene_Name", "mmulatta")
# Merge individual species together

# 

complexity_spread <- merge(species_df[["hsapiens"]], species_df[["mmusculus"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["ggallus"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["celegans"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["dmelanogaster"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["cintestinalis"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["trubripes"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["xtropicalis"]], by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, species_df[["mmulatta"]], by = c("Human_Gene_Name"), all.x = TRUE)


# Pearsons Correlation 
ctn_df <- read.table(text = "Species 	CTN
                     hsapiens	171
                     mmusculus	157
                     ggallus	150
                     celegans	28.5
                     dmelanogaster	60
                     cintestinalis	71
                     trubripes	114
                     xtropicalis	121
                     mmulatta	171
                     ", header = T)  
# How many unique genes should there be:
# unique(complexity_score[, c("Human_Gene_Name")]) # answer =  20,263 
library(tidyverse) 
functional_diversity <- complexity_score[, c("Human_Gene_Name", "Species", "Functional_Diversity")]
gene_corr <- functional_diversity %>%
  left_join(ctn_df) %>%
  group_by(Human_Gene_Name) %>%
  summarise(cor(Functional_Diversity, CTN)) 
# Add to correlation to converted table
complexity_spread <- merge(complexity_spread, gene_corr, by = c("Human_Gene_Name"), all.x = TRUE)
complexity_spread <- merge(complexity_spread, complexity_species[["hsapiens"]][, c("Human_Gene_Name", "Ensembl_Gene_ID")], all.x = TRUE)
complexity_spread <- complexity_spread[,c("Ensembl_Gene_ID","Human_Gene_Name", "celegans", "dmelanogaster", "cintestinalis", "trubripes", "xtropicalis", "ggallus", "mmusculus","mmulatta", "hsapiens", "cor(Functional_Diversity, CTN)")]




