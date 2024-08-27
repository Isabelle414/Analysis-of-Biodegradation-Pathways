#' ---
#' title: Global Analysis of CAP Degradation Pathway Using GMGC
#' date: 27.08.2024
#' ---


#### Prepare the data from the GMGC search for analysis ####
# Check working directory
getwd()

# Install necessary packages
pacman::p_load("tidyverse", "readxl", "janitor", "ggpubr", "BiocManager", "ape", "ggtree","ggalluvial", "pheatmap", "grid")

# Read in the biodegradation gene data
key <- read_excel("data/cap_degradation_gene_metadata_extended.xlsx")

# Delete the fifth column of the biodegradation gene metadata (was just an empty column in excel)
key <- key[, -5]

# Add a column to the biodegradation gene metadata associated with the steps Prot1 - Prot17
keymod <- key %>%
  dplyr::mutate(protnum = paste0("Prot_", 1:17))
view(keymod)

# Read in all of the spreadsheets from GMGC at once
files <- list.files("data/gmgc/New_naming", pattern = "20240312_IT_DataGmgcProt", full.names = T)
files

readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = T))
  ) %>%
  unnest(.,cols = c(file_contents)) %>%
  janitor::clean_names()
view(readin)
colnames(readin)

# Create new columns in the GMGC data in order to merge with the biodegradation gene metadata
datfix <- readin %>%
  dplyr::mutate(genus = word(sep = " ", taxon_predicted, 1)) %>%
  dplyr::mutate(protnum_temp = word(sep = "DataGmgc", filename, 2)) %>%
  dplyr::mutate(protnum = word(sep = "_v1.csv", protnum_temp, 1)) # for merge

# Combine "datfix" with the biodegradation gene data
comb <- datfix %>%
  dplyr::group_by(protnum) %>%
  dplyr::left_join(., keymod, by = "protnum")

# Filter and extract genus-level taxonomic data
comb_genus <- comb %>%
  filter (str_detect(taxon_predicted, "environmental", negate = TRUE)) %>% # removes environmental
  filter (str_detect(taxon_predicted, "superkingdom", negate = TRUE)) %>% # removes superkingdoms
  filter (str_detect(taxon_predicted, "order", negate = TRUE)) %>% # removes orders
  filter (str_detect(taxon_predicted, "phylum", negate = TRUE)) %>% # removes phyla
  filter (str_detect(taxon_predicted, "class", negate = TRUE)) %>% # removes classes
  filter (str_detect(taxon_predicted, "family", negate = TRUE)) %>% # removes families
  filter (str_detect(taxon_predicted, "cellular", negate = TRUE)) %>% # removes cellular organisms (no rank)
  mutate(
    genus = case_when(
      grepl("Candidatus", taxon_predicted) ~ stringr::word(taxon_predicted, sep = " ",2),
      
      .default = stringr::word(taxon_predicted, sep = " ", 1)
    ),
    genus = gsub("\\[|\\]","", genus), # remove the brackets of some genus names
    reaction_name = gsub("_", " ", reaction_name) # write the reaction name with a space bar
  )  

# There are many cases where there are multiple hits per genus - we filter to take one per genus & step
comb_genus_unique <- comb_genus %>%  
  dplyr::group_by(protnum, genus) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(step_ordered = as.numeric(word(protnum, sep = "Prot_", 2))) %>%
  dplyr::mutate(first_mid_last = word(step, sep = " ", 1))

# Check if there are any taxa which share all pathway steps
checktaxa <- comb_genus_unique %>%
  dplyr::group_by(genus) %>%
  summarize(n_steps = n()) %>%
  dplyr::arrange(desc(n_steps))
view(checktaxa)

# Take a subset of the top forty taxa of checktaxa
checktaxa_40 <- checktaxa[1:40,]

# Export a text file with the top forty taxa
dat_40 <- checktaxa_40$genus
write.table(dat_40, row.names = FALSE,
            col.names = FALSE,
            "data/gmgc/genus_gmgc")

# Manually import a text file with the curated top genera sharing the most pathway steps
dat_top <- read.table("data/gmgc/genus_check_gmgc", quote="\"")

# Change the column name to genus
colnames(dat_top) <- c("genus")

# Merge "checktaxa_40" with the top genera sharing most pathway steps to evaluate 
# their contribution to CAP degradation
checktaxa_top <- checktaxa_40 %>%
  dplyr::right_join(.,dat_top, by = "genus")

# Filter "comb_genus_unique" to include only the top genera sharing the most pathway steps, 
# generating a dataframe with complete metadata for further analysis
genus_top <- comb_genus_unique %>%
  filter(genus %in% checktaxa_top$genus)


##### Pathway plot for global analysis of CAP degradation ####
# Create a plot
plot_gmgc <- ggplot(genus_top, aes(x = step_ordered , y = genus, fill = first_mid_last)) +
  geom_tile() + 
  xlab ("Protein") +
  ylab ("Genus") +
  scale_fill_manual(name = "Steps of Degradation",
                    values = c("#E69F00", "#56B4E9", "#009E73"),
                    breaks = c("first", "middle", "last"),
                    labels = c("First steps CAP degradation\n 1-3: CAP amide tail hydrolysis\n 4-5: CAP deglycosylation",
                               "Middle steps CAP degradation\n 6-8: Ring opening\n 9: Uracil reduction\n 10: Amide hydrolysis\n 11: Deamination\n 12: Acetyl-CoA synthase\n 13: Deamination\n 14: Keto reduction",
                               "Last steps CAP degradation\n 15-17: Defluorination")) +
  theme_bw() +
  scale_x_continuous(breaks = c(1:17)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        panel.grid.major.x = element_line(size = 0.5),  # keep major vertical grid lines
        panel.grid.minor.x = element_blank(),  # remove minor vertical grid lines
        panel.grid.major.y = element_line(size = 0.5), # keep major horizontal grid lines
        panel.grid.minor.y = element_blank()) # remove minor horizontal grid lines)
plot_gmgc

# Write to file
ggsave(plot_gmgc, filename = "output/20240626_IT_conservationplotgmgc_v5.png", height = 8, width = 10)


#### Phylogenetic tree with heatmap to illustrate protein relatedness ####
# Export the top genera sharing most pathway steps for creating a phylogenetic tree (with NCBI common tree)
dat_tree <- genus_top$genus
write.table(dat_tree, row.names = FALSE,
            col.names = FALSE,
            "data/for_trees/genus_top_v1.txt")

# Read in the tree for the top genera
phy <- read.tree("data/for_trees/genus_top_tree_v1.phy")
ggtree <- ggtree(phy) +
  geom_tiplab() +
  xlim(NA, 23) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))  # margin
ggtree

# Save the tree
ggsave(ggtree, filename = "output/20240628_IT_tree_v3.png", height = 6, width = 9)

# Extract the genera from "genus_top" into a vector
newlev <- read.table("data/for_trees/genus_top_ordered_v1.txt", stringsAsFactors = FALSE) %>%
  pull(V1)

# Reorder the "genus" column in "genus_top" using "newlev" to align heatmap and phylogenetic tree
genus_top_ordered <- genus_top %>%
  mutate(genus = factor(genus, levels = newlev)) %>%
  arrange(genus)

# Prepare data frame from "genus_top_ordered" for heatmap visualization
heat <- data.frame(table(genus_top_ordered$genus, genus_top_ordered$step_ordered))
heat_spread <- reshape2::dcast(heat, Var1 ~ Var2, value.var = "Freq", fill = 0)
heat_fix <- heat_spread %>%
  mutate(rownames = Var1) %>%
  column_to_rownames(var = "rownames") %>%
  select(-Var1)

# Show the internal structure of "heat_fix"
str(heat_fix)

# Have only elements of 0 and 1 in the data frame
heat_fix[heat_fix == 2 | heat_fix == 3] <- 1

# Create a heatmap with pheatmap
plot_heat_gmgc <- pheatmap(heat_fix,
                         labels_col = c("Prot 1", "Prot 2", "Prot 3", "Prot 4",
                                        "Prot 5", "Prot 6", "Prot 7", "Prot 8",
                                        "Prot 9", "Prot 10", "Prot 11", "Prot 12",
                                        "Prot 13", "Prot 14", "Prot 15", "Prot 16",
                                        "Prot 17"),
                         legend_breaks = c(0,1),
                         legend_labels = c("no hit","hit"),
                         col = c("#E69F00","#009E73"),
                         cluster_rows = FALSE,
                         cluster_col = FALSE,
                         angle_col = 90,
                         # grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16)),
                         # main = "Heatmap experimental proteins",
                         # dendrogram = "none",
                         gap_col = FALSE,
                         border_color = "black")

# Save heatmap                      
ggsave(plot_heat_gmgc, filename = "output/20240628_IT_heat_v1.png", height = 6, width = 4)


#### Alluvial plot for habitat distribution ####
# Create a data frame containing genus and habitat
data_alluvial <- data.frame(
  genus = genus_top$genus,
  habitat = genus_top$habitat)

# Modify "data_alluvial"
data_alluvial <- data_alluvial %>%
  filter (habitat != "-") %>% # delete all entries with unknown environment (-)
  separate_rows (habitat, sep = ",") # when a genus is found in several environments: add a new row for each habitat

# Check the data type
class(data_alluvial$habitat)
class(data_alluvial$genus)

# Generate an alluvial plot
plot_habitat <- ggplot(data = data_alluvial, aes(axis1 = genus, axis2 = habitat)) +
  geom_alluvium(aes(fill = genus)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  theme_void() + 
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = margin(1,12,1,1)) +
  labs(fill = "Host Genus of Protein")
plot_habitat

# Save the alluvial plot
ggsave(plot_habitat, filename = "output/20240628_IT_habitat_v4.png", height = 7, width = 8)
 file









