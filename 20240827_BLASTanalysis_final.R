#' ---
#' title: Global Analysis of CAP Degradation Pathway Using Protein BLAST
#' date: 27.08.2024
#' ---


#### Prepare the data from the Protein BLAST search for analysis ####
# Check your working directory
getwd()

# Next load/install necessary packages
pacman::p_load("tidyverse", "readxl", "janitor", "ggpubr")

# Read in the biodegradation gene metadata as gene key
key <- read_excel("data/cap_degradation_gene_metadata_extended.xlsx") %>%
  dplyr::select(-5) # delete the fifth column

# Add a column associated with the steps Prot1 to Prot17 to the gene key
keymod <- key %>%
  dplyr::mutate(protnum = paste0("Prot_", 1:17))

# Read in all of the spreadsheets from the BLAST search at once
files <- list.files("data/Blast/New_numbering", pattern = "20240408_IT_DataBlastProt_", full.names = T)

readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_delim(delim = ";", file.path(.), col_names = T, col_types = c("c", "c", "c", rep("c", times = 9))))
  ) %>%
  unnest(.,cols = c(file_contents)) %>%
  janitor::clean_names()
view(readin)

# Create new columns in the BLAST data in order to merge with the gene key
datfix <- readin %>%
  dplyr::mutate(genus = word(sep = " ", scientific_name, 1)) %>%
  dplyr::mutate(protnum_temp = word(sep = "DataBlast", filename, 2)) %>%
  dplyr::mutate(protnum = word(sep = "_v1.csv", protnum_temp, 1))

# Combine the BLAST data ("datfix") with the gene key ("keymod") by a left-join
comb <- datfix %>%
  dplyr::group_by(protnum) %>% 
  dplyr::left_join(., keymod, by = "protnum") %>% 
  filter (!is.na(genus))

# (b) Prepare a new column with taxa information, including genera and other taxonomic levels
comb_genus <- comb %>%
  mutate(
    genus = case_when(
      grepl("ncultured|nclassified|Candidatus", scientific_name) ~ stringr::word(scientific_name, sep = " ",2),
      .default = stringr::word(scientific_name, sep = " ", 1)
    ),
    genus = gsub("\\[|\\]", "", genus)
  )

# There are many cases where there are multiple hits per taxa - we filter to take one per taxa & step
comb_genus_unique <- comb_genus %>%  
  dplyr::group_by(protnum, genus) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(step_ordered = as.numeric(word(protnum, sep = "Prot_", 2))) %>%
  dplyr::mutate(first_mid_last = word(step, sep = " ", 1))

# Check whether there are any taxa that encode all proteins essential for CAP degradation
checktaxa <- comb_genus_unique %>%
  dplyr::group_by(genus) %>%
  dplyr::summarize(n_steps = n()) %>%
  dplyr::arrange(desc(n_steps))
view(checktaxa)

# Make a subset of the top forty taxa sharing most pathway steps
checktaxa_40 <- checktaxa[1:40,]

# Export a text file with these top forty taxa
data_tree <- checktaxa_40$genus
write.table(data_tree, row.names = FALSE,
            col.names = FALSE,
            "data/for_trees/genus_blast")

# Import a text file with the manually curated top genera
dat_top <- read.table("data/for_trees/genus_check_blast", quote="\"") %>%
  rename(genus = V1)

# Merge "dat_top" with "checktaxa_40" to assess the enzymatic contribution of the 
# top genera to CAP degradation
checktaxa_top <- checktaxa_40 %>%
  dplyr::right_join(.,dat_top, by = "genus")

# Filter "comb_genus_unique" to include only the top genera sharing the most pathway steps,
# generating a data frame with metadata for further analysis
genus_top <- comb_genus_unique %>%
  filter(genus %in% checktaxa_top$genus)


#### Pathway plot ####
# Create a pathway plot
plot_genus_top <- ggplot(genus_top, aes(x = step_ordered , y = genus, fill = first_mid_last)) +
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
        panel.grid.major.x = element_line(size = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.y = element_blank())
plot_genus_top

# Write to file
ggsave(plot_genus_top, filename = "output/20240625_IT_checktaxa_top_blast_v1.png", height = 8, width = 10)