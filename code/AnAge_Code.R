###### Script for the manuscript Morimoto and Pietras (2025) 
###### NOte*: The analysis are not presented in the order of the manuscript, but they have been signposted here accordingly.
set.seed(12310)
## packages
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)
library(seqinr)
library(ape)
library(phangorn)
library(tidyr)
library(purrr)
library(BiocManager)
library(Biostrings)
library(patchwork)
library(phytools)
library(picante)
library(lmerTest)
library(PCAtest)
## AMino acid datasets
amino_acids_codon <- data.frame(aa_letter = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
                                numb_codons = c(4, 6, 2, 3, 2, 2, 2, 4, 2, 3, 6, 2, 1, 2, 4, 6, 4, 1, 2, 4))

### costs
# Define the data
amino_acids_costs <- data.frame(aa_letter =  c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), 
                                costATP = c(11.7, 24.7, 12.7, 15.3, 52, 11.7, 38.3, 32.3, 30.3, 27.3, 34.3, 14.7, 20.3, 16.3, 27.3, 11.7, 18.7, 23.3, 74.3, 50),
                                decay_invtime = c(1, 30, 9, 5, 4, 1, 14, 2, 8, 2, 13, 10, 3, 8, 4, 6, 6, 2, 12, 7),
                                cost_ATPtime = c(12, 741, 114, 77, 208, 12, 536, 65, 242, 55, 446, 147, 61, 130, 109, 70, 112, 47, 892, 350)) %>%
  right_join(., amino_acids_codon) %>%
  mutate(cost_ATPtime_codon = cost_ATPtime/numb_codons)



### amino acid information
# Your existing data frame
# Creating a single data frame with all columns, including full amino acid names
aa_data <- data.frame(
  'aa_letter' = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                  'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'),
  'aa_name' = c('Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid', 'Phenylalanine', 'Glycine', 'Histidine', 'Isoleucine', 
                'Lysine', 'Leucine', 'Methionine', 'Asparagine', 'Proline', 'Glutamine', 'Arginine', 'Serine', 'Threonine', 
                'Valine', 'Tryptophan', 'Tyrosine'),
  'mol_weight' = c(89.1, 121.2, 133.1, 147.1, 165.2, 75.1, 155.2, 131.2, 146.2, 131.2, 149.2, 132.1, 115.1, 146.2, 174.2, 105.1, 119.1, 117.1, 204.2, 181.2))


# Create data frame
amino_acids <- right_join(aa_data, amino_acids_codon, by = join_by('aa_letter')) %>%
  right_join(., amino_acids_costs)
amino_acids 




















### Analysis 1: Testing if amino acid profiles vary
## Output path
exome_processed_path <- "..."


### Storing all file names
exome_all <- list.files(path = exome_processed_path, 
                        pattern="\\.csv$")

### Listing all files in the folder and pasting the file name
exome_datasets_all <- lapply(paste(exome_processed_path ,
                                   exome_all,
                                   sep = "/"), 
                             read.csv, 
                             header = TRUE)

## joining aa information
exome_all_full <- bind_rows(exome_datasets_all) %>%
  select(-X.1, -X) %>%
  right_join(., amino_acids, by = join_by("aa_letter")) %>%
  mutate()


# getting the total sum of aa and proportions per species
total_count_overall_sp <- exome_all_full %>%
  group_by(species) %>%
  summarise(total_aa_overall_sp = sum(aa_total_count)) %>%
  left_join(., exome_all_full) %>%
  mutate(relativfreq_aa_sp = aa_total_count/total_aa_overall_sp,
         relativfreq_aa_sp_codon = relativfreq_aa_sp/numb_codons)### HERE





### PCA by class

### AnAge dataset

anage <- read.csv("...",
                  header = TRUE, stringsAsFactors = TRUE, na.strings = "")


### only using "acceptable" data quality
anage_acceptable <- subset(anage, Data.quality == "acceptable" | Data.quality == "high")

## creating new species
anage_acceptable$species <- with(anage_acceptable, paste(Genus, Species))


## how many species in both datasets
table(anage_acceptable$species %in% total_count_overall_sp$species)

## ensuring that the datais clean
anage_acceptable_subset <- anage_acceptable[anage_acceptable$species %in% total_count_overall_sp$species,]
anage_acceptable_subset_clean <- subset(anage_acceptable_subset, !is.na(Maximum.longevity..yrs.)) ## no entries where this column is NA, which is good


## merging with aa proportions
combined_anage <- merge(total_count_overall_sp, anage_acceptable_subset_clean)


## principal component analysis
prcomp_data_anage <-  combined_anage %>%
  select(lifespan = Maximum.longevity..yrs.,
         weight = Adult.weight..g.,
         class,
         order, 
         family,
         genus = Genus,
         species,
         aa_name,
         total_aa_overall_sp,
         aa_total_count,
         numb_codons) %>%
  group_by(species, class, order, family, genus, aa_name, numb_codons, lifespan, weight, proteome_size = total_aa_overall_sp) %>%
  summarise(aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons) %>% ## standardising by codons
  ungroup() %>%
  select(-numb_codons) %>%
  pivot_wider(., 
              names_from = aa_name,
              values_from = aa_freq_codon) %>%
  filter(class != "Saccharomycetes",
         class != "Hyperoartia") %>%
  na.omit()


prcomp_data_anage$class <- ifelse(prcomp_data_anage$class == "Lepidosauria", "Reptillia", paste(prcomp_data_anage$class))

## PCA on amino acids
pccomp_anage_vars <- prcomp(prcomp_data_anage[9:28], scale. = FALSE)
#PCA_test_anage <- PCAtest(prcomp_data_anage[9:28]) ## PCA test
summary(pccomp_anage_vars)
#PCA_test_anage <- PCAtest(prcomp_data_anage[9:28])
#summary(PCA_test_anage)

## data for plotting
pccomp_anage_df <- data.frame(prcomp_data_anage[1:8],
                              as.data.frame(pccomp_anage_vars$x))
pccomp_anage_df %>%
  group_by(class) %>%
  summarise(n())

## plot
PC13_plot <- ggplot(pccomp_anage_df, aes(PC1, PC3, color=class)) +
  geom_point(size=2, alpha = 0.5) +
  stat_ellipse() + 
  theme_light() + 
  theme(panel.grid = element_blank())+ 
  scale_color_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                     values = c('steelblue2', 'goldenrod1', 'firebrick1',
                                'chartreuse2', "grey20", "purple"))

PC12_plot <-ggplot(pccomp_anage_df, aes(PC1, PC2, color=class)) +
  geom_point(size=2, alpha = 0.5) + 
  stat_ellipse() + 
  theme_light() + 
  theme(panel.grid = element_blank())+ 
  scale_color_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                     values = c('steelblue2', 'goldenrod1', 'firebrick1',
                                'chartreuse2', "grey20", "purple"))

PC12_plot + PC13_plot + plot_layout(guides="collect")



## Permanova
adonis2(
  pccomp_anage_df[ , c("PC1","PC2", "PC3")] ~ class,
  data = pccomp_anage_df,
  method = "euc"
)


### Revisions

#Reviewer 1: non-linear relationship between 


#Pplots
Comparison_linearity <- ggplot(pccomp_anage_df, aes(x = log(weight^0.5), y = log(lifespan))) +
  geom_point(size = 3, pch = 21, alpha = 8, fill = "grey80") +
  geom_smooth(method = 'lm', se = TRUE, col = "black", fill = "black") + 
  geom_smooth(method = 'gam', se = TRUE, col = "red", fill = "red") + 
  theme_linedraw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))+ 
  xlab("log(Weight^0.5)") +
  ylab("log(Lifespan)")
Comparison_linearity


### groups datasets
combined_anage_Mammals <- combined_anage %>%
  filter(Class == "Mammalia")


combined_anage_Fish <- combined_anage %>%
  filter(Class == "Teleostei" | Class == "Holostei" | Class == "Actinopterygii" | Class == "Chondrostei")

combined_anage_SharksRay <- combined_anage %>%
  filter(Class == "Chondrichthyes")

combined_anage_Birds <- combined_anage %>%
  filter(Class == "Aves")

combined_anage_Amphibia <- combined_anage %>%
  filter(Class == "Amphibia")

combined_anage_Reptiles <- combined_anage %>%
  filter(Class == "Reptilia")



## loading trees

### family level tree: https://zenodo.org/records/13909811
tettree_Mammals <- read.tree("...mammals.tre")
avgtettree_Mammals <- ls.consensus(tettree_Mammals)
plot(avgtettree_Mammals)

tettree_Fish <- read.tree("...fish.tre")
plot(tettree_Fish)

tettree_Birds <- read.tree("...birds.tre")
avgtettree_Birds <- ls.consensus(tettree_Birds)
plot(avgtettree_Birds)

tettree_Reptile <- read.tree("...reptiles.tre")
avgtettree_Reptile <- ls.consensus(tettree_Reptile)
plot(avgtettree_Reptile )

tettree_Amphibia <- read.tree("...amphibian.tre")
avgtettree_Amphibia <- ls.consensus(tettree_Amphibia)
plot(avgtettree_Amphibia)

tettree_Sharks <- read.tree("...sharkray.tre")
avgtettree_Sharks <- ls.consensus(tettree_Sharks)
plot(avgtettree_Sharks)


## getting residual lifespan
## getting residual lifespan
combined_anage_Mammals_short <- combined_anage_Mammals %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs.,
         weight = Adult.weight..g., species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit() %>%
  droplevels()


## removing asterisks
#avgtettree_Mammals$tip.label <- str_replace(avgtettree_Mammals$tip.label, "\\*", "")
table(avgtettree_Mammals$tip.label %in% combined_anage_Mammals_short$species)

##cleaning tree
avgtettree_Mammals_clean <- drop.tip(avgtettree_Mammals, avgtettree_Mammals$tip.label[!avgtettree_Mammals$tip.label %in% combined_anage_Mammals_short$species])
avgtettree_Mammals_clean$tip.label %in% unique(combined_anage_Mammals_short$species)

## dropping tips
combined_anage_Mammals_short2 <- combined_anage_Mammals_short[combined_anage_Mammals_short$species %in% avgtettree_Mammals_clean$tip.label,]

length(unique(combined_anage_Mammals_short2$species))

## controlling for body weight
mammals_model_LSBM <- nlme::gls(log(lifespan) ~ log(weight^0.5), 
                                data = combined_anage_Mammals_short2, 
                                correlation=corBrownian(1, avgtettree_Mammals_clean, form = ~species))

summary(mammals_model_LSBM)


## getting residual lifespan
combined_anage_Mammals_short2$residual_lifespan <- as.numeric(residuals(mammals_model_LSBM))





## checking 
mammals_model_LSBM_check <- nlme::gls(residual_lifespan ~ log(weight^0.5), 
                                      data = combined_anage_Mammals_short2, 
                                      correlation=corBrownian(1, avgtettree_Mammals_clean, form = ~species))
summary(mammals_model_LSBM_check)

## proteome size
mammals_model_LSBM_proteome <- nlme::gls(residual_lifespan ~ log(proteome_size), 
                                         data = combined_anage_Mammals_short2, 
                                         correlation=corBrownian(1, avgtettree_Mammals_clean, form = ~species))

summary(mammals_model_LSBM_proteome)
combined_anage_Mammals_short2$pred_residual_lifespan <- as.numeric(predict(mammals_model_LSBM_proteome, combined_anage_Mammals_short2))

## Mammals per amino acid
### analysis of individual AAs (called "RAW" as opposed to the PCA transformed: but see manuscript for detailed methods of controlling for confounding variables)
mammals_wide_aa <- combined_anage_Mammals_short2 %>% 
  select(Class, species, aa_full, aa_freq) %>%
  pivot_wider(., names_from = aa_full,
              values_from = aa_freq)

## model raw and confit
models_mammals_aa <- combined_anage_Mammals_short2  %>%
  nest(data = -aa_full) %>%
  mutate(mdout = map_df(data, function(.x){
    md <- nlme::gls(residual_lifespan ~ scale(aa_freq), data=.x, corBrownian(1, avgtettree_Mammals, form = ~species))
    out <- data.frame(lwr = as.vector(confint(md, level = 0.9975)[2,1]),
                      slope = as.vector(md$coefficients[2]),
                      upr = as.vector(confint(md, level = 0.9975)[2,2]))
    return(out)
  })) %>%
  unnest(mdout) %>%
  select(-data)



models_mammals_aa %>% arrange(aa_full)

## plot birds
ggplot(combined_anage_Mammals_short2, mapping = aes(x = scale(aa_freq_codon), y = residual_lifespan,
                                                    col = aa_full, fill = aa_full)) + 
  geom_point() + 
  geom_smooth(combined_anage_Mammals_short2, mapping = aes(x = scale(aa_freq_codon), y = residual_lifespan,
                                                           col = aa_full, fill = aa_full),
              method = 'lm') + 
  facet_wrap(~aa_full, ncol = 4, scales = "free") + 
  scale_fill_viridis_d('Amino acid') + 
  scale_color_viridis_d('Amino acid') + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "grey90")) + 
  labs(x = "Frequency", y = "Residual lifespan")


#
unique(combined_anage_Mammals_short2$species)
unique(combined_anage_Birds_short$species)
#








### * Birds 

## getting residual lifespan
combined_anage_Birds_short <- combined_anage_Birds %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs.,
         weight = Adult.weight..g., species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit()

#avgtettree_Mammals$tip.label <- str_replace(avgtettree_Mammals$tip.label, "\\*", "")
##cleaning tree
## removing asterisks
avgtettree_Birds$tip.label <- str_replace(avgtettree_Birds$tip.label, "\\*", "")
avgtettree_Birds_clean <- drop.tip(avgtettree_Birds, avgtettree_Birds$tip.label[!avgtettree_Birds$tip.label %in% combined_anage_Birds_short$species])



avgtettree_Birds_clean$tip.label %in% unique(combined_anage_Birds_short$species)

## dropping tips
combined_anage_Birds_short2 <- combined_anage_Birds_short[combined_anage_Birds_short$species %in% avgtettree_Birds_clean$tip.label,]


length(avgtettree_Birds_clean$tip.label)


## controlling for body weight
birds_model_LSBM <- nlme::gls(log(lifespan) ~ log(weight^0.5), 
                              data = combined_anage_Birds_short2, 
                              correlation=corBrownian(1, avgtettree_Birds_clean, form = ~species))

summary(birds_model_LSBM)

## getting resicual lifespan
combined_anage_Birds_short2$residual_lifespan <- as.numeric(residuals(birds_model_LSBM))


## checking 
birds_model_LSBM_check <- nlme::gls(residual_lifespan ~ log(weight^0.5), 
                                    data = combined_anage_Birds_short2, 
                                    correlation=corBrownian(1, avgtettree_Birds_clean, form = ~species))
summary(birds_model_LSBM_check)

## proteome size
birds_model_LSBM_proteome <- nlme::gls(residual_lifespan ~ log(proteome_size), 
                                       data = combined_anage_Birds_short2, 
                                       correlation=corBrownian(1, avgtettree_Birds_clean, form = ~species))
summary(birds_model_LSBM_proteome)

combined_anage_Birds_short2$pred_residual_lifespan <- as.numeric(predict(birds_model_LSBM_proteome, combined_anage_Birds_short2))


## birds per amino acid
### analysis of individual AAs (called "RAW" as opposed to the PCA transformed: but see manuscript for detailed methods of controlling for confounding variables)
birds_wide_aa <- combined_anage_Birds_short2 %>% 
  select(Class, species, aa_full, aa_freq) %>%
  pivot_wider(., names_from = aa_full,
              values_from = aa_freq)

## model raw and confit
models_birds_aa <- combined_anage_Birds_short2  %>%
  nest(data = -aa_full) %>%
  mutate(mdout = map_df(data, function(.x){
    md <- nlme::gls(residual_lifespan ~ scale(aa_freq), data=.x, correlation=corBrownian(1, avgtettree_Birds_clean, form = ~species))
    out <- data.frame(lwr = as.vector(confint(md, level = 0.9975)[2,1]),
                      slope = as.vector(md$coefficients[2]),
                      upr = as.vector(confint(md, level = 0.9975)[2,2]))
    return(out)
  })) %>%
  unnest(mdout) %>%
  select(-data)

models_birds_aa %>% arrange(aa_full)

## plot birds
ggplot(combined_anage_Birds_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                  col = aa_full, fill = aa_full)) + 
  geom_point() + 
  geom_smooth(combined_anage_Birds_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                         col = aa_full, fill = aa_full),
              method = 'lm') + 
  facet_wrap(~aa_full, ncol = 4, scales = "free") + 
  scale_fill_viridis_d('Amino acid') + 
  scale_color_viridis_d('Amino acid') + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "grey90")) + 
  labs(x = "Frequency", y = "Residual lifespan")







### * Amphibia 

## getting residual lifespan
combined_anage_Amphibia_short <- combined_anage_Amphibia %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs., weight = Adult.weight..g.,
         species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit()







### * Fish 

## getting residual lifespan
combined_anage_Fish_short <- combined_anage_Fish %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs.,  weight = Adult.weight..g.,
         species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit()

table(tettree_Fish$tip.label %in% unique(combined_anage_Fish_short$species))


avgtettree_Fish <- tettree_Fish
#avgtettree_Mammals$tip.label <- str_replace(avgtettree_Mammals$tip.label, "\\*", "")
##cleaning tree
## removing asterisks
avgtettree_Fish$tip.label <- str_replace(avgtettree_Fish$tip.label, "\\*", "")
avgtettree_Fish_clean <- drop.tip(avgtettree_Fish, avgtettree_Fish$tip.label[!avgtettree_Fish$tip.label %in% combined_anage_Fish_short$species])



avgtettree_Fish_clean$tip.label %in% unique(combined_anage_Fish_short$species)

## dropping tips
combined_anage_Fish_short2 <- combined_anage_Fish_short[combined_anage_Fish_short$species %in% avgtettree_Fish_clean$tip.label,]


length(avgtettree_Fish_clean$tip.label)


## controlling for body weight
Fish_model_LSBM <- nlme::gls(log(lifespan) ~ log(weight^0.5), 
                             data = combined_anage_Fish_short2, 
                             correlation=corBrownian(1, avgtettree_Fish_clean, form = ~species))

summary(Fish_model_LSBM)

## getting resicual lifespan
combined_anage_Fish_short2$residual_lifespan <- as.numeric(residuals(Fish_model_LSBM))


## checking 
Fish_model_LSBM_check <- nlme::gls(residual_lifespan ~ log(weight^0.5), 
                                   data = combined_anage_Fish_short2, 
                                   correlation=corBrownian(1, avgtettree_Fish_clean, form = ~species))
summary(Fish_model_LSBM_check)

## proteome size
Fish_model_LSBM_proteome <- nlme::gls(residual_lifespan ~ log(proteome_size), 
                                      data = combined_anage_Fish_short2, 
                                      correlation=corBrownian(1, avgtettree_Fish_clean, form = ~species))
summary(Fish_model_LSBM_proteome)

AIC(Fish_model_LSBM_proteome)


combined_anage_Fish_short2$pred_residual_lifespan <- as.numeric(predict(Fish_model_LSBM_proteome, combined_anage_Fish_short2))



## birds per amino acid
### analysis of individual AAs (called "RAW" as opposed to the PCA transformed: but see manuscript for detailed methods of controlling for confounding variables)
Fish_wide_aa <- combined_anage_Fish_short2 %>% 
  select(Class, species, aa_full, aa_freq) %>%
  pivot_wider(., names_from = aa_full,
              values_from = aa_freq)

## model raw and confit
models_Fish_aa <- combined_anage_Fish_short2  %>%
  nest(data = -aa_full) %>%
  mutate(mdout = map_df(data, function(.x){
    md <- nlme::gls(residual_lifespan ~ scale(aa_freq), data=.x, correlation=corBrownian(1, avgtettree_Fish_clean, form = ~species))
    out <- data.frame(lwr = as.vector(confint(md, level = 0.9975)[2,1]),
                      slope = as.vector(md$coefficients[2]),
                      upr = as.vector(confint(md, level = 0.9975)[2,2]))
    return(out)
  })) %>%
  unnest(mdout) %>%
  select(-data)

models_Fish_aa %>% arrange(aa_full)

## plot birds
ggplot(combined_anage_Fish_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                 col = aa_full, fill = aa_full)) + 
  geom_point() + 
  geom_smooth(combined_anage_Fish_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                        col = aa_full, fill = aa_full),
              method = 'lm') + 
  facet_wrap(~aa_full, ncol = 4, scales = "free") + 
  scale_fill_viridis_d('Amino acid') + 
  scale_color_viridis_d('Amino acid') + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "grey90")) + 
  labs(x = "Frequency", y = "Residual lifespan")

















### * Reptiles 

## getting residual lifespan
combined_anage_Reptiles_short <- combined_anage_Reptiles %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs., weight = Adult.weight..g.,
         species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit()


avgtettree_Reptile$tip.label %in% unique(combined_anage_Reptiles_short$species)




##cleaning tree
## removing asterisks
avgtettree_Reptile$tip.label <- str_replace(avgtettree_Reptile$tip.label, "\\*", "")
avgtettree_Reptile_clean <- drop.tip(avgtettree_Reptile, avgtettree_Reptile$tip.label[!avgtettree_Reptile$tip.label %in% combined_anage_Reptiles_short$species])



avgtettree_Reptile_clean$tip.label %in% unique(combined_anage_Reptiles_short$species)

## dropping tips
combined_anage_Reptile_short2 <- combined_anage_Reptiles_short[combined_anage_Reptiles_short$species %in% avgtettree_Reptile_clean$tip.label,]


length(avgtettree_Reptile_clean$tip.label)


## controlling for body weight
Reptile_model_LSBM <- nlme::gls(log(lifespan) ~ log(weight^0.5), 
                                data = combined_anage_Reptile_short2, 
                                correlation=corBrownian(1, avgtettree_Reptile_clean, form = ~species))

summary(Reptile_model_LSBM)

## getting resicual lifespan
combined_anage_Reptile_short2$residual_lifespan <- as.numeric(residuals(Reptile_model_LSBM))


## checking 
Reptile_model_LSBM_check <- nlme::gls(residual_lifespan ~ log(weight^0.5), 
                                      data = combined_anage_Reptile_short2, 
                                      correlation=corBrownian(1, avgtettree_Reptile_clean, form = ~species))
summary(Reptile_model_LSBM_check)

## proteome size
Reptile_model_LSBM_proteome <- nlme::gls(residual_lifespan ~ log(proteome_size), 
                                         data = combined_anage_Reptile_short2, 
                                         correlation=corBrownian(1, avgtettree_Reptile_clean, form = ~species))
summary(Reptile_model_LSBM_proteome)


combined_anage_Reptile_short2$pred_residual_lifespan <- as.numeric(predict(Reptile_model_LSBM_proteome, combined_anage_Reptile_short2))



## birds per amino acid
### analysis of individual AAs (called "RAW" as opposed to the PCA transformed: but see manuscript for detailed methods of controlling for confounding variables)
Reptile_wide_aa <- combined_anage_Reptile_short2 %>% 
  select(Class, species, aa_full, aa_freq) %>%
  pivot_wider(., names_from = aa_full,
              values_from = aa_freq)

## model raw and confit
models_Reptile_aa <- combined_anage_Reptile_short2  %>%
  nest(data = -aa_full) %>%
  mutate(mdout = map_df(data, function(.x){
    md <- nlme::gls(residual_lifespan ~ scale(aa_freq), data=.x, correlation = corBrownian(1, avgtettree_Reptile_clean, form = ~species))
    out <- data.frame(lwr = as.vector(confint(md, level = 0.9975)[2,1]),
                      slope = as.vector(md$coefficients[2]),
                      upr = as.vector(confint(md, level = 0.9975)[2,2]))
    return(out)
  })) %>%
  unnest(mdout) %>%
  select(-data)

models_Reptile_aa %>% arrange(aa_full)

## plot birds
ggplot(combined_anage_Fish_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                 col = aa_full, fill = aa_full)) + 
  geom_point() + 
  geom_smooth(combined_anage_Fish_short2, mapping = aes(x = scale(aa_freq), y = residual_lifespan,
                                                        col = aa_full, fill = aa_full),
              method = 'lm') + 
  facet_wrap(~aa_full, ncol = 4, scales = "free") + 
  scale_fill_viridis_d('Amino acid') + 
  scale_color_viridis_d('Amino acid') + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "grey90")) + 
  labs(x = "Frequency", y = "Residual lifespan")

















### * Sharks 

## getting residual lifespan
combined_anage_Sharks_short <- combined_anage_SharksRay %>% 
  mutate(aa_freq = aa_total_count/total_aa_overall_sp,
         aa_freq_codon = (aa_total_count/total_aa_overall_sp)/numb_codons,
         species = str_replace(species, " ", "_")) %>%
  select(Class, Genus, Species, lifespan = Maximum.longevity..yrs., weight = Adult.weight..g.,
         species, proteome_size = total_aa_overall_sp,
         aa_freq, aa_freq_codon, aa_full) %>%
  na.omit()


avgtettree_Sharks$tip.label %in% unique(combined_anage_Sharks_short$species)







## whole dataset/


combined_subset <- bind_rows(combined_anage_Mammals_short2, combined_anage_Birds_short2,
                             combined_anage_Fish_short2,  combined_anage_Reptile_short2) %>%
  ungroup() %>%
  select(Class, species, aa_full, aa_freq) %>%
  arrange(aa_full) %>%
  pivot_wider(., names_from = aa_full, values_from = aa_freq)

unique(combined_subset$Class)
combined_subset$Class2 <- ifelse(combined_subset$Class == "Chondrostei" | 
                                   combined_subset$Class == "Holostei" | 
                                   combined_subset$Class == "Teleostei", 
                                 "Actinopteri",
                                 paste(combined_subset$Class))

combined_subset$Class2 <- ifelse(combined_subset$Class2 == "Reptilia",
                                 "Reptillia",
                                 paste(combined_subset$Class2))


combined_subset_barplotdt <- combined_subset %>% pivot_longer(Alanine:Valine, names_to = "Amino_acid", values_to = "Frequency")


ggplot(combined_subset_barplotdt, aes(x = Class2, y = Frequency, fill = Class2)) + 
  facet_wrap(~Amino_acid, ncol = 4) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  scale_fill_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                    values = c('steelblue2', 'goldenrod1', 'firebrick1',
                               'chartreuse2', "grey20", "purple")) +
  labs(x = "Class", y = "Frequency")




summary(pccomp_anage_MammalsBirds)
## PCA on amino acids
pccomp_anage_MammalsBirds <- prcomp(combined_subset[3:22], scale. = FALSE)
#PCA_test_anage <- PCAtest(prcomp_data_anage[9:28]) ## PCA test

#PCA_test_anage <- PCAtest(prcomp_data_anage[9:28])
#summary(PCA_test_anage)

## data for plotting
pccomp_anage_MammalsBirds_df <- data.frame(combined_subset[c(1,2,23)],
                                           as.data.frame(pccomp_anage_MammalsBirds$x))


as.data.frame(pccomp_anage_MammalsBirds$rotation) %>% arrange(PC2)
## plot
PC13_plot_subset <- ggplot(pccomp_anage_MammalsBirds_df, aes(PC1, PC3, color=Class2)) +
  geom_point(size=2, alpha = 0.5) +
  stat_ellipse() + 
  theme_light() + 
  theme(panel.grid = element_blank())+ 
  scale_color_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                     values = c('steelblue2', 'goldenrod1', 'firebrick1',
                                'chartreuse2', "grey20", "purple"))


PC12_plot_subset <-ggplot(pccomp_anage_MammalsBirds_df, aes(PC1, PC2, color=Class2)) +
  geom_point(size=2, alpha = 0.5) +
  stat_ellipse() + 
  theme_light() + 
  theme(panel.grid = element_blank())+ 
  scale_color_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                     values = c('steelblue2', 'goldenrod1', 'firebrick1',
                                'chartreuse2', "grey20", "purple"))

PC12_plot + PC13_plot + PC12_plot_subset + PC13_plot_subset + plot_layout(guides="collect")



## Permanova
adonis2(
  pccomp_anage_MammalsBirds_df[ , c("PC1","PC2", "PC3")] ~ Class2,
  data = pccomp_anage_MammalsBirds_df,
  method = "euc"
)


### slopws as barplots
slopsdf <- bind_rows(cbind(models_Reptile_aa, Class =  "Reptillia"),
                     cbind(models_Fish_aa, Class =  "Actinopteri"),
                     cbind(models_birds_aa, Class =  "Aves"),
                     cbind(models_mammals_aa, Class =  "Mammalia")) 

ggplot(slopsdf, aes(x = Class, y = slope, fill = Class)) + 
  facet_wrap(~aa_full, ncol = 4) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  scale_fill_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                    values = c('steelblue2', 'goldenrod1', 'firebrick1',
                               'chartreuse2', "grey20", "purple")) +
  labs(x = "Class", y = "Slope")








## proteome size and residual lifespan
combined_ProteomeSize <- bind_rows(combined_anage_Fish_short2, combined_anage_Reptile_short2, combined_anage_Birds_short2,combined_anage_Mammals_short2) %>%
  select(-aa_full) %>%
  unique()

combined_ProteomeSize$Class2 <- ifelse(combined_ProteomeSize$Class == "Chondrostei" | 
                                         combined_ProteomeSize$Class == "Holostei" | 
                                         combined_ProteomeSize$Class == "Teleostei", 
                                       "Actinopteri",
                                       paste(combined_ProteomeSize$Class))

combined_ProteomeSize$Class2 <- ifelse(combined_ProteomeSize$Class2 == "Reptilia",
                                       "Reptillia",
                                       paste(combined_ProteomeSize$Class2))

ggplot() + 
  #  geom_hline(yintercept = 0, lty = 1, col = "black", linewidth = 0.2) +
  geom_point(data = combined_ProteomeSize, mapping = aes(x = log(proteome_size), y = residual_lifespan, col = Class2), alpha = 0.008, size = 2) + 
  geom_line(data = combined_ProteomeSize, mapping = aes(x = log(proteome_size), y = pred_residual_lifespan, col = Class2), linewidth = 1) +
  theme_light() + 
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 0),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  scale_colour_manual('Class', breaks = c("Actinopteri", "Amphibia", "Mammalia", "Reptillia", "Chondrichthyes", "Aves"),
                      values = c('steelblue2', 'goldenrod1', 'firebrick1',
                                 'chartreuse3', "grey20", "purple")) +
  labs(x = "log(Proteome size)", y = "Residual lifespan")


##
DataS1 <- combined_ProteomeSize %>%
  mutate(class = Class2) %>%
  select(-c(Class2, aa_freq, aa_freq_codon)) %>%
  unique() %>%
  select(class, species, residual_lifespan, proteome_size, lifespan, weight_g = weight)
DataS1
