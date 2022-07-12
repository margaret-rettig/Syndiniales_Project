
#combined = left_join(taxa_df, otus_df, by = "ASV" ) 
#combine data frames 
#head(combined)
taxa_df %>%
  filter(!is.na(Species))

print(otus_df)

sums <- rowSums(otus_df)
print(sums)

rowSums(otus_df) %>% 
  summary
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       6      32    1321     218   87077 

otu_intm <- rownames_to_column(otus_df, "ASV")
otu_data <- otu_intm %>% 
  pivot_longer(cols = -ASV, names_to = "otu", values_to = "raw_counts")

statistics <- otu_data %>% 
  group_by(ASV) %>% 
  summarize(mean(raw_counts), min(raw_counts), max(raw_counts))

#Relative abundances
sum <- otu_data %>% 
  group_by(otu) %>% 
  summarize(sum=sum(raw_counts))

otu_table <- otu_data %>% 
  group_by(otu) %>%
  mutate(relative_abundance = raw_counts / sum(raw_counts)) %>% 
  select(ASV, otu, relative_abundance)

taxa_intm <- taxa_df

colnames(otu_table)
colnames(taxa_intm)

merged <- left_join(otu_table, taxa_intm)


merged %>% 
  group_by(Order) %>% 
  summarize(median = median(relative_abundance)) %>% 
  arrange(desc(median))


sample_info <- rownames_to_column(samples, "otu")
sample_and_taxa <- left_join(sample_info, merged)

mean_abund <- sample_and_taxa %>% 
  group_by(otu, Order) %>% 
  summarize(mean_rel_abund = mean(relative_abundance)) %>% 
  ungroup() %>% 
  arrange(desc(mean_rel_abund))

top_order <- mean_abund %>% 
  group_by(Order) %>% 
  arrange(desc(mean_rel_abund)) %>% 
  top_n(25, mean_rel_abund) %>%
  pull(Order)

mean_abund %>%  #this isn't working quite right, revisit later
  ggplot(aes(x=Order, y=mean_rel_abund, color = otu)) + geom_boxplot()

#Phyloseq Package section
ntaxa(ASV_phyloseq)
nsamples(ASV_phyloseq)
sample_names(ASV_phyloseq)
sample_variables(ASV_phyloseq)

dino_group_I <- subset_taxa(ASV_phyloseq, Order == "Dino-Group-I")
plot_bar(dino_group_I, fill = "Genus", facet_grid = ~date)

dino_group_II <- subset_taxa(ASV_phyloseq, Order = "Dino-Group-II")
plot_bar(dino_group_II, fill = "Genus", facet_grid = ~date) 

ps_abundance <- transform_sample_counts(ASV_phyloseq, function(x) x/sum(x))

dino_group_III <- subset_taxa(ps_abundance, Order == "Dino-Group-III")
plot_bar(dino_group_III, "Family", fill = "Genus", facet_grid = ~date) 
#this plot is fairly decent look into this more

#ORDINATIONS

#plot_ordination(ps_abundance, type = "taxa", color = "Genus", title = "taxa")
#does not work(confused)

plot.richness


#Erin's file
GP.chl <- subset_taxa(ASV_phyloseq, Class =="Syndiniales")
plot_tree(GP.chl, color="SampleType", shape="Order", label.tips="Genus", size="Abundance")

plot_tree(ASV_phyloseq)

#Tree trial 
ntaxa(ASV_phyloseq)
physeq = prune_taxa(taxa_names(ASV_phyloseq)[1:50], ASV_phyloseq)
plot_tree(physeq)

#Ordination (Erin's File)
pcoa <- ordinate(ASV_phyloseq, method="MDS", distance="euclidean")
eigen_vals <- pcoa$values$Eigenvalues 

plot_ordination(ASV_phyloseq, pcoa, color="layer_epoch") +
  scale_color_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  scale_fill_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  geom_point(aes(colour=factor(layer_epoch))) +
  #geom_point(aes(shape=factor(type)),size=3) + 
  labs(col="layer_epoch") + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = depth_layer))+
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  theme_bw()+
  theme(legend.title = element_blank())
#okay so like what does this graph mean

syn_nmds <- ordinate(ASV_phyloseq, method="NMDS", distance="jaccard")
syn_nmds$stress #0.06 good fit
plot_ordination(ASV_phyloseq, syn_nmds, color="layer_epoch") +
  scale_color_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  scale_fill_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  geom_point(aes(colour=factor(layer_epoch))) +
  geom_point(size=3)+
  labs(col="layer_epoch") + 
  theme_bw()+
  theme(legend.title = element_blank())

#Alpha Diversity Measure Plot 
plot_richness(ASV_phyloseq, x="julianD", color="depth",  measures=c("Chao1", "Shannon")) + 
  geom_point(size=2) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  theme(legend.title = element_blank())

#Heatmap (come back to this later)
theme_set(theme_bw())

gpt <- subset_taxa(ASV_phyloseq, Class == "Syndiniales")
gpt <- prune_taxa(names(sort(taxa_sums(gpt), TRUE)[1:841]), gpt)
plot_heatmap(gpt)

gpac <- subset_taxa(ASV_phyloseq, Order == "Dino-Group-I")
plot_heatmap(gpac)

(p <- plot_heatmap(gpt, "NMDS", "bray", "depth", "Order"))


sample_and_taxa %>% 
  group_by(depth_layer)
  