# Assuming id_df contains your individual-level data with Island, Species, sex, feeding

# 1. Species distribution by Island
species_island_dist <- full_id_df %>%
  group_by(Island, Species) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Count / sum(Count), TotalIsland = sum(Count))

print("Species Distribution by Island:")
print(species_island_dist)

ggplot(species_island_dist, aes(x = Island, y = Proportion, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Species Composition by Island", y = "Proportion of Individuals", fill = "Species") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Interpretation:
# - If one species is only on one island, you won't be able to separate island effects from species effects for that species.
# - If species are unevenly distributed, island-level random effects become very important.

# 2. Feeding distribution by Island and Species
feeding_island_species_dist <- full_id_df %>%
  group_by(Island, Species, feeding) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Count / sum(Count), TotalGroup = sum(Count))

print("Feeding Distribution by Island and Species:")
print(feeding_island_species_dist)

ggplot(feeding_island_species_dist, aes(x = Island, y = Proportion, fill = feeding)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Species) +
  labs(title = "Feeding Regimes by Island and Species", y = "Proportion of Individuals", fill = "Feeding") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Interpretation:
# - If a specific species-feeding combination is heavily concentrated on one island, it suggests potential confounding between that interaction and island-specific effects.

# 3. Sex distribution by Island and Species
sex_island_species_dist <- id_df %>%
  group_by(Island, Species, sex) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Count / sum(Count), TotalGroup = sum(Count))

print("Sex Distribution by Island and Species:")
print(sex_island_species_dist)

ggplot(sex_island_species_dist, aes(x = Island, y = Proportion, fill = sex)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Species) +
  labs(title = "Sex Distribution by Island and Species", y = "Proportion of Individuals", fill = "Sex") +
  theme_minimal() +
  scale_fill_manual(values = c("female" = "pink", "male" = "lightblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Interpretation:
# - Similar to feeding, check if sex ratios or distributions are highly skewed for certain species/islands.