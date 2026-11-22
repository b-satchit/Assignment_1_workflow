###### **********************
## GROUP 6 - TRICHINELLA SP. T8
##
## Group 6 - Ruhi Sajabi, Satchit Bhattacharya, Saniya Chadha, Phuc Pham, and Grace Enike.
##
## 2026-11-02
##
##***************************

## _ Packages used -------
rm(list = ls())
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(stringr)

#trich_main is the trichinella sp. we will be comparing our assigned trich. group to
#trich_t8 is the trichinella sp. our group was assigned

trich_main = read_csv("../data/GCA_001447595.1.csv")
trich_main |> view()

trich_t8 = read_csv("../data/GCA_001447745.1.csv")
trich_t8


## Identifying number of unique classes in the trichinella t8 sp.

trich_t8$matching_class |> unique() |> length()

## Identifying unique classes in the trichinella spiralis

trich_main$matching_class |> unique() |> length()

## Unique classes by matching class

trich_t8 |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n)) 

trich_main |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n))

## Identify the retrotransposons in the trich_t8 class. From Copilot, everything starting with LTR

trich_t8_retro = trich_t8 |> 
  filter(str_detect(matching_class, "^(LTR|LINE|SINE|PLE)"))

str_detect(trich_t8$matching_class, "^(LTR|LINE|SINE|PLE)")

## Identify the retrotransposons in the trich_main class. From Copilot, everything starting with LTR

trich_main_retro = trich_main |> 
  filter(str_detect(matching_class, "^(LTR|PLE|SINE|LINE)"))

str_detect(trich_main$matching_class, "^(LTR|PLE|SINE|LINE)")

unique(trich_t8_retro$matching_class)


## Percentage deletions for trich_t8_retro

trich_t8_retro |> 
  ggplot(aes(x = matching_class, y = perc_del)) +
  geom_boxplot() +
  labs(title = "Trich T8 Retrotransposons", y ="Percent Deletions (%)" , x = "Retrotransposons class")

## Percentage deletions for trich_main_retro

trich_main_retro |> 
  ggplot(aes(x = matching_class, y = perc_del)) +
  geom_boxplot() +
  labs(title = "T. spiralis Retrotransposons", y ="Percent Deletions (%)" , x = "Retrotransposons class")
 
## Percentage insertions for trich_t8_retro

trich_t8_retro |> 
  ggplot(aes(x = matching_class, y = perc_insert)) +
  geom_boxplot() +
  labs(title = "Trich T8 Retrotransposons", y ="Percent Insertions (%)" , x = "Retrotransposons class")

## Percentage insertions for trich_main_retro

trich_main_retro |> 
  ggplot(aes(x = matching_class, y = perc_insert)) +
  geom_boxplot() +
  labs(title = "T. spiarlis Retrotransposons", y ="Percent Insertions (%)" , x = "Retrotransposons class")

## Percentage deletion and insertion in trich_te_retro

trich_t8_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_point(aes(color = matching_class)) +
  geom_smooth()

## 2nd Attempt at easier to interpret visual via gpplot

trich_t8_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_boxplot(aes(color = matching_class)) +
  geom_smooth()
 
## Percentage deletion and insertion in trich_main

trich_main_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_point(aes(color = matching_class)) +
  geom_smooth()

## Alternative visuals 
trich_main_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_boxplot(aes(color = matching_class)) +
  geom_smooth()


## Comparing retrotransposons insertion and deletion patterns between both species

retro_combined = bind_rows(trich_main_retro |> mutate(species = "T. spiralis"),
                           trich_t8_retro |> mutate(species = "T8"))

retro_combined |>
  ggplot(aes(x = perc_del, y = perc_insert, color = species)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Insertion and Deletion Patterns in Retrotransposons", x = "Percentage deletions (%)", y = "Percentage insertions (%)", color = "Organism")

summary(trich_t8_retro$perc_insert)
summary(trich_t8_retro$perc_del)


############################################### Show individual classes
retro_combined |>
  ggplot(aes(x = matching_class,
             y = perc_insert,
             fill = species)) + #fill species colours the box plot
  labs(
    title = "Insertion Percentages Across Retrotransposon Classes",
    x = "Retrotransposon Class",
    y = "Percentage Insertions (%)"
  ) +
  geom_boxplot() +
  facet_wrap(~species) + #make separate panels for each species instead of one graph 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# angle 45 = rotate the x axis to fit 
#h_just = aligns the labels 

###
retro_combined |>
  ggplot(aes(x = matching_class,
             y = perc_del,
             fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  labs(
    title = "Deletion Percentages Across Retrotransposon Classes",
    x = "Retrotransposon Class",
    y = "Percentage Deletions (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Grace's code template stops here ----
## Novelty - Retrotransposon Composition - Stacked Bar Chart ----
retro_combined |>
  count(species, matching_class) |>
  group_by(species) |> #Normalize counts into percentages for each genome
  mutate(percent = n / sum(n) * 100) |>
  ggplot(aes(x = species,
             y = percent,
             fill = matching_class)) +
  geom_col() +
  labs(
    title = "Retrotransposon Composition",
    x = "Genome",
    y = "Percentage of Annotations",
    fill = "Retrotransposon Class"
  ) +
  theme_minimal()

## Novelty - Retrotransposon Composition - Percentage Differences ----
retro_combined |>
  count(species, matching_class) |>
  group_by(species) |>
  mutate(percent = n / sum(n) * 100) |>
  select(species, matching_class, percent) |> #Place genome percentages side by side
  pivot_wider(
    names_from = species,
    values_from = percent,
    values_fill = 0
  ) |>
  mutate(difference = T8 - `T. spiralis`) |> #Calculate percentage-point differences (T8 minus T. spiralis)
  arrange(desc(abs(difference)))

## Novelty - Statistical Analysis - Chi-square Test ----
class_counts <- retro_combined |>
  count(matching_class, species) |>
  pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  )# Count annotations in each class for both genomes

count_matrix <- class_counts |>
  select(`T. spiralis`, T8) |>
  as.matrix()
# Creates a table for the statistical test

chi_result <- chisq.test(count_matrix)
chi_result
# Run the chi-square test and view results

##diana's proposed edits 

##edit 1 - when you run the original chi-square test, it gives a warning message. The group also noted the low frequency of some retrotransposon classes means the chi-squared test is not fully reliable 
#--Wilcoxon tests for insertion/deletion--
#The original compared species with boxplots only, so "similar" had no test.
#Wilcoxon rank-sum tests are used because the percentages are skewed (many zeros, some outliers).
wilcox.test(perc_insert ~ species, data = retro_combined)
wilcox.test(perc_del ~ species, data = retro_combined)


##edit 2 - improved visualization 
#in the original boxplot, a single extreme outlier (~91% in T8 LINE/CR1) stretched the original y-axis and flattened all boxes near 0%. 
#coord_cartesian() zooms into 0-15% visually without deleting the outlier from the data.


retro_combined |>
  ggplot(aes(x = matching_class, y = perc_insert, fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    title    = "Insertion Percentages Across Retrotransposon Classes",
    subtitle = "Y-axis zoomed to show box distributions (single T8 LINE/CR1 outlier out of view)",
    x        = "Retrotransposon Class",
    y        = "Percentage Insertions (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##edit 3 -- Sample Size & Data Sparsity Check--
# Added this table to count annotations per class across both genomes to flag sparse classes
# (n < 5 in at least one species; 8/13 total) that violate chi-square test assumptions and can bias percentage comparisons.


low_annotation_counts <- retro_combined |>
  count(matching_class, species) |>
  complete(matching_class, species, fill = list(n = 0)) |>
  group_by(matching_class) |>
  mutate(
    total_n       = sum(n),
    sample_status = ifelse(min(n) < 5, "Sparse (n < 5)", "Adequate")
  ) |>
  pivot_wider(names_from = species, values_from = n) |>
  arrange(total_n)

low_annotation_counts

