library(tidyverse)
library(umap)
library(jsonlite)
library(mirt)
library(dbscan)
library(ggrepel)

## Data input
# Path
path_fs <- 'C:/.../_results/'
# Props
prop_files <- list.files( path = path_fs, pattern = 'props*' )
# Rslts
rslt_files <- list.files( path = path_fs, pattern = 'rslt*' )

# Process the files
data_rslt <- tibble(compound_id = c(NA), target_id = c(NA), pa_pi = c(NA), target_name = c(NA), type = c(NA) )
for (i in seq(1:length(rslt_files))) {
	rslt <- read_json( str_glue('{path_fs}{rslt_files[i]}') )
	rslt_d <- rslt$directi |> enframe() |>
					 unnest_wider(col = value) |>
					 mutate(type = 'direct') |>
					 select(-name) |>
					 rename(target_id = tid, pa_pi = prob, target_name = pname) |>
					 mutate(compound_id = str_glue('cmpnd_{i}'))
	rslt_i <- rslt$indericti |> enframe() |>
					 unnest_wider(col = value) |>
					 mutate(type = 'indirect') |>
					 select(-name) |>
					 rename(target_id = tid, pa_pi = prob, target_name = pname) |>
					 mutate(compound_id = str_glue('cmpnd_{i}'))
	data_rslt <- bind_rows(data_rslt, rslt_d, rslt_i)
}

data_prop <- tibble(compound_id = c(NA), target_id = c(NA), go = c(NA), class = c(NA) )
for (i in seq(1:length(prop_files))) {
	prop <- read_json( str_glue('{path_fs}{prop_files[i]}') ) |>
				enframe() |>
				unnest_wider(col = value) |>
				select(-name, -type) |>
				rename( target_id = tid, class = classif ) |>
				mutate(compound_id = str_glue('cmpnd_{i}'))
	data_prop <- bind_rows(data_prop, prop)
}
data_prop <- data_prop |> filter(!is.na(target_id)) |>
			    separate_longer_delim(go, delim = ',')
data_rslt <- data_rslt |> filter(!is.na(target_id)) |>
			    		  mutate(pa_pi = as.double(pa_pi))

## Get the predicted interactions (cmpnd-target)
predicted_ints <- data_rslt |> select(compound_id, target_id) |> distinct()

## Prepare the files for output
# Predictions
prediction_output <- data_rslt |>
					 arrange(desc(pa_pi))
# GO terms
go_output <- data_prop |> semi_join( predicted_ints, by = join_by(compound_id, target_id) ) |>
						select(compound_id, go) |>
						group_by(compound_id, go) |>
						summarize(term_count = n(), .groups = 'keep') |>
						arrange(desc(term_count)) |>
						arrange(compound_id)
# Classification
class_output <- data_prop |> semi_join( predicted_ints, by = join_by(compound_id, target_id) ) |>
						select(compound_id, target_id, class) |>
						distinct() |>
						select(compound_id, class) |>
						group_by(compound_id, class) |>
						summarize(term_count = n(), .groups = 'keep') |>
						arrange(desc(term_count)) |>
						arrange(compound_id)

## Output the results
write_tsv(prediction_output, str_glue('{path_fs}8-steroids_PT_predictions_gathered.tab'))
write_tsv(go_output, str_glue('{path_fs}8-steroids_PT_GOterms_gathered.tab'))
write_tsv(class_output, str_glue('{path_fs}8-steroids_PT_Classes_gathered.tab'))

## Basic stats
n_targets <- data_rslt |> pull(target_id) |> unique() |> length()
n_d_targets <- data_rslt |> filter(type == 'direct') |> pull(target_id) |> unique() |> length()
n_i_targets <- data_rslt |> filter(type == 'indirect') |> pull(target_id) |> unique() |> length()

## Extended stats
# Number of targets
n_targets <- data_rslt |> group_by(compound_id) |> summarize(n_targets = n())
# Pa-Pi distributions
pa_pi_distr_plot <- ggplot(data_rslt, aes(pa_pi)) +
					geom_histogram() +
					geom_vline(xintercept = .2) +
					theme_minimal() +
					facet_wrap(~compound_id)
pa_pi_distr_plot
# Save the plot
ggsave(str_glue('{path_fs}Pa-Pi_distr.png'), width = 7, height = 7, units = 'in', dpi = 300)


### Do the further analysis using only targets with Pa-Pi > .2
## Filter and process
data_rslt_f <- data_rslt |> filter(pa_pi > .2)
n_targets_f <- data_rslt_f |> pull(target_id) |> unique() |> length()

data_rslt_f <- data_rslt_f |> mutate( target_id = str_glue('{target_id}_{type}') ) |>
						select(compound_id, target_id, pa_pi) |>
						mutate(pa_pi = 1) |>
						pivot_wider( names_from = target_id, values_from = pa_pi, values_fill = 0 )
# Get the col summ
f_colsum <- data_rslt_f |> select(-compound_id) |>
			summarise_all(sum)
# SO, there are columns without the variance, they should be omitted

# Get the colnames with zero and near zero variance
cols <- f_colsum |> pivot_longer(cols = everything(), names_to = 'clmn', values_to = 'val') |>
					filter(val < 6) |>
					arrange(desc(val)) |>
					pull(clmn)

# Select columns from data_rslt_f
data_rslt_f <- data_rslt_f |> select(all_of(cols)) |> relocate(any_of(cols))

# Prepare the data_rslt_f to work with it latter
data_rslt_f_id <- data_rslt_f |> mutate(compound_id = str_glue('cmpnd_{row_number()}'))
# Try basic UMAP with the predicted activities
set.seed(11)
umap_basic_zerotwo <- umap(data_rslt_f, n_neighbors = 7, preserve.seed = TRUE)
data_rslt_f_id$x <- umap_basic_zerotwo$layout[,1]
data_rslt_f_id$y <- umap_basic_zerotwo$layout[,2]
data_rslt_f_id <- data_rslt_f_id |>
				  unite(allkey, starts_with('CHEMBL'), sep = '', remove = FALSE)
# Plot the data
plot_basic <- ggplot( data_rslt_f_id, aes(x, y, label = compound_id) ) +
				coord_fixed(ratio = 1) +
				geom_point() +
				geom_text_repel() +
				theme_minimal()
plot_basic

# Save the plot
ggsave(str_glue('{path_fs}basic_UMAP.png'), width = 7, height = 7, units = 'in', dpi = 300)

## Try UMAP with the neighboring points
# Prepare the available points
availableSpace <- data_rslt_f_id |> select(-allkey, -compound_id, -x, -y) |> as.matrix()
mode(availableSpace) <- 'integer'
# Generate the whole space
featureSpace <- thetaComb(theta = c(0,1), nfact = ncol(availableSpace), intercept = FALSE)
mode(featureSpace) <- 'integer'

# Create the shared path
path <- matrix( data = NA, nrow = nrow(availableSpace)^2 * ncol(availableSpace), ncol = ncol(availableSpace) )
# Populate the shared path
path_counter <- 0
for (i in seq(1:nrow(availableSpace))) {
	start_point <- availableSpace[i,]
	inter_point <- availableSpace[i,]
	for (k in seq(1:nrow(availableSpace))) {
		end_point <- availableSpace[k,]
		for (q in seq(1:length(end_point))) {
			if (start_point[q] != end_point[q]) {
				inter_point[q] <- end_point[q]
			}
			path_counter <- path_counter + 1
			path[path_counter,] <- inter_point
		}
	}
}
# Create the first row of the nearest neighbors for shared path
# 'distinct' could be placed at another point with totally different effect
sampleSpace_first <- kNN(featureSpace, ncol(path), path)
set.seed(43)
neighbors_first_rows <- sampleSpace_first$id |> as_tibble() |> 
				 pivot_longer(cols = everything(), names_to = 'clmn', values_to = 'val') |>
				 distinct() |>
				 sample_frac(size = .5)
neighbors_rows <- neighbors_first_rows |> pull(val) 
# Prepare the actual points
availableSpace <- availableSpace |> as_tibble() |>
					unite(allkey, sep = '', remove = FALSE) |>
					mutate(compound_id = str_glue('cmpnd_{row_number()}'))
# Prepare the neighboring points
neighborSpace <- featureSpace[neighbors_rows,] |> as_tibble() |>
			   unite(allkey, sep = '', remove = FALSE) |>
			   left_join( data_rslt_f_id |> select(allkey, compound_id) ) |>
			   mutate( compound_id = if_else(is.na(compound_id), 'background', compound_id) )
colnames(neighborSpace) <- colnames(availableSpace)
# Bind points
sampleSpace <- bind_rows(availableSpace, neighborSpace) |> distinct()

# UMAP the points along with their surroundings
umap_zerotwo <- umap(sampleSpace |> select(-compound_id, -allkey), n_neighbors = 25, preserve.seed = TRUE)
sampleSpace$x <- umap_zerotwo$layout[,1]
sampleSpace$y <- umap_zerotwo$layout[,2]

# Plot the data
plot_cmpnds <- ggplot( sampleSpace |> filter(compound_id != 'background'), aes(x, y, label = compound_id) ) +
				coord_fixed(ratio = 1) +
				geom_point() +
				geom_text_repel() +
				theme_minimal()
plot_cmpnds
# Save the plot
ggsave(str_glue('{path_fs}extended_UMAP.png'), width = 7, height = 7, units = 'in', dpi = 300)