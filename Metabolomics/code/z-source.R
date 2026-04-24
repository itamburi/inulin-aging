

library(tidyverse)
library(here)
library(purrr)
library(broom)
library(cowplot)
library(ComplexHeatmap)
library(grid)
library(ggplot2)
# ========== set_tissue_levels() ==========
# --
set_tissue_levels <- function(df) { df %>% dplyr::mutate( tissue = factor(tissue, levels=c("Intestine","Colon","Liver","Spleen","Head","Leg","Kidney","Skin","Heart","Lung")))}

# ========== select_top_isomer3() ==========
# -- improved algorithm that identifies RT clusters and selects best cluster based on n samples where that cluster was Max IC

select_top_isomer3 <- function(df, rt_tol_min = 20/60) {
  # df columns expected: name, mz, rt.min, ion.count,
  # optional: raw.file, sequence.id, formula, ion.mode
  
  # DEBUG ONLY (commented out so it doesn't override df arg)
  #df <- batch_list[[1]]; rt_tol_min <- 20/60
  
  cat(paste0("\nCalling select_top_isomer3() on ", deparse(substitute(df)), "\n"))
  
  # --- Normalize column names to lower case ---
  df <- df %>% dplyr::rename_with(~tolower(.))
  
  needed <- c("name", "rt.min", "ion.count")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  has_seqid <- "sequence.id" %in% names(df)
  
  # Ensure numeric
  df <- df %>%
    dplyr::mutate(
      rt.min    = as.numeric(rt.min),
      ion.count = as.numeric(ion.count)
    )
  
  # --- Step 1: Per annotation (name, rt.min), calculate median ion count across samples ---
  rt_level <- df %>%
    dplyr::filter(!is.na(rt.min)) %>%
    dplyr::group_by(name, rt.min) %>%
    dplyr::summarise(
      median_ic = median(ion.count, na.rm = TRUE),
      n_obs     = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(median_ic = dplyr::if_else(is.na(median_ic), 0, median_ic))
  
  # --- Step 2: Assign clusters based on RT. Clusters defined by RT gap ≤ rt_tol_min ---
  rt_clustered <- rt_level %>%
    dplyr::arrange(name, rt.min) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(
      gap         = (rt.min - dplyr::lag(rt.min, default = dplyr::first(rt.min))),
      new_cluster = dplyr::if_else(dplyr::row_number() == 1, TRUE, gap > rt_tol_min),
      cluster_id  = cumsum(new_cluster)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-gap, -new_cluster)
  
  # --- Step 3a: Per cluster, calculate median IC, # RTs and min RT ---
  # ... initially we selected annotations based on this criteria alone. However step 3b is more robust
  cluster_stats <- rt_clustered %>%
    dplyr::group_by(name, cluster_id) %>%
    dplyr::summarise(
      representative_rt = min(rt.min, na.rm = TRUE),
      cluster_median_ic = median(median_ic, na.rm = TRUE),
      n_rts             = dplyr::n(),
      .groups           = "drop"
    )
  
  # --- Step 3b: Per cluster, determine the number of samples (seq.id) where the IC in that cluster was the max (n_max) ---
  if (has_seqid) {
    # Per (name, rt.min): how many samples had this RT as max IC?
    n_max_rt <- df %>%
      dplyr::group_by(name, sequence.id) %>%
      dplyr::filter(ion.count == max(ion.count, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::count(name, rt.min, name = "n_max") %>%
      dplyr::arrange(dplyr::desc(n_max))
    
    # Per (name, cluster_id): how many samples had ANY RT in this cluster as max IC?
    n_max_cluster <- df %>%
      dplyr::group_by(name, sequence.id) %>%
      dplyr::filter(ion.count == max(ion.count, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(name, sequence.id, rt.min) %>%
      dplyr::inner_join(
        rt_clustered %>% dplyr::select(name, rt.min, cluster_id),
        by = c("name", "rt.min")
      ) %>%
      dplyr::distinct(name, cluster_id, sequence.id) %>%
      dplyr::count(name, cluster_id, name = "n_max_cluster")
  } else {
    # No sequence.id → no vote information; fall back to zeros
    n_max_rt <- dplyr::tibble(name = character(), rt.min = numeric(), n_max = integer())
    n_max_cluster <- dplyr::tibble(name = character(), cluster_id = integer(), n_max_cluster = integer())
  }
  
  # --- Step 3c: Rank clusters, prioritizing n_max_cluster ---
  cluster_stats <- cluster_stats %>%
    dplyr::left_join(n_max_cluster, by = c("name", "cluster_id")) %>%
    dplyr::mutate(n_max_cluster = tidyr::replace_na(n_max_cluster, 0L)) %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(
      # *** Rank order criteria:
      dplyr::desc(n_max_cluster),
      dplyr::desc(cluster_median_ic),
      dplyr::desc(n_rts),
      representative_rt,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      cluster_rank = dplyr::row_number(),
      name_isomer  = dplyr::if_else(cluster_rank == 1,
                                    name,
                                    paste0(name, " : isomer ", cluster_rank - 1))
    ) %>%
    dplyr::ungroup()
  
  # --- Step 4: Map every (name, rt.min) to its cluster and label ---
  rt_map <- rt_clustered %>%
    dplyr::inner_join(cluster_stats, by = c("name", "cluster_id")) %>%
    dplyr::select(
      name, rt.min, cluster_id, cluster_rank, name_isomer, representative_rt
    ) %>%
    dplyr::left_join(n_max_rt, by = c("name", "rt.min")) %>%
    dplyr::mutate(n_max = tidyr::replace_na(n_max, 0L))
  
  # --- Step 5: Attach cluster labels back to the raw dataframe ---
  df_tagged <- df %>%
    dplyr::inner_join(rt_map, by = c("name", "rt.min")) %>%
    dplyr::mutate(name = name_isomer) %>%
    dplyr::select(-name_isomer)
  
  # --- Summary outputs ---
  primary_clusters <- cluster_stats %>%
    dplyr::filter(cluster_rank == 1) %>%
    dplyr::distinct(name, cluster_id) %>%
    dplyr::mutate(primary = TRUE)
  
  # ** total RT ranges **
  name_rt_ranges <- rt_level %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      rt.range = max(rt.min, na.rm = TRUE) - min(rt.min, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::left_join(primary_clusters[, c("name", "cluster_id")], by = "name") %>%
    dplyr::arrange(dplyr::desc(rt.range))
  
  # ** cluster ranges **
  cluster_rt_ranges <- rt_clustered %>%
    dplyr::left_join(primary_clusters, by = c("name", "cluster_id")) %>%
    dplyr::group_by(name, cluster_id, primary) %>%
    dplyr::summarise(
      clust_rt_range = max(rt.min, na.rm = TRUE) - min(rt.min, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    dplyr::arrange(name, cluster_id, dplyr::desc(clust_rt_range))
  
  repeats <- name_rt_ranges %>% dplyr::filter(rt.range > rt_tol_min)
  cat("\n", paste0(length(unique(repeats$name))," metabolites contain multiple RTs greater than 20 seconds apart\n"))
  
  all_others <- name_rt_ranges %>%
    dplyr::filter(rt.range <= rt_tol_min)
  cat("\n", paste0(length(unique(all_others$name))," metabolites contain multiple RTs within 20 seconds (or are singletons)\n"))
  
  return(df_tagged)
}
# ========== Function: 'select_top_isomer' (updated Sept 2025) ==========
# --
# ****  ****
# *** Background:
# This function improves on rm_repeat_names() but instead of removing the replicate name annotations altogether,
# it labels them as name : isomer 1, name : isomer 2, etc.
#
# *** Pre requisites:
# Run this function by group splitting the melted data into a list of batches: batch_list = exp_data  %>% group_by( Raw.File ) %>% group_split(.)
# lapply this function to the list of batches
#
# *** Logic:
# First, establish the annotated metabolites which have variant RTs across the batch
# We will assert that if the RT range for a metabolite in the batch is > 0.3 min then we need to determine which RT is most precise
# "divergent.repeats"
# Otherwise, if the RT range of a metabolite is <= 1, it is just RT drift. We will group these together and sum these
# "all.others"
select_top_isomer = function(df){
  
  #df = batch_list[[2]] # Debug line
  cat(paste0("\nCalling select_top_isomer() on ", deparse(substitute(df)),"\n"))
  
  df = df %>% rename_with(~ tolower(.))
  
  # (1) Establish the metabolite names that are > 0.3 RT range across the batch
  repeats = df %>%
    group_by( name ) %>%
    summarise(
      rt.range = max(rt.min, na.rm = TRUE)-min(rt.min, na.rm = TRUE),
      # range.mzcloud = max(mzcloud.conf, na.rm = TRUE) - min(mzcloud.conf, na.rm = TRUE)
      ## ...moving forward i'm going to ignore the mzcloud confidence for these data bc it doesnt seem to be variant between diff RTs of the same NAME
    ) %>% 
    filter(
      rt.range > 0.3 # we need to understand how changing this affects the output
      # when > 1.0 we have 485 metabolites; when > 0.3 we have 641 metabolites;
      # there are 1070 unique metabolites in the dataset so about half have repeat/isomer annotations
    ) 
  # these are the only Names to which we need to apply selective criteria
  n=length(unique(repeats$name))
  cat("\n",paste0(n," metabolites contatin multiple RTs greater than 20 seconds apart\n"))
  
  # (2) Select only the ICs where the metabolite range > 0.3,
  #     and determine the [Name, MZ, RT.min] with the highest median IC of the batch
  #     Label the rest with " : isomer n"
  divergent.repeats = df %>%
    filter(
      name %in% repeats$name,
      is.na(ion.count) == FALSE
    ) %>%
    mutate(
      value = ion.count
    ) %>%
    group_by( name, mz, rt.min) %>% # Multiple RT/MZ values in the batch are associated with the same name. For each, calculate the median
    summarise(
      
      # There is probably an even more robust set of selection criteria we can make depending on how the IC data are distributed
      
      # ** consider the delta mass
      # ** if second highest median is far away, keep both but label the metab as #2
      
      # For now we are just going to pick the max median
      batch.median.ic = median(value, na.rm = TRUE),
      
      
      # min= min(value),
      # max = max(value),
      # IQR = format(IQR(value),  scientific = TRUE, digits = 2),
      # var.mean.coef = sd(value, na.rm = TRUE)/mean(value),
      # shapiro.test = shapiro.test(log10(value))$p.value,
      # count = n()
      
    ) %>%
    group_by(name) %>%
    mutate(
      rank = rank(-batch.median.ic),
      name.isomer = ifelse(rank != 1, paste0(name, " : isomer ", rank-1), name),
      rank = NULL
    )
  
  
  # (3) Make a reduced data frame of all other metabolites where the RT range was < 0.3
  all.others = df %>%
    filter(
      !(name %in% repeats$name)
    ) %>%
    mutate( name.isomer = name) %>%
    group_by( name, name.isomer, mz, rt.min) %>% # Name and Name.isomer are the same for this df. Just need the .isomer column to bind to divergent.repeats
    summarise()
  n=length(unique(all.others$name))
  cat("\n",paste0(n," metabolites contatin multiple RTs within 20 seconds\n"))
  
  # (3) Bind the two data frames together to yeild a df of only the metabolite annotations we care about
  # - At this point we have:
  #     > divergent.repeats = df with top isomer and lower abundance "isomers" labeled for repeat entries > 20 sec apart in RT
  #     > all.others = df where repeated names have RT range within 20 seconds
  df2 = rbind( divergent.repeats[,!names(divergent.repeats) %in% c("batch.median")], all.others )
  
  
  # (4) Inner joint the metabolite annotations we care about with the original data.
  #  ... before in rm_repeat_names this would have reduced the # of rows vs exp_data2
  # ... but now since were just labeling repeats as isomer, we expect equal # of rows
  df3 = inner_join(df, df2)
  
  
  # (5) Lastly, now that we have dealt with the RT ranges > 0.3, we address those RT ranges < 0.3:
  #     Repeat entries with RT range < 0.3 still remain from binding all.others.
  #     In these instances, I am making the assertion that this is simply RT drift
  #     Therefore, on a per-sample basis, for every metabolite NAME I will take the median MZ and RT values,
  #               and sum the Ion Count vales. 
  
  # ** still needs work; not the perfect approach because we do the grouping by samples
  # ** Technically, between samples we could have the same metabolite Name with different RTs between samples
  # **    ... (the spirit of this all is to have unified metabolite annotations)
  
  df4 = df3 %>%
    group_by(raw.file, sequence.id, name.isomer, formula) %>%
    summarise(
      mz = min(mz),
      rt.min = min(rt.min),
      ion.count = sum(ion.count),
      #Ion.Mode = Ion.Mode
    ) %>%
    ungroup() %>%
    mutate(
      name = name.isomer,
      name.isomer = NULL
    ) %>%
    select(raw.file, sequence.id, name, formula, mz, rt.min, ion.count)
  
  n = length(unique(df4$name))
  cat("\n",paste0(n," metabolites including isomers"))
  
  filt = subset(df4, grepl("isomer", df4$name) == FALSE )
  n = length(unique(filt$name))
  cat("\n",paste0(n," metabolites isomers removed"))
  
  # metabs = unique( df4[, c("Name", "Ion.Mode")] )
  
  # length(unique(df3$Name))
  # length(unique(df4$Name))
  # Still same # metabolites but we removed redundancies
  
  return(df4)
}

# ========== Function: 'rm_repeat_names()' ==========
# --
######## ****  ****
# *** Pre requisites:
# Run this function by group splitting the melted data into a list of batches: batch_list = exp_data  %>% group_by( Raw.File ) %>% group_split(.)
# lapply this function to the list of batches

# *** Logic:
# First, establish the annotated metabolites which have variant RTs across the batch
# We will assert that if the RT range for a metabolite in the batch is > 1 min then we need to determine which RT is most precise
# "divergent.repeats"
# Otherwise, if the RT range of a metabolite is <= 1, it is just RT drift. We will group these together and sum these
# "all.others"

rm_repeat_names = function(df){

  #df= batch_list[[1]] # Debug line
  
  
  # (1) Establish the metabolite names that are > 1 RT range across the batch
  repeats = df %>%
    group_by( Name ) %>%
    summarise(
      rt.range = max(RT.min, na.rm = TRUE)-min(RT.min, na.rm = TRUE),
      # range.mzcloud = max(mzcloud.conf, na.rm = TRUE) - min(mzcloud.conf, na.rm = TRUE)
      ## ...moving forward i'm going to ignore the mzcloud confidence for these data bc it doesnt seem to be variant between diff RTs of the same NAME
    ) %>% 
    filter(
      rt.range > 1.0 # we need to understand how changing this affects the output
      # 
    ) 
  # these are the only Names to which we need to apply selective criteria
  
  
  
  # (2) Select only the ICs where the metabolite range > 1,
  #     and extract only the [Name, MZ, RT.min] with the highest median IC of the batch
  
  divergent.repeats = df %>%
    filter(
      Name %in% repeats$Name,
      is.na(Ion.Count) == FALSE
    ) %>%
    mutate(
      value = Ion.Count
    ) %>%
    group_by( Name, MZ, RT.min) %>% # Multiple RT/MZ values in the batch are associated with the same name. For each, calculate the median
    summarise(
      
      # There is probably an even more robust set of selection criteria we can make depending on how the IC data are distributed
      
      # ** consider the delta mass
      # ** if second highest median is far away, keep both but label the metab as #2
      
      # For now we are just going to pick the max median
      batch.median = median(value),
      
      
      # min= min(value),
      # max = max(value),
      # IQR = format(IQR(value),  scientific = TRUE, digits = 2),
      # var.mean.coef = sd(value, na.rm = TRUE)/mean(value),
      # shapiro.test = shapiro.test(log10(value))$p.value,
      # count = n()
      
    ) %>%
    group_by(Name) %>%
    # mutate(
    #   rank = rank(-batch.median),
    #   Name = ifelse(rank != 1, paste0(Name, " : isomer ", rank-1), Name),
    #   rank = NULL
    #   )
    top_n(1, batch.median)
  
  
  # (3) Make a reduced data frame of all other metabolites where the RT range was < 1
  all.others = df %>%
    filter(
      !(Name %in% repeats$Name)
    ) %>%
    group_by( Name, MZ, RT.min) %>%
    summarise()
  
  # (3) Bind the two data frames together to yeild a df of only the metabolite annotations we care about
  df2 = rbind( divergent.repeats[,!names(divergent.repeats) %in% c("batch.median")], all.others )
  
  
  # (4) Inner joint the metabolite annotations we care about with the original data:
  #     Yields a data frame of ONLY IC values in the batch that match the annotations we want! Excludes all others.
  df3 = inner_join(df, df2)
  
  
  # (5) Lastly, now that we have dealt with the RT ranges >1, we address those RT ranges <1:
  #     Repeat entries with RT range <1 still remain from binding all.others.
  #     In these instances, I am making the assertion that this is simply RT drift
  #     Therefore, on a per-sample basis, for every metabolite NAME I will take the median MZ and RT values,
  #               and sum the Ion Count vales. 
  
  # ** still needs work; not the perfect approach because we do the grouping by samples
  # ** Technically, between samples we could have the same metabolite Name with different RTs between samples
  # **    ... (the spirit of this all is to have unified metabolite annotations)
  
  df4 = df3 %>%
    # mutate(
    #   Ion.Mode = NA,
    #   Ion.Mode = ifelse( grepl("^[P]_", Sequence.ID) == TRUE, "Pos", Ion.Mode),
    #   Ion.Mode = ifelse( grepl("^[N]_", Sequence.ID) == TRUE, "Neg", Ion.Mode),
    #   Sequence.ID = sub("^[PN]_", "", Sequence.ID)
    # ) %>%
    group_by(Raw.File, Sequence.ID, Name, Formula) %>%
    summarise(
      MZ = median(MZ),
      RT.min = median(RT.min),
      Ion.Count = sum(Ion.Count)
      #Ion.Mode = Ion.Mode
    )
  
  # metabs = unique( df4[, c("Name", "Ion.Mode")] )
  
  # length(unique(df3$Name))
  # length(unique(df4$Name))
  # Still same # metabolites but we removed redundancies
  
  
  
  return(df4)
}




# ========== Function: 'concat_RTs()'  ==========
# --
# scans data frame for near identical RT & MW entries and concatenates them into one entry according to the default or specified tolerances
# designed to operate on one data frame for a SINGLE sample. Wrap into a loop to operate on data frames for multiple samples
# ... e.g. write all individual sample data frames into a list of data frames and loop over the list, applying this function each iteration
concat_RTs = function(df, RT.diff = 0.5, MZ.diff = 0.00001){
  
  #cat(paste0( "Operating on df: ", deparse(substitute(df)), "\n"))
  
  df = arrange(df, MZ, RT.min ) %>% mutate( Ion.Count = ifelse(is.na(Ion.Count), 0, Ion.Count) )
  
  i=1
  while( i <= nrow(df)-1 ){
    
    #cat(paste0("on row ", i, " of ", nrow(df), "\n"))
    
    j = i+1
    
    while( j <= nrow(df) ){
      if( abs(df$MZ[j] - df$MZ[i]) < MZ.diff  && abs(df$RT.min[j] - df$RT.min[i]) < RT.diff ){
        
        df$Ion.Count[i] = sum(df$Ion.Count[i], df$Ion.Count[j])
        df = df[-c(j), ]
        
      }else{
        if( j < nrow(df) ){
          if( abs(df$MZ[j+1] - df$MZ[i]) < MZ.diff ){ j = j+1 }else{ break }
        }else{ break }
        
      }
    }
    i = i+1
  }
  return(df)
  
}





# ========== get_legend2() ==========
#cowplot get_legend() currently not working
get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


# ========== Z-OLD: 'rm_repeat_names()' ==========
# --
# Function: 'rm_repeat_names()'
# Retains the repeat entry with the highest average value across ALL samples (only addresses repeats with a "Name")
# ...designed to operate on a data frame containing ALL samples
Z_OLD_rm_repeat_names = function(df){
  
  avgs = df %>%
    mutate(
      # set NAs to zero to properly compute mean
      value = ifelse(is.na(Ion.Count), 0, Ion.Count)
    ) %>%
    group_by( Name, Calc.MW, RT.min) %>%
    summarise(
      metab_average = mean(value)
    ) %>%
    group_by(Name) 
  
  
  top_avgs = rbind(
    avgs %>% top_n(1, metab_average),
    subset(avgs, Name == "")
  )
  
  df = inner_join( df, top_avgs) %>% select(-c("metab_average"))
  
  return(df)
}



# ========== Z-OLD: Function: 'select_top_isomer'==========
# --
# *** Background:
# This function improves on rm_repeat_names() but instead of removing the replicate name annotations altogether,
# it labels them as name : isomer 1, name : isomer 2, etc.
#
# *** Pre requisites:
# Run this function by group splitting the melted data into a list of batches: batch_list = exp_data  %>% group_by( Raw.File ) %>% group_split(.)
# lapply this function to the list of batches
#
# *** Logic:
# First, establish the annotated metabolites which have variant RTs across the batch
# We will assert that if the RT range for a metabolite in the batch is > 0.3 min then we need to determine which RT is most precise
# "divergent.repeats"
# Otherwise, if the RT range of a metabolite is <= 1, it is just RT drift. We will group these together and sum these
# "all.others"

Z_select_top_isomer = function(df){
  
  # df= batch_list[[1]] # Debug line
  
  
  # (1) Establish the metabolite names that are > 0.3 RT range across the batch
  repeats = df %>%
    group_by( Name ) %>%
    summarise(
      rt.range = max(RT.min, na.rm = TRUE)-min(RT.min, na.rm = TRUE),
      # range.mzcloud = max(mzcloud.conf, na.rm = TRUE) - min(mzcloud.conf, na.rm = TRUE)
      ## ...moving forward i'm going to ignore the mzcloud confidence for these data bc it doesnt seem to be variant between diff RTs of the same NAME
    ) %>% 
    filter(
      rt.range > 0.3 # we need to understand how changing this affects the output
      # when > 1.0 we have 485 metabolites; when > 0.3 we have 641 metabolites;
      # there are 1070 unique metabolites in the dataset so about half have repeat/isomer annotations
    ) 
  # these are the only Names to which we need to apply selective criteria
  
  
  
  # (2) Select only the ICs where the metabolite range > 0.3,
  #     and determine the [Name, MZ, RT.min] with the highest median IC of the batch
  #     Label the rest with " : isomer n"
  divergent.repeats = df %>%
    filter(
      Name %in% repeats$Name,
      is.na(Ion.Count) == FALSE
    ) %>%
    mutate(
      value = Ion.Count
    ) %>%
    group_by( Name, MZ, RT.min) %>% # Multiple RT/MZ values in the batch are associated with the same name. For each, calculate the median
    summarise(
      
      # There is probably an even more robust set of selection criteria we can make depending on how the IC data are distributed
      
      # ** consider the delta mass
      # ** if second highest median is far away, keep both but label the metab as #2
      
      # For now we are just going to pick the max median
      batch.median = median(value),
      
      
      # min= min(value),
      # max = max(value),
      # IQR = format(IQR(value),  scientific = TRUE, digits = 2),
      # var.mean.coef = sd(value, na.rm = TRUE)/mean(value),
      # shapiro.test = shapiro.test(log10(value))$p.value,
      # count = n()
      
    ) %>%
    group_by(Name) %>%
    mutate(
      rank = rank(-batch.median),
      Name.isomer = ifelse(rank != 1, paste0(Name, " : isomer ", rank-1), Name),
      rank = NULL
    )
  
  
  # (3) Make a reduced data frame of all other metabolites where the RT range was < 0.3
  all.others = df %>%
    filter(
      !(Name %in% repeats$Name)
    ) %>%
    mutate( Name.isomer = Name) %>%
    group_by( Name, Name.isomer, MZ, RT.min) %>% # Name and Name.isomer are the same for this df. Just need the .isomer column to bind to divergent.repeats
    summarise()
  
  # (3) Bind the two data frames together to yeild a df of only the metabolite annotations we care about
  df2 = rbind( divergent.repeats[,!names(divergent.repeats) %in% c("batch.median")], all.others )
  
  
  # (4) Inner joint the metabolite annotations we care about with the original data.
  #  ... before in rm_repeat_names this would have reduced the # of rows vs exp_data2
  # ... but now since were just labeling repeats as isomer, we expect equal # of rows
  df3 = inner_join(df, df2)
  
  
  # (5) Lastly, now that we have dealt with the RT ranges > 0.3, we address those RT ranges < 0.3:
  #     Repeat entries with RT range < 0.3 still remain from binding all.others.
  #     In these instances, I am making the assertion that this is simply RT drift
  #     Therefore, on a per-sample basis, for every metabolite NAME I will take the median MZ and RT values,
  #               and sum the Ion Count vales. 
  
  # ** still needs work; not the perfect approach because we do the grouping by samples
  # ** Technically, between samples we could have the same metabolite Name with different RTs between samples
  # **    ... (the spirit of this all is to have unified metabolite annotations)
  
  df4 = df3 %>%
    group_by(Raw.File, Sequence.ID, Name.isomer, Formula) %>%
    summarise(
      MZ = min(MZ),
      RT.min = min(RT.min),
      Ion.Count = sum(Ion.Count),
      #Ion.Mode = Ion.Mode
    ) %>%
    ungroup() %>%
    mutate(
      Name = Name.isomer,
      Name.isomer = NULL
    ) %>%
    select(Raw.File, Sequence.ID, Name, Formula, MZ, RT.min, Ion.Count)
  
  # t3 = subset(df4a, grepl(": isomer", df4a$Name.isomer) == FALSE )
  # metabs = unique( df4[, c("Name", "Ion.Mode")] )
  
  # length(unique(df3$Name))
  # length(unique(df4$Name))
  # Still same # metabolites but we removed redundancies
  
  
  
  return(df4)
}
