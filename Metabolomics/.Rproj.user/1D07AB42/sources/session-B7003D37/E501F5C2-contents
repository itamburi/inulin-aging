# ****** Inulin Aging Mouse Tissue (Liver and Muscle)
# ****** Spring 2026 - Ian T, Jessie K, Randa L


# 
here::i_am("code/01 Process tissue metabolomics.R")
source(here::here("code/z-source.R"))
#

# ========== 1.0 - Convert xlsx to csv ==========
# --
here::i_am("code/01 convert xlsx to csv.R")

library(tidyverse)
library(here)
library(readxl)

# Loop through each .xlsx file and convert to .csv
CD_files <- list.files(here("data/raw/compound discoverer"), pattern = "\\.xlsx$", full.names = TRUE)

for (CD_file in CD_files) {

  file_name = tools::file_path_sans_ext(basename(CD_file))   # Extract the file name without extension

  data = readxl::read_xlsx(CD_file)   # Read .xlsx file
  
  csv_file = file.path(here("data/raw"), paste0(file_name, ".csv"))   # Set the output .csv file path

  write.csv(data, csv_file, row.names = FALSE)  # Write data to .csv file
  cat("Converted", CD_file, "to", csv_file, "\n")
}

# ========== 2.0 - Format Compound Disvoverer Data ==========
# --

# List all CD .csv files in the raw directory
csv_files = list.files(here("data/raw/compound discoverer"), pattern = "\\.csv$", full.names = TRUE)

# Loop through each .csv file and read into a list of data frames
csv_data_list = list()
for (csv_file in csv_files) {
  csv_data = read.csv(csv_file)
  csv_data_list[[csv_file]] = csv_data
}


# Function to tidy each data frame. Remove/rename columns, and reshape the structure
tidy_df_list <- function(filename) {
  
  #tidy = csv_data_list[[1]] %>%
  tidy = csv_data_list[[filename]] %>%
    select(
      # remove columns containing these strings....
      #-contains("Annot..Source"),
      -c("Annot..Source..Predicted.Compositions","X..mzCloud.Results", "mzCloud.Best.Match", "mzCloud.Best.Sim..Match"),
      -contains("Mass.List."), -contains("X..Adducts"), -contains("MS2"), -contains("Reference.Ion")
      # With all these parameters in place, we end up retaining just the columns indicated in the pivot_longer() line
    ) %>% 
    # delete these strings from columns containing them
    rename_with(~str_remove(., "Area..")) %>% 
    rename_with(~str_remove(., ".raw..*")) %>%
    rename(
      mz = m.z,
      calc.mw = Calc..MW,
      rt.min = RT..min.,
      annot.source.mzcloud = Annot..Source..mzCloud.Search,
      annot.source.masslist = Annot..Source..MassList.Search,
      mzcloud.conf = mzCloud.Best.Match.Confidence,
      delta.mass.ppm = Annot..DeltaMass..ppm.
    ) %>%
    # melt to ion count observation per row
    pivot_longer(.,
                 cols = !c("Name", "Formula", "calc.mw", "mz", "rt.min", "annot.source.masslist","annot.source.mzcloud","mzcloud.conf", "delta.mass.ppm"),
                 names_to = "sequence.id",
                 values_to = "ion.count") %>%
    mutate(
      Formula = gsub(" ","", Formula),
      raw.file = gsub(".*/", "", filename) # make a new column referencing the source CD file
    ) %>%
    rename_with(tolower)
  
  print(paste0("Completed ", filename))
  
  return(tidy)
  
}

# for every dataframe within the list of dfs, lapply() the tidy function above
tidy_list = lapply(names(csv_data_list), tidy_df_list)
tidy_df = bind_rows(tidy_list)


# We realized that LDL15_2 should actually be LDL16. Previously we just excluded it in 06. Now we will make the change
t = tidy_df %>%
  filter(grepl("ldl15_2",sequence.id)) %>%
  group_by( sequence.id ) %>%
  summarise(n())

# change P_CS_060_ldl15_2 to ldl16
tidy_df$sequence.id <- ifelse(tidy_df$sequence.id == "P_CS_060m_ldl15_2", "P_CS_060m_ldl16", tidy_df$sequence.id)

# save all features:
write.csv(tidy_df, here("data/processed/ion counts/raw CD observations combined and melted.csv"), row.names = FALSE)

# Apply initial filtering
tidy_df2 = tidy_df %>%
  filter( 
    # (1) remove unannotated names
    is.na(name) == FALSE,
    # (2) remove similar to compounds
    grepl("\\[Similar to:", name) == FALSE,
    # (3) remove blanks
    grepl("_blank", sequence.id) == FALSE,
    # (4) Remove delta mass outside of +/- 10 ppm
    abs(delta.mass.ppm) <= 10
  )


# save annotated compounds:
write.csv(tidy_df2, here("data/processed/ion counts/filtered CD observations combined and melted.csv"), row.names = FALSE)

# *** note that these samples do not exist: ***
# S_060m_ldl08, E_000m_ldl16, E_060m_ldl11


# ========== 3.0 - Check Annotations Against Reference List, Apply mz Confidence Score Filtering ==========

# -- load annotated ion count data and make a df of unique compound annotations
expdata  = read.csv(here("data/processed/ion counts/filtered CD observations combined and melted.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate(name.exp = tolower(name))

cpds = expdata %>%
  distinct(raw.file, name, name.exp, formula, annot.source.mzcloud, annot.source.masslist, delta.mass.ppm, calc.mw, mz, rt.min, mzcloud.conf) %>%
  select(name, name.exp, everything())
length(unique(cpds$name))

# -- 1) Check against reference list.
# -- > If a name matches our reference list but the RT is > 3 min away from the reference list, we will flag it for removal
reference_list = read.csv(here("data/libraries/metabolites list_230424_ForCDprocessing.csv"), fileEncoding = "UTF-8-BOM") %>%
  #filter(is.na(RT) == FALSE) %>%
  rename_with(~ tolower(.)) %>%
  select( "hmdb.name", "rt", "formula") %>%
  rename( c("rt.ref" = "rt", "name.ref" = "hmdb.name")) %>%
  mutate( name.ref = trimws(name.ref) )

intersect = inner_join( cpds, reference_list, by = c("name.exp" = "name.ref"), keep = TRUE) %>%
  select("name.exp", "name.ref", "mz", "rt.min", "rt.ref") %>%
  mutate(
    rt.diff = abs(rt.min - rt.ref),
    rt.match = ifelse( rt.diff < 3, TRUE, FALSE)
  )

length(unique(subset(intersect, rt.match == TRUE)$name.exp)) # 792*
length(unique(subset(intersect, rt.match == FALSE)$name.exp)) # 52*

# -- These are the CD metabolties that matched our reference list but the RT was > 3 minutes different. We will remove these! 
discard = intersect %>%
  filter( rt.match == FALSE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(discard$name.exp))
# write.csv(discard, here("data/processed/discarded/summary of annotations removed from expdata due to non-matching RTs.csv"), row.names = FALSE )

# -- this is a summary of the CD metabolite-RT combinations that DID match our reference list. We will keep these!
keep = intersect %>%
  filter( rt.match == TRUE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(keep$name.exp))

# -- anti_join() will remove the observations in 'discard' from 'expcpds'
cpds2 = anti_join(cpds, discard) %>%
  left_join(., keep[,c("name.exp","mz","rt.min","rt.match")])

# -- mz cloud filtering among the cases that didnt match our ref list
mz_filt = cpds2 %>%
  filter(
    is.na(rt.match), # if rt.match == TRUE we keep regardless of mzcloud.conf. FALSE cases were already removed
    mzcloud.conf < 40
  ) 

# -- anti_join() will remove the observations from 'mz_filt'
cpds3 = anti_join(cpds2, mz_filt)


# -- Finally, Keep only the full expdata with annotations in cpds3
expdata2 = inner_join(expdata, cpds3)
length(unique(expdata2$name))
setdiff(expdata$name, expdata2$name)

write.csv(expdata2, here("data/processed/ion counts/CD data filtered by mismatches to ref and mz cloud conf.csv"),row.names = FALSE)



# ========== 4.0 - Resolve redundant compound entries / select top isomer ==========
# -- Manage repeated names with different RTs (but not yet consolidating by ion mode)

expdata = read.csv(here("data/processed/ion counts/CD data filtered by mismatches to ref and mz cloud conf.csv")) %>%
  mutate(
    cohort = case_when(
      grepl("ldl0[0-9]|ldl10", sequence.id) ~ "3mo",
      grepl("ldl[11-20]", sequence.id) ~ "5mo"
    )
  )

unique(expdata$cohort)
expdata %>% filter(is.na(cohort)) %>% distinct(sequence.id) %>% pull(sequence.id)

# --- How many replicate annotations for same name? ---
n = expdata %>%
  distinct(raw.file, name, mz, formula, rt.min) %>%
  group_by(raw.file, name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(n > 1)
# .. many metabolites still have multiple annotations per ion mode (ion mode is denoted by raw file)

# --- (A) Apply a function select_top_isomer() ---
# When multiple features share a compound name, select_top_isomer() selects the preferred peak and labels others as isomers.
# Because 3mo and 5mo cohorts were acquired in separate LC–MS sequences but processed together in Compound Discoverer,
# the function is run separately by cohort and ion mode to preserve within-batch consistency, and concordance between cohorts is assessed after.


# -- Make a list of data frames for each ion mode, per cohort. "raw.file" corresponds to posi/nega
batch_list =  expdata %>%
  filter(cohort == "3mo") %>%
  group_by( raw.file ) %>%
  group_split(.)

batch_3mo = lapply(batch_list, select_top_isomer3)

batch_list =  expdata %>%
  filter(cohort == "5mo") %>%
  group_by( raw.file ) %>%
  group_split(.)

batch_5mo = lapply(batch_list, select_top_isomer3)

# -- data with secondary isomers annotated
expdata2 = bind_rows(batch_3mo, batch_5mo)

# -- data frame with secondary isomers removed
expdata3 = expdata2 %>% filter(grepl(": isomer", name) == FALSE)

# -- consolidate the RT/MZ for the primary cluster and select the best ion count in the cluster for the sample
expdata4 = expdata3 %>%
  filter(!(name == "9-Decenoic acid" & formula == "C10H18N4O6S")) %>% # misannotation of 9-Decenoic acid that makes an issue later
  group_by(cohort, raw.file, sequence.id, name, formula) %>%
  arrange(rt.min, .by_group = TRUE) %>%   # ensures rows are ordered by rt.min
  summarise(
    rt = first(rt.min),           # == representative_rt for the primary cluster
    mz = first(mz),                   # corresponding MZ for the first RT
    ion.count = max(ion.count, na.rm = TRUE)  # Largest annotation PER SAMPLE
  ) %>%
  ungroup()


# --- (B) Evaluate the annotation concordance between cohorts ---
# * note mz can be different bw cohorts, but here we are evaluating by RT consensus *
anno = expdata4 %>%
  distinct(raw.file, name, rt, cohort) %>% 
  pivot_wider(names_from = cohort, values_from = rt) %>%
  mutate(
    diff =  abs(`3mo` - `5mo`)
  ) %>%
  filter(
    diff > 0
  ) %>%
  arrange(desc(diff))
# ... 363 annotations did not align between cohorts; select_top_isomer() chose different annotations
# After manual checking of this list and examining the chromatograms produced in '00b Isomer analysis.R' proceeding as follows:
# - RT difference of < 1.0 min between batches are reasonable RT drift and these can be kept/merged
# - RT diff > 1.0 lack consensus and so are unreliable -> remove
# - nonsense names, isomers from cmpd disc. will be removed
# ... remove:
isomers = anno %>% filter(grepl("isomer",name)) %>% pull(name) 
nonsense = c("4,7-dimethylpyrazolo[5,1-c][1,2,4]triazine-3-carbonitrile","4-(1H-1,2,4-triazol-1-yl)aniline","4,7-dimethylpyrazolo[5,1-c][1,2,4]triazine-3-carbonitrile",
             "5-methyl-1-(1,3,5-trimethyl-1h-pyrazol-4-yl)-1h-1,2,3-triazole","n6-(delta2-isopentenyl)-adenine")
drift = anno %>% filter(diff > 1.0) %>% pull(name)

remove = anno %>% filter( name %in% c(isomers, nonsense, drift) )
keep = anno %>% filter( !name %in% c(isomers, nonsense, drift) )

expdata5 = expdata4 %>%
  anti_join(., remove) %>%
  group_by(raw.file, name) %>%
  mutate(
    # Track min RT and corresponding m/z per group to flag residual duplicates
    min_rt = min(rt, na.rm = TRUE), # among RT drift <1, use the min RT and corresponding MZ
    min_mz = mz[which.min(rt)],
    # *** intitially we did the following to unify the RT drift < 1 MZ and RT values. But its better data hygine to keep the original values. We can come back later and unify if needed
    # rt = if_else(name %in% keep, min_rt, rt),
    # mz = if_else(name %in% keep, min_mz, mz)
  ) %>%
  ungroup()


# --- (C) check the number of occurences of each metabolite Name ---
t = expdata5 %>%
  group_by( name, raw.file, formula, min_mz, min_rt ) %>%
  summarise(median.ic = median(ion.count))
n = t %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
# ...Should 1 or 2 for everything, depending on n ion modes for the metabolite, if we did everything correctly

t = expdata5 %>%
  group_by( name, raw.file, formula, mz, rt ) %>%
  summarise(median.ic = median(ion.count))
n = t %>%
  group_by(name) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(drift = ifelse(name %in% keep$name, "RT drift < 1min bw cohorts",NA))
# ...any n>2 should only be due to small RT drift between the cohorts that we decided to keep anyways

names(expdata5)

expdata6 = expdata5 %>% select(cohort, raw.file, sequence.id, name, formula, rt, mz, ion.count)

# save dereplicated data
write.csv(expdata6, here("data/processed/ion counts/annotated and dereplicated CD obs without repeat isomers.csv"), row.names = FALSE)
write.csv(remove, here("data/processed/ion counts/names removed due to cohort RT drift or nonsense name.csv"), row.names = FALSE)

# ========== 5.0 - Consolidate Negative and Positive Ion Modes ==========
# -- Among metabolites detected in +&- modes, select the annotation that has the highest median ion count

# -- read in dereplicated data
expdata = read.csv(here("data/processed/ion counts/annotated and dereplicated CD obs without repeat isomers.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate(
    ion.mode = case_when(
      raw.file == "LDL1_2_neg_230831.csv" ~ "nega",
      raw.file == "LDL1_2_pos_230913.csv" ~ "posi"
    ),
    sample = gsub("^N_|^P_","",sequence.id)
  ) 

# -- confirm the n annotations per mode as in previous step
n_modes = expdata %>%
  distinct( name, ion.mode, formula, mz, rt ) %>%
  group_by(name, ion.mode) %>%
  mutate(n.per.mode = n()) %>%
  ungroup() %>%
  arrange(desc(n.per.mode), name, ion.mode, rt)
# anything with n.per.mode == 2 had an acceptable RT drift bw the cohorts.
# anything greater than 2 would be an error that we would go back and address in prior step
# 2025-12-28: following major reprocessing in step 4 we are now good to proceed!

# -- Determine the median of every metabolite 'Name' per ion mode
ic_medians = expdata %>%
  group_by( name, ion.mode, formula) %>%
  summarise(median.ic = median(ion.count)) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(n = n())

# -- Extract the annotation with the larger median across all the samples
# ** notice how we are not grouping by cohort; we want the ion modes per metabolite to match across cohorts
top_mode = ic_medians %>% group_by(name) %>% top_n(1, median.ic) 
length(unique(top_mode$name))
nrow(top_mode)
# rows == # names, we know we selected one each

# -- Select from expdata only the top median annotations
expdata2 = inner_join(expdata, top_mode %>% select(-median.ic, -n) )

# -- check that we only have 1 mode per metabolite
n_modes2 = expdata2 %>%
  distinct( name, ion.mode, formula) %>%
  group_by(name) %>%
  summarise(n.obs = n())


names(expdata2)
expdata3 = expdata2 %>% select(cohort, raw.file, sequence.id, ion.mode, sample, name, formula, mz, rt, ion.count)

# -- save the data consolidated by ion mode
write.csv(expdata3, here("data/processed/ion counts/annotated and processed CD obs consolidated by top ion mode.csv"), row.names = FALSE)
#


# ========== 5.0 - Merge and Consolidate CD and Maven Observations ==========
# --

# -- Skip to (B) if (A) is complete
# -- (A) Build Maven data frames and dfs of unique maven and CD annotations --
mav_files1 = list.files(here("../../Metabolites/Serum/Maven/serum data"), pattern = "\\.csv$", full.names = TRUE)
mav_files2 = list.files(here("../../Metabolites/Serum/Maven/BA"), pattern = "\\.csv$", full.names = TRUE)
mav_files3 = list.files(here("../../Metabolites/Serum/Maven/Bicine"), pattern = "\\.csv$", full.names = TRUE)
mav_files = c(mav_files1, mav_files2, mav_files3)

csv_data_list = list()

for (csv_file in mav_files) {
  
  #csv_file = mav_files[1] # Debug line
  csv_data = read.csv(csv_file)
  filename = sub(".*/", "", csv_file)
  
  csv_long = csv_data %>%
    filter(
      !adductName == "" # BA data form animal 11 had isotopologues whcih was creating erros; this makes sure we only extract the parent ions and ignore isotopologues
    ) %>% 
    pivot_longer(., cols = grep(".*_ldl\\d+", names(csv_data), value = TRUE), names_to = "sequence.id", values_to = "ion.count") %>%
    select("sequence.id", "compound", "formula", "medRt", "medMz","ion.count") %>%
    mutate(
      sample = sub("^[PN]_", "", sequence.id),
      ion.mode = ifelse( grepl("_p", filename) != TRUE | grepl("_BA", filename) == TRUE, "nega", "posi"),
      raw.file = filename
    )
  
  csv_data_list[[csv_file]] = csv_long
}

# -- Combine and save all maven ion count data
mav_obs = bind_rows(csv_data_list) %>%
  mutate(
    name = compound,
    sequence.id = ifelse(sequence.id == "P_A_060m_ldl18_","P_A_060m_ldl18", sequence.id),
    sample = ifelse(sample == "A_060m_ldl18_","A_060m_ldl18", sample)
  ) %>%
  rename_with(tolower) %>%
  select( raw.file, sequence.id, ion.mode, sample, name, formula, medrt, medmz, ion.count )
write.csv(mav_obs, here("data/processed/ion counts/maven merger/all maven observations.csv"), row.names = FALSE)

# -- Create a data frame of the unique Maven compounds
mav_compounds = mav_obs %>%
  group_by(ion.mode, name, formula) %>%
  summarise(
    source = "Maven",
    med.ic = median(ion.count),
    rt = median(medrt),
    mz = median(medmz),
    cv = sd(ion.count, na.rm = TRUE)*100/mean(ion.count, na.rm=TRUE),
    n.samples = n()
  ) %>%
  select(source, ion.mode, name, formula, mz, rt, med.ic, n.samples, cv)

# -- Create a data frame of the unique CD compounds
cd_compounds = read.csv(here("data/processed/ion counts/annotated and processed CD obs consolidated by top ion mode.csv")) %>%
  group_by(ion.mode, name, formula, mz, rt.min) %>%
  summarise(
    source = "CD",
    med.ic = median(ion.count),
    n.samples = n(),
    cv = sd(ion.count, na.rm = TRUE)*100/mean(ion.count, na.rm=TRUE)
  ) %>%
  select(source, ion.mode, name, formula, mz, rt.min, med.ic, n.samples, cv) %>%
  rename(rt = "rt.min")

names(mav_compounds)
names(cd_compounds)

all_cpds = rbind(mav_compounds, cd_compounds)
write.csv(all_cpds, here("data/processed/ion counts/maven merger/unique maven and cd compounds before consolidation.csv"), row.names = FALSE)


# (B) -- BY HAND/MANUAL -- Compare Maven and CD compounds  --

# *** At this point, we combine the dataframes of unique compounds from CD and Maven in Excel (dfs produced above)
# *** In excel, we note the source (MAV or CD) then sort the dataframe by RT, then formula
# *** We manually go through by row and select the preferred annotation between MAV and CD
# *** Usually we defer to the MAVEN annotation. However, there is some discretion involved
# *** In the end, we select the annotations to delete or rename as indicated by "remove" and "rename" columns

# *** We worked through this process for this LDL data initially with SJ, HB, and I.T. on 10/18/23 
# *** This document was saved as "how to replace metabolites in CD to Maven - final action_240108.csv"
# *** For data hygiene purposes, I was consolidated the ion count processing scripts into one script,
# *** ... and during this I repeated this the MAV/CD manual consolidation process for a second time on 10-21-2025
# *** ... I compared the annotations with the original "final action_240108" file to confirm that I was selecting the same annotations as before


# *** Having completed this manual compariosn, we can move to the next step of merging the annotations


# -- load the CD and Maven complete data and combine --
cd_data = read.csv(here("data/processed/ion counts/annotated and processed CD obs consolidated by top ion mode.csv")) %>%
  mutate(
    source = "CD"
  )

mav_data = read.csv(here("data/processed/ion counts/maven merger/all maven observations.csv"), fileEncoding = "UTF-8-BOM") %>%
  rename(
    rt = medrt,
    mz = medmz
  ) %>%
  mutate(
    source = "Maven"
  ) %>%
  group_by(name, formula) %>%
  mutate(
    # recall that for maven-picked compounds, maven provides a unique rt and mz per sample
    rt = median(rt),
    mz = median(mz),
    ion.mode = case_when(ion.mode =="N" ~ "nega", ion.mode == "P" ~ "posi"),
    cohort = case_when(
      grepl("ldl0[0-9]|ldl10", sequence.id) ~ "3mo",
      grepl("ldl[11-20]", sequence.id) ~ "5mo"
    )
  )

x = names(cd_data)
y = names(mav_data)
setdiff(union(x,y), intersect(x,y))

expdata = rbind(cd_data, mav_data)
names(expdata)

# -- Load the hand-picked comparisons of metabolites to delete/rename from I.T. 2025-10-21:
action = read.csv(here("data/processed/ion counts/maven merger/final maven and cd compounds to consolidate IJT 2025-10-21.csv"), fileEncoding = "UTF-8-BOM") %>%
  select(!c("mz","rt")) %>% # recall maven produces unique mz/rts per sample so we wont merge on these numbers
  mutate(
    ion.mode = case_when(ion.mode =="N" ~ "nega", ion.mode == "P" ~ "posi"),
    remove = ifelse(remove =="",NA,remove),
    rename = ifelse(rename =="",NA,rename)
  )

# -- Merge acitons with expdata, rename and filter
expdata2 = expdata %>%
  left_join(action) %>%
  mutate(
    name = ifelse( is.na(rename) == FALSE, rename, name)
  ) %>%
  filter(
    is.na(remove),
    !(name == "Glycine" & source == "CD" & ion.mode == "nega"),
    !(name == "Bicine" & source == "CD" & ion.mode == "nega"),
  ) %>%
  select(
    cohort, raw.file, source, sequence.id, sample, ion.mode, name, formula, mz, rt, ion.count
  )


# -- check the n annotations
n = expdata2 %>%
  distinct(name, rt, mz, cohort, source, ion.mode) %>%
  group_by(name, cohort) %>%
  mutate(n.annotations=n()) %>%
  arrange(desc(n.annotations), name, cohort)
# all 1 obs, good!
n2 = expdata2 %>%
  distinct(name, rt, mz) %>%
  group_by(name) %>%
  mutate(n.annotations=n()) %>%
  arrange(desc(n.annotations),name)
# anyhting == 2 was acceptable cohort RT drift <1 min from earlier

# -- save the consolidate CD + Mav data --
write.csv(expdata2, here("data/processed/ion counts/consolidated CD plus MAV obs.csv"), row.names = FALSE)


# ========== 6.0 - Remove Flagged Compounds, Add Metadata, and Perform Normalization ==========
# --
# === (A) Remove Flagged compounds ===

# -- Load Ion Count data
expdata = read.csv(here("data/processed/ion counts/consolidated CD plus MAV obs.csv"))
length(unique(expdata$name))

# -- Load the compounds we flagged to remove from the dataset
flag = read.csv(here("data/processed/ion counts/consolidated CD plus MAV compounds flagged for filtering 2025-10-21.csv")) %>%
  select(name, flag) %>%
  mutate(
    flag = case_when(flag == "flag" ~ "remove", flag == "" ~ "keep")
  )

expdata2 = expdata %>%
  left_join(., flag) %>%
  filter( flag != "remove" ) %>%
  select(-flag)


# === (B) Add Metadata ===

# -- make metadata
meta = expdata2 %>%
  distinct(sequence.id, sample) %>%
  filter(
    grepl("DDA_", sample) == FALSE,
    grepl("Blank", sample) == FALSE,
    grepl("ldl19_2", sample) == FALSE,
    grepl("ldl19_2", sample) == FALSE,
    grepl("N_B_000m_ldl14", sample) == FALSE,
  ) %>%
  separate(sequence.id, into = c("ion.mode","tissue", "time", "animal"), sep = "_", remove = FALSE) %>%
  mutate(
    # add cohort 3mo v 5mo
    cohort = case_when(
      grepl("ldl0[0-9]|ldl10", animal) ~ "3mo",
      grepl("ldl[11-20]", animal) ~ "5mo"
    ),
    # HL was only collected in pigs 9 and 10
    tissue = gsub("Hl", tissue, replacement = "HL"), 
    # H and HS are effectively the same. Change all annotations to HS
    tissue = gsub("Hs", tissue, replacement = "HS"),
    tissue = gsub("\\bH\\b", tissue, replacement = "HS"),
    
    ion.mode = case_when(ion.mode =="N" ~ "nega", ion.mode == "P" ~ "posi"),
  )

# -- Add metadata to expdata
expdata3 = inner_join(expdata2, meta)



# === (C) Combine labeled Valine ICs with expdata ===
# ** see Z-OLD section at end of this script for old assembly of Valine internal standard data **
# -- recall each cohort used a different labeled standard (N15 v C13N15) and each was run in different batches

valic = read.csv(here("data/raw/maven valine/consolidated valine internal standard data.csv")) %>%
  separate(sequence.id, into =c("ion.mode","sample"), sep="_", extra = "merge", remove = FALSE) %>%
  mutate(
    ion.mode = case_when(ion.mode =="N" ~ "nega", ion.mode == "P" ~ "posi")
  ) %>%
  # Note this file contians "N_S_000m_ldl19" & N_S_000m_ldl19_2". Using the first inner_join with metadata
  inner_join(., meta) %>%
  # mark potential valine outliers
  group_by(ion.mode, tissue, cohort) %>%
  mutate(
    Q1 = quantile(val.ic, 0.25, na.rm = TRUE),
    Q3 = quantile(val.ic, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is.outlier = val.ic < (Q1 - 1.5 * IQR) | val.ic > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup() %>%
  distinct(sequence.id, animal, is.outlier, val.ic)

expdata4 = left_join(expdata3, valic)

z = expdata4 %>% filter(is.na(val.ic)) %>% distinct(sequence.id,val.ic)
# *** 4 OF THESE STILL HAS A NA for VALINE


plot_df <- expdata4 %>%
  distinct(ion.mode, sample, tissue, time, animal, cohort, val.ic, is.outlier) %>% 
  tidyr::replace_na(list(is.outlier = FALSE)) %>%    # in case any NAs
  filter(!is.na(val.ic))                              # drop missing values

ggplot(plot_df, aes(tissue, val.ic)) +
  facet_wrap(~ion.mode) +
  geom_jitter(aes(color = cohort), width = 0.2, alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(.x, is.outlier),
    aes(label = animal),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  labs(
    title = "Labeled valine ion abundances across tissues and cohorts",
    x = "Tissue", y = "Valine ion count (val.ic)"
  ) +
  theme_bw()

# -- compute valine normalization
expdata5 = expdata4 %>%
  group_by(ion.mode, cohort, tissue) %>%
  mutate(
    med.val = median(val.ic, na.rm = TRUE),
    ic.val.norm = (ion.count /val.ic ) * med.val
  ) %>%
  ungroup() %>%
  mutate(
# *** need see if we can impute values based on sequence position? for now this works
    ic.val.norm = ifelse( is.na(val.ic), ion.count, ic.val.norm)
  )

# -- compute median normalization
expdata6 = expdata5 %>%
  group_by(tissue, time, cohort, name) %>% 
  mutate(
    med.metab = median(ion.count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    `log10(ic/med.metab)` = log10(ic.val.norm/med.metab)
  )%>%
  group_by(sample) %>%
  # is it necessary to group the norm also by ion mode here?
  # I would argue No since the distributions of pos and neg are in the same range
  mutate(
    med.sample = median(`log10(ic/med.metab)`, na.rm = TRUE),
    ic.med.norm = ( 10^(`log10(ic/med.metab)` - med.sample)*med.metab )
  ) %>%
  ungroup() %>%
  mutate(
    # *** need to address
    ic.med.norm = ifelse( ic.med.norm == "NaN", ic.val.norm, ic.med.norm)
  ) %>%
  select(
    sequence.id, sample, ion.mode, source, name, formula, mz, rt, tissue, time, animal, cohort, ion.count, ic.val.norm, ic.med.norm,
    is.outlier
  ) %>%
  rename(ic.raw = "ion.count")

  

# (E) === PCA compare raw IC vs normalized / basic QA ===

# pivot wide for raw ion.count
icmat_raw <- expdata6 %>%
  select(sample, name, ic.raw) %>%
  pivot_wider(
    names_from = sample,
    values_from = ic.raw,
    values_fill = 0
  ) %>%
  column_to_rownames(var="name") %>%
  log1p() %>%  # log(1 + x)
  t()          # transpose so rows = samples

# pivot wide for val.adj.ic
icmat_val <- expdata6 %>%
  select(sample, name, ic.val.norm) %>%
  pivot_wider(
    names_from = sample,
    values_from = ic.val.norm,
    values_fill = 0
  ) %>%
  column_to_rownames(var="name") %>%
  log1p() %>%  # log(1 + x)
  replace(is.na(.), 0) %>% 
  t()          # transpose so rows = samples

# pivot wide for val.adj.ic
icmat_med <- expdata6 %>%
  select(sample, name, ic.med.norm) %>%
  pivot_wider(
    names_from = sample,
    values_from = ic.med.norm,
    values_fill = 0
  ) %>%
  column_to_rownames(var="name") %>%
  log1p() %>%  # log(1 + x)
  replace(is.na(.), 0) %>% 
  t()          # transpose so rows = samples

# run pca
pca_raw <- prcomp(icmat_raw, scale. = TRUE, center = TRUE)
pca_val <- prcomp(icmat_val, scale. = TRUE, center = TRUE)
pca_med <- prcomp(icmat_med, scale. = TRUE, center = TRUE)

# convert to data frame with metadata
plot_pca <- function(pca_obj, title = "PCA") {
  
  pca_obj = pca_val
  
  mm = meta %>% distinct(sample, tissue, time, animal, cohort)
  
  scores <- as.data.frame(pca_obj$x[, 1:2]) %>%
    rownames_to_column(var="sample") %>%
    left_join(., mm, by = "sample")
  
  # centroids + counts per tissue
  tissue_centroids <- scores %>%
    group_by(tissue) %>%
    summarise(
      PC1 = median(PC1, na.rm = TRUE),
      PC2 = median(PC2, na.rm = TRUE),
      n   = n(),
      .groups = "drop"
    )
  
  scores %>%
  ggplot(., aes(x = PC1, y = PC2, color = tissue, shape = cohort)) +
    geom_point(size = 3, alpha = 0.85) +
    # optional: 68% ellipses per tissue
    # stat_ellipse(aes(group = tissue, color = tissue),
    #              linetype = "dashed", alpha = 0.5, linewidth = 0.6) +
    # #centroid labels (no inherited aes -> no 'cohort' required)
    # geom_text(
    #   data = tissue_centroids,
    #   mapping = aes(x = PC1, y = PC2, label = paste0(tissue, " (n=", n, ")")),
    #   inherit.aes = FALSE,
    #   color = "black", fontface = "bold", size = 4, vjust = -0.8
    # ) +
    #xlim(-50,50) + ylim(-20,20) +
    theme_classic(base_size = 14)

  }

plot_pca(pca_raw, title = "PCA of ion.count")
plot_pca(pca_val, title = "PCA of val.adj.ic")
plot_pca(pca_med, title = "PCA of med.adj.ic")

# generally normal no changes in variation pattern bw raw and normalized
# Urine are distinct from blood samples
# Strongest pattern of variation among blood appears to be cohort
# Note a handful of samples are very different from others eg. PC1 -200


# check that we didnt produce any duplicate annotations once more
t = expdata6 %>%
  ungroup() %>%
  distinct(name, mz) %>%
  group_by(name) %>%
  summarise(n())

ggplot( expdata6, aes(y= (log10(ic.raw)))) +
  geom_line()
  geom_histogram(bins = 100) +
  facet_wrap(~Animal)

# -- correlation
expdata6 %>%
    summarise(
      cor_raw_val = cor(ic.raw, ic.val.norm, use = "complete.obs"),
      cor_raw_med = cor(ic.raw, ic.med.norm, use = "complete.obs"),
      cor_val_med = cor(ic.val.norm, ic.med.norm, use = "complete.obs")
    )

# -- scatter
ggplot(expdata6, aes(x = ic.raw, y = ic.val.norm)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(x = "Raw ion count", y = "Value-normalized ion count")

# -- fold change distribution
fc = expdata6 %>%
  mutate(
    log2fc_val = log2(ic.val.norm / ic.raw),
    log2fc_med = log2(ic.med.norm / ic.raw)
  ) %>%
  pivot_longer(starts_with("log2fc"), names_to = "norm_type", values_to = "log2fc")

ggplot(fc, aes(x = log2fc, fill = norm_type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()

outlier = fc %>%
  filter( abs(log2fc) > .5, norm_type == "log2fc_val") %>%
  select(sequence.id, source, name, tissue, time, animal, cohort, ic.raw, ic.val.norm, log2fc)

z = outlier %>% group_by(sequence.id) %>% summarise(n())

# * check median ion for that sample or for those maven samples

unique(outlier$sequence.id)


# -- global intensity shift
expdata6 %>%
  group_by(sample) %>%
  summarise(
    mean_raw = mean(ic.raw, na.rm = TRUE),
    mean_val = mean(ic.val.norm, na.rm = TRUE),
    mean_med = mean(ic.med.norm, na.rm = TRUE)
  ) %>%
  pivot_longer(-sample, names_to = "type", values_to = "mean_ic") %>%
  ggplot(aes(x = sample, y = mean_ic, color = type, group = type)) +
  geom_line() + geom_point() +
  scale_y_log10() + theme_minimal()


write.csv(expdata6, here("data/processed/ion counts/tidy and normalized ion count data.csv"), row.names = FALSE)


# *** Notes and observations
# > in the future:
# -- move the hematocrit normalization to this section and make a master IC df of all adjustments
# -- Examine PCA across all normalization approaches (raw, valine, median, hematocrit)
# -- calculate variation from fixed and random effects?

# --- wide format
library(RefMet)
try(httr::GET("https://www.metabolomicsworkbench.org/rest/refmet/"))
mycpds = unique(expdata6$name)
mapped = refmet_map_df(mycpds)

names(expdata6)
expdata_wide = expdata6 %>%
  select(sample, ion.mode, source, name, ic.med.norm) %>%
  pivot_wider(names_from = sample, values_from = ic.med.norm)
write.csv(expdata_wide, here("data/processed/ion counts/tidy ion count data wide 2025-12-28.csv"), row.names = FALSE)


#### 2.0 - convert to wide format for manuscript revisions, August 2024


exp = read.csv( here("data/processed/ion counts/finalized and median normalized observations.csv"), fileEncoding = "UTF-8-BOM") %>%
  select(
    Sample, Name, Formula,  RT, MZ, Adj.ic
  ) %>%
  pivot_wider(names_from = Sample, values_from = Adj.ic)

write.csv(exp, here("data/processed/ion counts/finalized and median normalized observations_wide.csv"), row.names = FALSE)





# ========== z-OLD Code ==========
# --


# === (C) Tidy Valine Internal standard data ===
# *** Note: 10-29-2025 ***
# --- Initially the valine internal standard data were aggregated from el-maven export files and processed with the code below
# --- The process was kind of confusing due to different internal standards by batch and due bc files were split up by animal/ion mode
# --- I actually ended up managing and getting it correct!
# --- ... However, I realized HB & JL already consolidated the maven internal standard data by hand in the file "valine_neg_pos.xlsx"
# --- I reorganized the consolidated information into the csv file "data/raw/maven valine/consolidated valine internal standard data.csv"
# --- I checked that the ion counts for each sequence.id were equal between the version I assembled and in the file above^^
# --- "consolidated valine internal standard data.csv" had a sample that was missing before!
# --- For simplicity, I ended up using "consolidated valine internal standard data.csv" in the processing, but am retaining the old code here:
# *** End note ***
#
# the valine data files were kind of messy so went in and manually edited them to deal with duplicates
# kept the original files in raw/maven valine/z-old

val_files = list.files(here("data/raw/maven valine"), pattern = "\\.csv$", full.names = TRUE)

csv_data_list = list()

for (csv_file in val_files) {
  
  # csv_file = val_files[8] # Debug line
  csv_data = read.csv(csv_file)
  filename = sub(".*/", "", csv_file)
  
  csv_long = csv_data %>%
    pivot_longer(., cols = grep("*_ldl\\d+", names(csv_data), value = TRUE), names_to = "sequence.id", values_to = "ion.count") %>%
    select("sequence.id", "compound", "isotopeLabel","formula", "medRt", "medMz","ion.count") %>%
    mutate(
      sample = sub("^[PN]_", "", sequence.id),
      ion.mode = ifelse( grepl("_p", filename) != TRUE | grepl("_BA", filename) == TRUE, "N", "P"),
      raw.file = filename
    ) %>%
    mutate_all(as.character)
  
  csv_data_list[[csv_file]] = csv_long
}

# === (D) Compute Valine and Median Normalization ===
# -- some of the pigs used N15 valine, some used C13N15 valine
valic = bind_rows(csv_data_list) %>%
  rename_with(tolower) %>%
  select(compound, isotopelabel, sequence.id, ion.count) %>%
  filter(
    isotopelabel %in% c("C13N15-label-5-1","N15-label-1"),
    # urine maven for two valine label types were done in the same file which makes this necessary:
    !(grepl("_U_",sequence.id) & ion.count == 0)
  ) %>%
  separate(sequence.id, into =c("ion.mode","sample"), sep="_", extra = "merge", remove = FALSE) %>%
  inner_join(., meta) %>%
  mutate(
    compound = "Valine-labeled",
    ion.count = as.numeric(ion.count),
    #ion.count = ifelse(ion.count == 0, NA, ion.count),
    val.ic = ion.count
  ) %>%
  # -- mark valine outliers
  group_by(ion.mode, tissue, cohort) %>%
  mutate(
    Q1 = quantile(val.ic, 0.25, na.rm = TRUE),
    Q3 = quantile(val.ic, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is.outlier = val.ic < (Q1 - 1.5 * IQR) | val.ic > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup() %>%
  distinct(sequence.id, is.outlier, val.ic)

