## Pig Tissue Cohort Processing - Winter 2026

here::i_am("code/02 Tissue and Mass Action/01 process tissue ion count data.R")
source(here::here("code/z-source.R"))
#

suppressPackageStartupMessages({
  library(readxl)
})

# ========== 0.0 - Enviornment ==========
# --
raw_cd_dir <- here("data/raw/tissue/CD")
raw_maven_dir <- here("data/raw/tissue/maven")
processed_dir <- here("data/processed/tissue/ion counts")
processed_maven_dir <- file.path(processed_dir, "maven merger")

first_existing_path <- function(paths) {
  hit <- paths[file.exists(paths)][1]

  if (is.na(hit)) {
    stop("Missing required file. Checked: ", paste(paths, collapse = ", "))
  }

  hit
}

clean_header <- function(x) {
  x %>%
    str_replace("^\ufeff", "") %>%
    str_replace_all("[^A-Za-z0-9]+", ".") %>%
    str_replace_all("\\.+", ".") %>%
    str_replace("^\\.", "") %>%
    str_replace("\\.$", "") %>%
    tolower()
}

read_metadata_csv <- function(path, ...) {
  df <- read.csv(path, check.names = FALSE, ...)
  cleaned_names <- clean_header(names(df))
  empty_idx <- which(is.na(cleaned_names) | cleaned_names == "")

  if (length(empty_idx)) {
    cleaned_names[empty_idx] <- paste0("unnamed.", empty_idx)
  }

  names(df) <- make.unique(cleaned_names, sep = ".")
  df
}

detect_cohort_label <- function(values) {
  cohort_labels <- values %>%
    str_extract("[A-Za-z]+\\d+$") %>%
    str_remove("\\d+$")

  cohort_labels <- unique(cohort_labels[!is.na(cohort_labels) & cohort_labels != ""])

  if (!length(cohort_labels)) {
    return(NA_character_)
  }

  if (length(cohort_labels) > 1) {
    warning("Multiple cohort labels detected: ", paste(cohort_labels, collapse = ", "), ". Using the first.")
  }

  cohort_labels[[1]]
}

# ========== 1.0 - Converd xlsx to csv ==========
# -- Convert the compound discoverer xlsx files to csv files


# Loop through each .xlsx file and convert to .csv
CD_files <- list.files(raw_cd_dir, pattern = "\\.xlsx$", full.names = TRUE)

for (CD_file in CD_files) {

  file_name = tools::file_path_sans_ext(basename(CD_file))   # Extract the file name without extension

  data = readxl::read_xlsx(CD_file)   # Read .xlsx file

  csv_file = file.path(raw_cd_dir, paste0(file_name, ".csv"))   # Set the output .csv file path

  write.csv(data, csv_file, row.names = FALSE)  # Write data to .csv file
  cat("Converted", CD_file, "to", csv_file, "\n")

}



# ========== 2.0 - Format Compound Discoverer Data ==========
# -- Format Compound Disc files into tidy long form data frames

# List all CD .csv files in the raw directory
csv_files = list.files(raw_cd_dir, pattern = "\\.csv$", full.names = TRUE)

# Loop through each .csv file and read into a list of data frames
csv_data_list = list()
for (csv_file in csv_files) {
  csv_data = read.csv(csv_file)
  csv_data_list[[csv_file]] = csv_data
}

# Function to tidy each data frame. Remove/rename columns, and reshape the structure
tidy_df_list <- function(filepath) {
  
  
  #tidy = csv_data_list[[2]] %>% # debug
  tidy = csv_data_list[[filepath]] %>%
    select(
      # remove columns containing these strings....
      #-c("X..mzCloud.Results", "mzCloud.Best.Match", "mzCloud.Best.Sim..Match"),
      #-contains("Annot..Source"),
      -c("Annotation.MW"),
      -contains("Mass.List."), -contains("X..Adducts"), -contains("MS2"), -contains("Reference.Ion"), -contains("Group.Area..")
      # With all these parameters in place, we end up retaining just the columns indicated in the pivot_longer() line
    ) %>% 
    # delete these strings from columns containing them
    rename_with(~str_remove(., "Area..")) %>% 
    rename_with(~str_remove(., ".raw..*")) %>%
    rename(
      MZ = m.z,
      Calc.MW = Calc..MW,
      RT.min = RT..min.,
      annot.source.predicted = Annot..Source..Predicted.Compositions,
      annot.source.mzcloud = Annot..Source..mzCloud.Search,
      annot.source.masslist = Annot..Source..MassList.Search,
      mzcloud.n.hits = X..mzCloud.Results,
      mzcloud.best.match = mzCloud.Best.Match,
      mzcloud.best.sim.match = mzCloud.Best.Sim..Match,
      mzcloud.conf = mzCloud.Best.Match.Confidence,
      delta.mass.ppm = Annot..DeltaMass..ppm.
      ) %>%
    # melt to ion count observation per row
    pivot_longer(
      cols = !c("Name", "Formula", "Calc.MW", "MZ", "RT.min", "annot.source.predicted","annot.source.mzcloud","annot.source.masslist",
                "mzcloud.n.hits", "mzcloud.best.match","mzcloud.best.sim.match","mzcloud.conf","delta.mass.ppm"),
      names_to = "Sequence.ID",
      values_to = "Ion.Count"
    ) %>%
    mutate(
      Formula = gsub(" ","", Formula),
      Raw.File = gsub(".*/", "", filepath) # make a new column referencing the source CD file
    ) %>%
    rename_with(~ tolower(.))
  
  print(paste0("Completed ", filepath))
  
  return(tidy)
  
}


# for every dataframe within the list of dfs, lapply() the tidy function above
tidy_list = lapply(names(csv_data_list), tidy_df_list)
tidy_df = bind_rows(tidy_list)

# save all features:
write.csv(tidy_df, file.path(processed_dir, "raw CD observations combined and melted.csv"), row.names = FALSE)

# Apply filtering
tidy_df2 = tidy_df %>%
  filter( 
    # (1) remove unannotated names
    is.na(name) == FALSE,
    name != "",
    # (2) remove similar to compounds
    grepl("\\[Similar to:", name) == FALSE,
    # (3) remove blanks
    grepl("_blank", sequence.id) == FALSE,
    # (4) Remove delta mass outside of +/- 10 ppm
    abs(delta.mass.ppm) <= 10
  )


names(tidy_df2)
unique(tidy_df2$name)
unique(tidy_df2$sequence.id)
unique(tidy_df2$raw.file)

# save annotated compounds:
write.csv(tidy_df2, file.path(processed_dir, "filtered CD observations combined and melted.csv"), row.names = FALSE)


# ========== 3.0 - Check Annotations Against Reference List, Apply Filtering ==========
# -- Correct mismatches with reference lists

# load annotated data
expdata  = read.csv(file.path(processed_dir, "filtered CD observations combined and melted.csv")) %>%
  mutate(name.exp = tolower(name))

expcpds = expdata %>%
  distinct(raw.file, name, name.exp, formula, annot.source.predicted, annot.source.mzcloud,annot.source.masslist, delta.mass.ppm, calc.mw, mz, rt.min, mzcloud.n.hits, mzcloud.best.match, mzcloud.conf) %>%
  select(name, name.exp, everything())
length(unique(expcpds$name))

# -- 1) Check against reference list.
# -- > If a name matches our reference list but the RT is > 3 min away from the reference list, we will flag it for removal
reference_list = read_metadata_csv( here("data/libraries/metabolites list_230424_ForCDprocessing.csv") ) %>%
  select( "cj", "hmdb.name", "rt", "formula", "detected.in") %>%
  rename( c("rt.ref" = "rt", "name.ref" = "hmdb.name")) %>%
  mutate(
    name.ref = trimws(name.ref)
  ) %>%
  filter(
    !(name.ref == "gentisic acid" & rt.ref == 2.38) #one of the few compounds with multiple annotations. One matches our data; removing the other for tracability purposes
  )

intersect = inner_join( expcpds, reference_list, by = c("name.exp" = "name.ref"), keep = TRUE) %>%
  select("name.exp", "name.ref", "mz", "rt.min", "rt.ref") %>%
  mutate(
    rt.diff = abs(rt.min - rt.ref),
    rt.match = ifelse( rt.diff < 3, TRUE, FALSE)
  ) %>%
  arrange(desc(rt.diff))
length(unique(intersect$name.exp))

# These are the CD metabolties that matched our reference list but the RT was > 3 minutes different. We will remove these! 
discard = intersect %>%
  filter( rt.match == FALSE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(discard$name.exp))

# this is a summary of the CD metabolite-RT combinations that DID match our reference list. We will keep these!
keep = intersect %>%
  filter( rt.match == TRUE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(keep$name.exp))


# anti_join() will remove the observations from 'discard' from 'exp_data'
expdata_2 = anti_join(expdata, discard) %>%
  left_join(., keep[,c("name.exp","mz","rt.min","rt.match")]) %>%
  rename(rt.ref.match = "rt.match") %>%
  select(!c("name.exp"))

n = expdata_2 %>%
  group_by(name, mz, rt.min, annot.source.predicted, annot.source.mzcloud, annot.source.masslist, mzcloud.n.hits, mzcloud.best.match, mzcloud.conf,rt.ref.match) %>%
  summarise(
    count = n()
  )
length(unique(n$name))

# -- mz cloud filtering: discard metabolites
mz_filt = n %>% filter( is.na(rt.ref.match), mzcloud.conf < 40 ) %>% select(!c("count"))
length(unique(mz_filt$name))

# anti_join() will remove the observations from 'mz_filt' from 'exp_data2'
expdata_3 = anti_join(expdata_2, mz_filt)
length(unique(expdata_3$name))

write.csv(expdata_3, file.path(processed_dir, "CD data filtered by mismatches to ref and mz cloud conf.csv"), row.names = FALSE)

# **** in the future this section could be a good place to indicate matches to our RT library


# ========== 4.0 - De-replicate repeated names ==========
# -- Manage repeated names with different RTs (but not yet consolidating by ion mode)

expdata_3 = read.csv(file.path(processed_dir, "CD data filtered by mismatches to ref and mz cloud conf.csv"))

# Apply a function select_top_isomer2()
# *** Among multiple name annotations, it determines which is the greatest based on median,
#     and then labels subsequent annotations as isomer 1..2..3 etc.

# make a list of data frames for every sequenced batch (Raw.File)
# recall Raw.File corresponds to Pos or Neg mode
batch_list = expdata_3 %>%
  group_by( raw.file ) %>%
  group_split(.)

batch_list2 = lapply(batch_list, select_top_isomer3)
# ... more splitting of major peaks in these resulting chromatograms

# data frame with secondary isomers annotated
expdata_4 = bind_rows(batch_list2)

# -- remove secondary isomers --
expdata_5= subset(expdata_4, grepl(": isomer", expdata_4$name) == FALSE)
length(unique(expdata_5$name))

# -- Consolidate the compounds with several RTs and MZs for consistent annotations
# .. at this point any duplicate annotations are within RT range of 20 seconds
# .. previously this occurred in select_top_isomer() step (5)
expdata_6 = expdata_5 %>%
  group_by(raw.file, sequence.id, name, formula) %>%
  arrange(rt.min, .by_group = TRUE) %>%   # ensures rows are ordered by rt.min
  summarise(
    rt.min = first(rt.min),               # smallest RT in group
    mz = round(first(mz),4),                       # corresponding MZ for that RT
    #ion.count = sum(ion.count, na.rm = TRUE)  # sum ICs if multiple annotations
    #** revise to pickk the largest annotation PER SAMPLE
    ion.count = max(ion.count, na.rm = TRUE)
  ) %>%
  ungroup()


# -- lastly, remove specific metabolites flagged by BC and SJ --
filt_metabs = read_metadata_csv(here("data/libraries/081225_list of unique metabolites_Bryan_SJ.csv")) %>%
  filter(flag == 1) # remove metabolites flagged with "1"

expdata_7 = expdata_6 %>%
  filter(!name %in% filt_metabs$name)
length(unique(expdata_7$name))

# -- check the number of occurences of each metabolite Name --
t = expdata_7 %>%
  mutate(sequence.id = sub("_", "@", sequence.id, fixed = TRUE)) %>%
  separate(sequence.id, into = c("ion.mode", "sample"), sep = "@") %>%
  group_by( name, ion.mode, formula, mz, rt.min ) %>%
  summarise(median.ic = median(ion.count))

t2 = t %>% group_by(name) %>%
  summarise(n=n()) %>%
  arrange(desc(n))# 1 or 2 for everything, depending on n ion modes for the metabolite
length(unique(t2$name))

# -- save dereplicated, with and without secondary isomers
write.csv(expdata_4, file.path(processed_dir, "annotated and dereplicated CD obs with repeat isomers.csv"), row.names = FALSE)
write.csv(expdata_7, file.path(processed_dir, "annotated and dereplicated CD obs without repeat isomers.csv"), row.names = FALSE)


# -- summarise filtered metabolites in a single data frame --

# April 15 2026 - low priority SKIP for now


# (A) RT matches to reference annotations
rtmatch = rbind(keep, discard)
df1 = expcpds %>%
  left_join(., rtmatch) %>%
  mutate(
    ref.match = case_when(
      rt.match == TRUE  ~ "Name + RT match to reference",
      rt.match == FALSE ~ "Name match, RT mismatch",
      is.na(rt.match)   ~ "Not in reference list"
    )
  ) %>%
  select(!c("rt.match"))

# (B) mz cloud filtering
df2 = df1 %>%
  mutate(
    mz.cloud.filt = case_when(
      mzcloud.conf < 40 ~ "Excluded, below threshold",
      mzcloud.conf >= 40 ~ "Passing, above threshold",
      is.na(mzcloud.conf) ~ "Passing, no score"
    )
  )

# (C) Isomers flagged from select_top_isomer
# *Recall if the annotations were within 20 seconds, we extracted the min RT and MZ for the group
# *This will make it so the topisomer results wont merge with df2 based on mz and rt.min
topisomer = expdata_4 %>%
  mutate(
    name.isom = name,
    name = gsub(" : isomer \\d+(\\.\\d+)?","",name),
    isomer = ifelse(grepl(": isomer", name.isom),"Secondary isomer","Primary annotation cluster"),
  ) %>%
  distinct(name, raw.file, formula, mz, rt.min, isomer)

df3 = df2 %>%
  left_join(., topisomer) %>%
  mutate(
    isomer = ifelse(ref.match == "Name match, RT mismatch" | mz.cloud.filt == "Excluded, below threshold", "Filtered in prior steps", isomer)
  )


# this came down to an error when writing out the file with a space  after the .csv extension ???
a = df2 %>% filter(name == "Caffeine") %>% select(name, raw.file, formula, mz, rt.min) %>% mutate(mz = as.numeric(mz))
b = topisomer %>%  filter(name == "Caffeine")%>% mutate(mz = as.numeric(mz))

join = a %>% left_join(.,b)

c = expdata %>% filter(name == "Caffeine") %>% distinct(name, raw.file, formula, mz, rt.min)
d = expdata_3 %>% filter(name == "Caffeine") %>% distinct(name, raw.file, formula, mz, rt.min)

a$mz[1] == b$mz[1] # F (expcpds > df1 > df2) != (expdata_4 > topisomer)
b$mz[1] == c$mz[1] # F (expdata_4 > topisomer) != expdata
a$mz[1] == c$mz[1] # T (expcpds > df1 > df2) == expdata
# therefore b - from expdata_4 -> topisomer is the problem

d$mz[1] == a$mz[1] # F expdata_3 != (expcpds > df1 > df2)
d$mz[1] == b$mz[1] # T expdata_3 == (expdata_4 > topisomer)
d$mz[1] == c$mz[1] # F expdata_3 != expdata 




# (E) Indicate the annotation RT that was adopted
min_annot = expdata_6 %>%
  distinct( name, rt.min, mz) %>%
  mutate(selected.rt = "selected annotation (rt+mz)")

df4 = df3 %>%
  left_join(., min_annot) %>%
  mutate(
    #selected.rt = ifelse( is.na(selected.rt), "consolidated into selected annotation", selected.rt),
    selected.rt = case_when(
      is.na(selected.rt) & (isomer != "Filtered in prior steps") & (mz.cloud.filt != "Excluded, below threshold") ~ "consolidated into selected annotation",
      is.na(selected.rt) & ((isomer == "Filtered in prior steps") | (mz.cloud.filt == "Excluded, below threshold")) ~ "Filtered in prior steps",
      TRUE ~ selected.rt
    )

  )

# (D) Flag to remove
df5 = df4 %>%
  mutate(
    manual.flag = ifelse(name %in% filt_metabs$name, "flagged for removal", "keep")
  )

# ** tidy and save the summary df
df6 = df5 %>%
  rename(rt="rt.min") %>%
  mutate(
    ion.mode = case_when(
      grepl("_Neg",raw.file) ~ "N",
      grepl("_Pos",raw.file) ~ "P",
    )
  ) %>%
  select(name, ion.mode, formula, mz, rt, ref.match, mz.cloud.filt, isomer, selected.rt, manual.flag)

# -- save processing summary
write.csv(df6, file.path(processed_dir, "summary of compound processing through step 4.csv"), row.names = FALSE)


# ========== 5.0 - Consolidate negative and positive modes ==========
# -- Select the ion mode with best quality annotation for metabolites with annotations in both neg and pos mode

# read in dereplicated data
exp_data = read.csv(file.path(processed_dir, "annotated and dereplicated CD obs without repeat isomers.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate( sequence.id = sub("_", "@", sequence.id, fixed = TRUE)) %>%
  separate( sequence.id, into = c("ion.mode", "sample"), sep = "@")
length(unique(exp_data$name))

# -- Determine the median of every metabolite 'Name' per mode
ic_medians = exp_data %>%
  group_by( name, ion.mode, formula, mz, rt.min ) %>%
  summarise(median.ic = median(ion.count))

# make sure that at most 1 mode is present
# ... should be max 2 per name if we did prior scripts correctly
n_modes = ic_medians %>% group_by(name) %>% summarise(n.modes = n()) 
# ... looks good

# -- Extract the annotation with the larger median across the samples
top_mode = ic_medians %>% group_by(name) %>% top_n(1, median.ic) %>% mutate( median.ic = NULL )
# 1500 rows; we have 1640 unique metabolite Names so we know we selected one each

# Select from exp_data only the top median annotations
exp_data2 = inner_join(exp_data, top_mode)

# check that we only have 1 mode per metabolite
t = exp_data2 %>%
  group_by( name, ion.mode, formula, mz, rt.min ) %>%
  summarise(median.ic = median(ion.count))

# should now only be 1 per name
n_modes2 = t %>% group_by(name) %>% summarise(n.modes = n()) 

exp_data3 = exp_data2 %>%
  mutate( sequence.id = paste0(ion.mode, "_", sample)) %>%
  select(raw.file, sequence.id, sample, ion.mode, everything())
length(unique(exp_data3$name))

# 
t2 = exp_data3 %>%
  distinct(name, mz) %>%
  group_by(name) %>%
  summarise(n())
# 9-decenoic acid has two formulas starting from step 2

n = exp_data3 %>%
  distinct(name, mz)

# -- save the data of top ion modes
write.csv(exp_data3, file.path(processed_dir, "annotated and processed CD obs consolidated by top ion mode.csv"), row.names = FALSE)

#

# ========== 6.0 - Merge Maven annotations with Compound Disc. Data ==========
# -- tidy/Add/replace the hand picked El-Maven observations

# (A) -- Tidy the maven data, and make data frames of unique Maven and CD compounds --

# Summer 2024 LCMS maven data were copied to this repo to "/data/raw/maven" from OneDrive "../LC-MS/Raw data/MAVEN"
mav_files = list.files(raw_maven_dir, pattern = "\\.csv$", full.names = TRUE)
csv_data_list = list()

for (csv_file in mav_files) {
  
  #csv_file = mav_files[12] # debug line
  csv_data = read.csv(csv_file)
  filename = sub(".*/", "", csv_file)
  
  csv_long = csv_data %>%
    # filter(
    #   !adductName == "" # in the LDL analysis this was necessary to avoid errors, but this time around it isn't needed.
    #        ) %>% 
    mutate(
      # Maven valine naming varies by cohort/export; normalize unlabeled/labeled valine consistently.
      compound = ifelse(compound %in% c("m37_Valine", "L-Valine"), "Valine", compound),
      compound = ifelse(compound == "Valine" & isotopeLabel == "C13N15-label-5-1", "Valine-C13N15", compound)
    ) %>%
    pivot_longer(., cols = grep("^[PN]_.+_[A-Za-z]+\\d+$", names(csv_data), value = TRUE), names_to = "Sequence.ID", values_to = "Ion.Count") %>%
    select("Sequence.ID", "compound", "formula", "medRt", "medMz","Ion.Count") %>%
    mutate(
      Sample = sub("^[PN]_", "", Sequence.ID),
      Ion.Mode = ifelse(grepl("^Pos", filename), "P", "N"),
      Raw.File = filename
    )
  
  csv_data_list[[csv_file]] = csv_long
}

mav_obs = bind_rows(csv_data_list) %>%
  mutate(
    Name = compound,
    Formula = formula
  ) %>%
  select( Raw.File, Sequence.ID, Ion.Mode, Sample, Name, Formula, medRt, medMz, Ion.Count) %>%
  rename_with(~ tolower(.)) %>% #lowercase colnames
  mutate(
    # Nov 19 2025 - rename L-Histidine to Histidine to be consistent with existing data, and match RT and MZ
    name = ifelse(name == "L-Histidine", "Histidine", name),
    name = ifelse(name == "Glycine (1)","Glycine",name)
  ) %>%
  group_by(raw.file, sequence.id, ion.mode, sample, name, formula, medrt, medmz)
write.csv(mav_obs, file.path(processed_maven_dir, "all maven observations.csv"), row.names = FALSE)


# -- Make df of maven compounds --
mav_compounds = mav_obs %>%
  group_by(ion.mode, name, formula) %>%
  summarise(
    mz = median(medmz),
    med.ic = median(ion.count),
    rt = median(medrt),
    n.samples = n()
  ) %>%
  mutate( source = "maven" ) %>%
  select(
    source, name, ion.mode, formula, mz, rt, med.ic, n.samples
  )

write.csv(mav_compounds, file.path(processed_maven_dir, "unique compounds and formulas from maven.csv"), row.names = FALSE)


# -- Create a data frame of the unique compounds found in the CD data --
cd_compounds = read.csv(file.path(processed_dir, "annotated and processed CD obs consolidated by top ion mode.csv")) %>%
  group_by(ion.mode, name, formula, rt.min) %>%
  summarise(
    mz = median(mz),
    med.ic = median(ion.count),
    n.samples = n()
  ) %>%
  mutate( source = "CD" ) %>%
  rename(rt = rt.min) %>%
  select(
    source, name, ion.mode, formula, mz, rt, med.ic, n.samples
  )

n = cd_compounds %>% distinct(name , mz)

write.csv(cd_compounds, file.path(processed_maven_dir, "unique compounds and formulas from compound discoverer.csv"), row.names = FALSE)

names(cd_compounds)
names(mav_compounds)



# -- inspect name overlaps between CD and Maven before defining duplicate annotations.
all_compounds = rbind(cd_compounds, mav_compounds) %>%
  arrange(name, formula, source, ion.mode)

shared_formulas = intersect(cd_compounds$formula, mav_compounds$formula)

formula_audit = all_compounds %>%
  filter(formula %in% shared_formulas) %>%
  arrange(mz, rt)

write.csv( formula_audit, file.path(processed_maven_dir, "formula inspection bw cd and mav.csv"), row.names = FALSE)
# ... maually flag keep/remove in excel
formula_audit2 = read.csv(file.path(processed_maven_dir, "formula inspection bw cd and mav_flagged.csv"))




# *** Having completed this manual compariosn, we can move to the next step of merging the annotations

# -- load the CD and Maven complete data and combine --
cd_data = read.csv(file.path(processed_dir, "annotated and processed CD obs consolidated by top ion mode.csv"), fileEncoding = "UTF-8-BOM") %>%
  rename( rt = rt.min) %>%
  mutate(
    source = "CD"
  )

mav_data = read.csv(file.path(processed_maven_dir, "all maven observations.csv"), fileEncoding = "UTF-8-BOM") %>%
  rename(
    rt = medrt,
    mz = medmz
  ) %>%
  mutate(
    source = "maven"
  ) %>%
  group_by(name, formula) %>%
  mutate(
    rt = median(rt),
    mz = median(mz)
  )

expdata = rbind(cd_data, mav_data) 

# -- Load the manual formula audit decisions.
# ... remove == 1 marks the specific source/name/mode/formula annotation to discard.

action = formula_audit2 %>%
  filter(remove == 1) %>%
  distinct(source, name, ion.mode, formula) %>%
  mutate(flag = "remove")

# -- merge with expdata --
expdata2 = expdata %>%
  left_join(action, by = c("source", "name", "ion.mode", "formula")) %>%
  filter(
    is.na(flag)
  ) %>%
  select(
    c("raw.file", "sequence.id", "sample", "ion.mode", "source", "name", "formula", "mz", "rt", "ion.count")
  )

t = expdata2 %>% distinct(name, formula, mz, rt, ion.mode) %>% group_by(name) %>% summarise(n()) %>% arrange(desc(`n()`))
t = expdata2 %>% distinct(name, formula, mz, rt, ion.mode) %>% group_by(name) %>% mutate(n=n()) %>% filter(n>1) %>% arrange(desc(n))
# just labeled valine remains as dupes which is what we want for normalization later :)

# -- save the CD + Mav data --
write.csv(expdata2, file.path(processed_dir, "annotated and processed CD plus MAV obs.csv"), row.names = FALSE)




# ========== 6.0 - Valine and Median Normalization ==========
# -- Valine and Median normalization of the processed data

expdata = read.csv(file.path(processed_dir, "annotated and processed CD plus MAV obs.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate(
    cohort = "NA",
    tissue = sub("_.*", "", sample),
    animal = str_extract(sample, "(?<=_).*")
    )

# check that there arent redundant names with different MZs
t = expdata %>%
  ungroup() %>%
  distinct(name, mz) %>%
  group_by(name) %>%
  summarise( n() )
# all 1 obs looks good!

# examine distributions by certain vars
ggplot( subset(expdata, source == "CD"), aes(x= log10(ion.count))) +
  geom_histogram(bins = 1000, aes(fill = factor(ion.mode))) +
  facet_wrap(~ion.mode) +
  ggtitle("Ion count distribution by ion mode")


# -- Valine Normalization

# normalize negative ion.mode metabs to labeled valine in neg; and pos with pos
valine = expdata %>%
  filter(name == "Valine-C13N15") %>%
  group_by(sequence.id) %>%
  summarise(
    val.ic = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

expdata2 = left_join(expdata, valine)

nrow(subset(expdata2, is.na(val.ic))) # 0, good!

expdata3 = expdata2 %>%
  group_by(ion.mode) %>% # 2026-01-28 Ion mode (w/o tissue) is correct grouping
  mutate(
    med.val = median(val.ic, na.rm = TRUE),
    val.adj.ic = (ion.count / val.ic) * med.val
  ) %>%
  ungroup()


# -- 3.0 - median normalization
expdata4 = expdata3 %>%
  group_by(tissue, name) %>% 
  mutate(
    med.metab = median(val.adj.ic, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    `log10(ic/med.metab)` = log10(val.adj.ic/med.metab)
  )%>%
  group_by(sample) %>%
  # is it necessary to group the norm also by ion mode here? I would argue No since the distributions of pos and neg are in the same range
  mutate(
    med.sample = median(`log10(ic/med.metab)`, na.rm = TRUE),
    med.adj.ic = ( 10^(`log10(ic/med.metab)` - med.sample)*med.metab )
  ) %>%
  select(
    "raw.file", "sequence.id", "ion.mode", "sample", "tissue", "cohort", "animal",
    "name", "formula", "source", 
    "mz", "rt", "ion.count", "val.adj.ic", "med.adj.ic"             
  ) %>%
  ungroup()


# -- 4.0 - lastly rename, remove, and add classes to a few things based on BC list and write out dfs

remove = c(
  "Ketorolac","(Â±)-Metalaxyl","5-Sulfosalicylic acid","Acetaminophen glucuronide","Acipimox","Phenol",
  "Aflatoxin B1","Barbital","Bisoprolol","Carbamazepine-O-quinone","Ceftiofur","Citicoline","Dedimethylchlorpromazine",
  "Dehydronorketamine", "Donepezil metabolite M4","Ephedrine","Eplerenone","Ketamine","Lomefloxacin",
  "Metharbital","Midodrine","Myristicin","Netilmicin", "Norketamine", "Oxybuprocaine", "Phenmetrazine", "Pyridostigmine",
  "Salbutamol","Terbutaline","clonazolam","laurolactam","Pramanicin isomer","trimethoprim impurity b isomer"
)
expdata5 = expdata4 %>%
  filter(
    !(name %in% remove) | grepl("4,13,14-trihydroxy-9-oxo-8,17-dioxatetracyclo", name)
  )


bc1 = read_metadata_csv(here("data/libraries/081825_list of unique metabolites_Bryan.csv")) %>%
  select(name, formula, new.name, delete, metabolite.classification) %>%
  mutate(
    new.name = ifelse(new.name == "",NA,new.name),
    
    #fix methylbutanol annotation rename
    name = ifelse(name == "(횂짹)-2-Methylbutanal", "(Â±)-2-Methylbutanal", name),
    new.name = ifelse(name == "(Â±)-2-Methylbutanal","2-Methylbutanal", new.name),
    
    # myoinositol make sure the class merges later
    name = ifelse(name == " myo-inositol", "Myo-inositol", name),
    name = ifelse(name == "Cholic Acid", "Cholic acid", name)
    
  ) %>%
  rename(
    class = "metabolite.classification"
  ) %>%
  rename(
    new.name1 = new.name,
    class1 = class,
    delete1 = delete
  )

# *** 2/16/26 - incorperate new updates to a few of the fields ***
# -- if a value is present in new name replace name with new name
# -- new class is the new classification field (will do in next step in winde format)
#bc2 =  read.csv(here("data/libraries/tidy_and_normalized_ion_count_data_na_-_wide_bryan_v3.csv")) %>%
bc2 = read_metadata_csv(first_existing_path(c(
  here("data/libraries/tidy and normalized ion count data AA - wide_Bryan_v3.csv"),
  here("data/libraries/tidy_and_normalized_ion_count_data_na_-_wide_bryan_v3.csv")
))) %>%
  select(name, formula, new.name, delete, new.class) %>%
  rename(
    new.name2 = new.name,
    class2 = new.class
  ) %>%
  mutate(
    new.name2 = ifelse( new.name2 == "" , NA, new.name2),
    # this is necessary to get the merge right
    name = ifelse(name == "(횂짹)-2-Hydroxy-4-(methylthio)butanoic acid", "(Â±)-2-Hydroxy-4-(methylthio)butanoic acid", name),
    name = ifelse(name == "Acetyl-棺-methylcholine", "Acetyl-β-methylcholine", name),
    name = ifelse(name == "2-Nitro-5-(phenylsulfonyl)phenol","2-nitro-5-(phenylsulfonyl)phenol", name)
    
  ) %>%
  rename(
    delete2 = delete
  )

expdata6 = expdata5 %>%
  # apply changes from bc1 (Bryan's first manual check of the wide data frame)
  left_join(., bc1) %>%
  mutate(
    name = coalesce(new.name1, name), # return name in first non NA col in the indicated order
    # # rename a few items by request
    name = ifelse(name == "N-Palmitoyl Glutamic acid", "N-Palmitoyl glutamate", name),
    name = ifelse(name == "3-keto Fusidic acid", "3-Keto fusidic acid", name),
    name = ifelse(name == "butyl isodecyl phthalate  isomer", "Butyl isodecyl phthalate isomer", name),
    name = ifelse(name == "phenylglycol 3-o-sulfate isomer", "Phenylglycol 3-o-sulfate isomer", name),
    name = ifelse(name == "aminoisophthalic acid isomer", "Aminoisophthalic acid isomer", name),
    name = ifelse(name == "phenylglycol 3-o-sulfate isomer", "Phenylglycol 3-o-sulfate isomer", name),
    name = ifelse(name == "2-aminophenol", "2-Aminophenol", name),
    name = ifelse(name == "luguine isomer", "Luguine isomer", name),
    name = ifelse(name == "De-O-methylsterigmatocystin", "De-O-Methylsterigmatocystin", name),
    
    # add these classes by request
    class1 = ifelse(name == "D-Fructose", "Carbohydrates and carbohydrate conjugates", class1),
    class1 = ifelse(name == "L-Allothreonine", "Amino acids and derivatives", class1),
    class1 = ifelse(name == "N-Palmitoyl glutamate", "Amino acids and derivatives", class1),
    class1 = ifelse(name == "Taurodeoxycholic acid", "Bile acids", class1),
    class1 = ifelse(name == "Cresol sulfate isomer", "Microbiota metabolites", class1)
    
  ) %>%
  filter(
    # requested to remove these
    is.na(delete1),
    !name %in% c("Gentisic acid","Cuminaldehyde"),
    !(name == "Glycoursodeoxycholic acid" & source == "CD"),
    !(name == "Cetrimonium") # duplicate entry metabolite
  ) %>%
  select(
    !c("new.name1","delete1")
  ) %>%
  # apply changes from bc2 (Bryan's second manual check of the wide data frame)
  left_join(., bc2) %>%
  # Cetrimonium is at ttwo RTs
  mutate(
    name = coalesce(new.name2, name),
    class = coalesce(class2, class1)
  ) %>%
  filter(
    is.na(delete2)
  ) %>%
  select(!c("class1","new.name2","class2","delete2"))


replacement_names <- c(
  "N-Palmitoyl glutamate",
  "3-Keto fusidic acid",
  "Butyl isodecyl phthalate isomer",
  "Phenylglycol 3-o-sulfate isomer",
  "Aminoisophthalic acid isomer",
  "Phenylglycol 3-o-sulfate isomer",
  "2-Aminophenol",
  "Luguine isomer"
)
replacement_names %in% expdata6$name

write.csv(expdata6, file.path(processed_dir, paste0("tidy and normalized ion count data - TISSUE NA.csv")), row.names = FALSE)


# -- wide format
expwide = expdata6 %>%
  group_by(sample, name, formula, ion.mode, source, mz, rt, class) %>%
  summarise(
    med.adj.ic = max(med.adj.ic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(
    sample, name, formula, ion.mode, source, mz, rt, class, med.adj.ic
  ) %>%
  pivot_wider(names_from = "sample", values_from = "med.adj.ic") %>%
  select(name, formula, ion.mode, source, mz, rt, class, everything())

write.csv(expwide, file.path(processed_dir, paste0("tidy and normalized ion count data TISSUE NA - wide.csv")))
