library(dplyr)

raw_amyloid_db <- read.csv("py_con_zero/amyloid_explorer_2026_10_28.csv")

# Extract standard PDB IDs and trim unusual entries
pdb_ids <- raw_amyloid_db %>%
  mutate(
    PDB.ID = substr(PDB.ID, 1, 4)
  ) %>%
  filter(!is.na(PDB.ID), nchar(PDB.ID) == 4) %>%
  distinct(PDB.ID) %>%
  pull(PDB.ID)

writeLines(
  paste(pdb_ids, collapse = ","),
  "py_con_zero/pdb_ids.txt"
)

unusual_pdb_entries <- raw_amyloid_db[nchar(raw_amyloid_db$PDB.ID) > 4, ]
unusual_pdb_ids <- unique(unusual_pdb_entries$PDB.ID) # Please be aware that these files contain two different amyloid sequences
writeLines(
  paste(unusual_pdb_ids, collapse = ","),
  "py_con_zero/to_be_manually_edited_pdb_ids.txt"
)
