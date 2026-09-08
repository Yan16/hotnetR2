# read TF list
library(openxlsx)

library(dplyr)
library(readr)
tf <- read.xlsx("data-raw/TF/TF-1600-mb1.xlsx", sheet = "Sheet2")
# readr::write_rds("TF1600.rds")

dim(tf)
# [1] 1639   32


colnames(tf)
#  [1] "gene.ID"
#  [2] "gene.Name"
#  [3] "gene.DBD"
#  [4] "Is.TF?"
#  [5] "f.TF.assessment"
#  [6] "f.Binding.mode"
#  [7] "f.Motif.status"
#  [8] "f.Notes"
#  [9] "f.Comments"
# [10] "n.Committee.notes"
# [11] "n.MTW.Notes"
# [12] "n.TRH.Notes"
# [13] "n.SL.notes"
# [14] "n.AJ.notes"
# [15] "cu.Disagree.on.Assessment"
# [16] "cu.Disagree.on.Binding"
# [17] "c1.Author1"
# [18] "c1.Assesment1"
# [19] "c1.Binding1"
# [20] "c1.Comment1"
# [21] "c1.Notes1"
# [22] "c2.Author2"
# [23] "c2.Assesment2"
# [24] "c2.Binding2"
# [25] "c2.Comment2"
# [26] "c2.Notes2"
# [27] "pc.Vaquerizas.2009.TF.classification"
# [28] "pc.CisBP.considers.it.as.a.TF?"
# [29] "pc.TFclass.considers.it.as.a.TF?"
# [30] "pc.TF-CAT.classification"
# [31] "pc.Is.a.GO.TF"
# [32] "pc.PDB"

table(tf$`Is.TF?`)
head(tf[, 1:3])
#           gene.ID gene.Name    gene.DBD
# 1 ENSG00000137203    TFAP2A        AP-2
# 2 ENSG00000008196    TFAP2B        AP-2
# 3 ENSG00000087510    TFAP2C        AP-2
# 4 ENSG00000008197    TFAP2D        AP-2
# 5 ENSG00000116819    TFAP2E        AP-2
# 6 ENSG00000116017    ARID3A ARID/BRIGHT

# save gene symbol
tf |>
    select(gene.Name) |>
    write_tsv("TF1600_symbol.tsv", col_names = FALSE)
# save ENSG
tf |>
    select(gene.ID) |>
    write_tsv("TF1600_ENSG.tsv", col_names = FALSE)
