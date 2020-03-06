
# Import
library(VennDiagram)

###################################################################################################
# FUNCTIONS
###################################################################################################

# Input
#in2FileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/1_Tables/table_pair.txt"
in2FileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/1_Tables/table_pair_rep.txt"
out2Loc = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/2_Venn2/"
table2 = read.table(in2FileName, header = FALSE, sep = "\t")
in3FileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/1_Tables/table_triple.txt"
out3Loc = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/2_Venn3/"
table3 = read.table(in3FileName, header = FALSE, sep = "\t")

###################################################################################################
# FUNCTIONS
###################################################################################################

# Pairwise Venn Function
venn2Diagram <- function(area1, area2, inters, label1, label2, outFileName){
  label1 = gsub("_", "\n", label1)
  label2 = gsub("_", "\n", label2)
  pdf(file = outFileName, width = 8, height = 8)
  par(mar = c(1,1,1,1))
  draw.pairwise.venn(area1, area2, inters, category = c(label1, label2), lty = rep("blank", 2), fill = c("gold3", "deepskyblue"), alpha = rep(0.5, 2), cex = 1.5,
                     cat.pos = c(-50, 50), cat.dist = rep(0.05, 2), cat.cex = rep(2.0, 2), scaled = TRUE, print.mode = c("raw", "percent"), sigdigs = 3, margin = c(0.1,0.1,0.1,0.1))
  dev.off()
}

# Triple Venn Function
venn3Diagram <- function(area1, area2, area3, n12, n23, n13, n123, label1, label2, label3, outFileName){
  pdf(file = outFileName, width = 8, height = 8)
  par(mar = c(1,1,1,1))
  draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c(label1, label2, label3), lty = rep("blank", 3), fill = c("gold3", "deepskyblue", "darksalmon"), alpha = rep(0.5, 3), cex = 1.5,
                     cat.pos = c(-40, 40, 180), cat.dist = rep(0.05, 3), cat.cex = rep(2.0, 3), scaled = TRUE, print.mode = c("raw", "percent"), sigdigs = 3, margin = c(0.1,0.1,0.1,0.1))
  dev.off()
}

###################################################################################################
# EXECUTION
###################################################################################################

# Execution

for(i in 1:nrow(table2)){

  name1 = as.character(table2[i,1])
  name2 = as.character(table2[i,2])
  area1 = table2[i,3]
  area2 = table2[i,4]
  inters = table2[i,5]
  outFileName = paste(out2Loc, name1, "_", name2, ".pdf", sep = "")

  venn2Diagram(area1, area2, inters, name1, name2, outFileName)

}

for(i in 1:nrow(table3)){

  name1 = as.character(table3[i,1])
  name2 = as.character(table3[i,2])
  name3 = as.character(table3[i,3])
  area1 = table3[i,4]
  area2 = table3[i,5]
  area3 = table3[i,6]
  n12 = table3[i,7]
  n23 = table3[i,8]
  n13 = table3[i,9]
  n123 = table3[i,10]
  outFileName = paste(out3Loc, name1, "_", name2, "_", name3, ".pdf", sep = "")

  venn3Diagram(area1, area2, area3, n12, n23, n13, n123, name1, name2, name3, outFileName)

}



#draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category =
#rep("", 3), rotation = 1, reverse = FALSE, euler.d =
#TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
#rep("solid", 3), col = rep("black", 3), fill = NULL,
#alpha = rep(0.5, 3), label.col = rep("black", 7), cex
#= rep(1, 7), fontface = rep("plain", 7), fontfamily =
#rep("serif", 7), cat.pos = c(-40, 40, 180), cat.dist =
#c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
#cat.cex = rep(1, 3), cat.fontface = rep("plain", 3),
#cat.fontfamily = rep("serif", 3), cat.just =
#list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
#= "outer", cat.prompts = FALSE, rotation.degree = 0,
#rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
#0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
#sigdigs = 3, direct.area = FALSE, area.vector = 0,
#...)


#draw.pairwise.venn(area1, area2, cross.area, category = rep("", 2),
#euler.d = TRUE, scaled = TRUE, inverted = FALSE,
#ext.text = TRUE, ext.percent = rep(0.05, 3), lwd =
#rep(2, 2), lty = rep("solid", 2), col = rep("black",
#2), fill = NULL, alpha = rep(0.5, 2), label.col =
#rep("black", 3), cex = rep(1, 3), fontface =
#rep("plain", 3), fontfamily = rep("serif", 3), cat.pos
#= c(-50, 50), cat.dist = rep(0.025, 2), cat.cex =
#rep(1, 2), cat.col = rep("black", 2), cat.fontface =
#rep("plain", 2), cat.fontfamily = rep("serif", 2),
#cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos
#= "outer", cat.prompts = FALSE, ext.pos = rep(0, 2),
#ext.dist = rep(0, 2), ext.line.lty = "solid",
#ext.length = rep(0.95, 2), ext.line.lwd = 1,
#rotation.degree = 0, rotation.centre = c(0.5, 0.5),
#ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop =
#NULL, print.mode = "raw", sigdigs = 3, ...)
