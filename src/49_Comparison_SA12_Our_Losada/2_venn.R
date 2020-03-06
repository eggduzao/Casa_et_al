
# Import
library(VennDiagram)

# Function
vennDiagram <- function(area1, area2, inters, label1, label2, outFileName){
  label1 = gsub("_", "\n", label1)
  label2 = gsub("_", "\n", label2)
  pdf(file = outFileName, width = 8, height = 8)
  par(mar = c(1,1,1,1))
  draw.pairwise.venn(area1, area2, inters, category = c(label1, label2), lty = rep("blank", 2), fill = c("gold3", "deepskyblue"), alpha = rep(0.5, 2), cex = 1.5,
                     cat.pos = c(-50, 50), cat.dist = rep(0.025, 2), cat.cex = rep(2.0, 2), scaled = TRUE, print.mode = c("raw", "percent"), sigdigs = 3, margin = c(0.1,0.1,0.1,0.1))
  dev.off()
}

# Input
inFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/20_Comparison_SA12_Our_Losada/1_Table/table.txt"
outLoc = "/usr/users/egadegu/Projects/Wendt_Stag/Results/20_Comparison_SA12_Our_Losada/2_Venn/"
table = read.table(inFileName, header = FALSE, sep = "\t")

# Execution
for(i in 1:nrow(table)){

  area1 = table[i,3] + table[i,4]
  inters = table[i,4]
  area2 = table[i,5] + table[i,4]
  name1 = as.character(table[i,1])
  name2 = as.character(table[i,2])
  outFileName = paste(outLoc, name1, "_", name2, ".pdf", sep = "")

  vennDiagram(area1, area2, inters, name1, name2, outFileName)

}


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
