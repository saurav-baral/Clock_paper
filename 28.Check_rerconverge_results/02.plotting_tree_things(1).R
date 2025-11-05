
# 
setwd("C:/Users/sauba/Desktop/20.RERConverge_trimmed/3.New_colias/61.New_alignment_colias/28.Check_rerconverge_results")

# for this I am using phytools
# http://blog.phytools.org/2024/03/function-for-plotting-discrete-andor.html

library(phytools)

library(ape); library(phytools); library(TreeTools)
packageVersion("phytools")


tree<-read.tree("tree_tree.nwk")
# x<-as.matrix(read.csv("x.csv",row.names=1))[,1]

if (is.null(tree$edge.length)) tree$edge.length <- rep(0.001, nrow(tree$edge))

tree <- compute.brlen(tree, method="Grafen")

pdf("State_tree_color.pdf", height = 20, width = 10)
state.data<-read.csv("two_state_species.csv",row.names=1)
fmode<-as.factor(setNames(state.data[,1],rownames(state.data)))
dotTree(tree,fmode,colors=setNames(c("blue","red"),
                                       c("diapause","no_diapause")),ftype="i",fsize=1)
dev.off()





library(ape)
library(phytools)

# --- Load tree ---
tree <- read.tree("tree_tree.nwk")
if (is.null(tree$edge.length)) tree$edge.length <- rep(0.001, nrow(tree$edge))
tree <- compute.brlen(tree, method = "Grafen")

# --- Load binary state data ---
state.data <- read.csv("two_state_species.csv", row.names = 1)
fmode <- as.factor(setNames(state.data[,1], rownames(state.data)))

# --- Load continuous data ---
cont.data <- read.csv("period_rers.csv", row.names = 1)
cont <- setNames(cont.data[,1], rownames(cont.data))
cont <- cont[tree$tip.label]  # align order with tree tips

# --- Handle missing values ---
valid <- !is.na(cont)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
cont.col <- rep("grey", length(cont))  # default for missing
cont.col[valid] <- colors[as.numeric(cut(cont[valid], breaks = 100))]

# --- Plot ---
pdf("State_tree_with_continuous_dot.pdf", height = 20, width = 10)

dotTree(
  tree,
  fmode,
  colors = setNames(c("blue", "red"), c("diapause", "no_diapause")),
  ftype = "i", fsize = 1
)

# --- Add continuous-value dots ---
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_x <- max(obj$xx) * 1.02
points(rep(tip_x, length(tree$tip.label)), obj$yy[1:length(tree$tip.label)],
       col = cont.col, pch = 19, cex = 1)

# --- Add continuous legend (continuous colorbar) ---
par(xpd = TRUE)
legend_x <- tip_x + (max(obj$xx) - min(obj$xx)) * 0.05
legend_ybottom <- min(obj$yy)
legend_ytop <- max(obj$yy)
color.legend(
  xleft = legend_x,
  ybottom = legend_ybottom,
  xright = legend_x + 0.2,
  ytop = legend_ytop,
  legend = round(seq(min(cont, na.rm = TRUE), max(cont, na.rm = TRUE), length.out = 5), 2),
  rect.col = colors,
  gradient = "y",
  align = "right",
  cex = 1,
  title = "Continuous trait"
)

dev.off()






library(ape)
library(phytools)

# --- Load tree ---
tree <- read.tree("tree_tree.nwk")
if (is.null(tree$edge.length)) tree$edge.length <- rep(0.001, nrow(tree$edge))
tree <- compute.brlen(tree, method = "Grafen")

# --- Load binary state data ---
state.data <- read.csv("two_state_species.csv", row.names = 1)
fmode <- as.factor(setNames(state.data[,1], rownames(state.data)))

# --- Load continuous data ---
cont.data <- read.csv("period_rers.csv", row.names = 1)
cont <- setNames(cont.data[,1], rownames(cont.data))
cont <- cont[tree$tip.label]  # align order with tree tips

# --- Handle missing values ---
valid <- !is.na(cont)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
cont.col <- rep("grey", length(cont))  # grey for missing
cont.col[valid] <- colors[as.numeric(cut(cont[valid], breaks = 100))]

# --- Plot ---
pdf("State_tree_with_continuous_dot.pdf", height = 20, width = 10)

dotTree(
  tree,
  fmode,
  colors = setNames(c("blue", "red"), c("diapause", "no_diapause")),
  ftype = "i", fsize = 1
)

# --- Add continuous dots to the right of labels ---
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_x <- max(obj$xx) * 1.05
points(rep(tip_x, length(tree$tip.label)), obj$yy[1:length(tree$tip.label)],
       col = cont.col, pch = 19, cex = 1)

# --- Draw continuous colorbar legend ---
legend_x <- tip_x + 0.02
legend_width <- 0.03
legend_height <- 5
ncolors <- length(colors)

for(i in 1:ncolors){
  ybottom <- min(obj$yy) + (i-1)/ncolors * legend_height
  ytop <- min(obj$yy) + i/ncolors * legend_height
  rect(legend_x, ybottom, legend_x + legend_width, ytop, col=colors[i], border=NA)
}

# Add numeric labels to legend
legend_labels <- round(seq(min(cont, na.rm=TRUE), max(cont, na.rm=TRUE), length.out=5), 2)
legend_y <- seq(min(obj$yy), max(obj$yy), length.out=5)
text(x = legend_x + legend_width*1.01, y = legend_y, labels = legend_labels, adj = 0)

# Legend title
text(x = legend_x + legend_width*1.5, y = max(obj$yy) + 0.05*legend_height,
     labels = "Continuous trait", adj = 0, font = 2)

dev.off()


library(ape)
library(phytools)

# --- Load tree ---
tree <- read.tree("tree_tree.nwk")
if (is.null(tree$edge.length)) tree$edge.length <- rep(0.001, nrow(tree$edge))
tree <- compute.brlen(tree, method = "Grafen")

# --- Load binary state data ---
state.data <- read.csv("two_state_species.csv", row.names = 1)
fmode <- as.factor(setNames(state.data[,1], rownames(state.data)))

# --- Load continuous data ---
cont.data <- read.csv("period_rers.csv", row.names = 1)
cont <- setNames(cont.data[,1], rownames(cont.data))
cont <- cont[tree$tip.label]

# --- Handle missing values ---
valid <- !is.na(cont)
colors <- colorRampPalette(c("blue","white", "red"))(100)
cont.col <- rep("grey", length(cont))
cont.col[valid] <- colors[as.numeric(cut(cont[valid], breaks = 100))]

# --- Plot ---
pdf("State_tree_with_continuous_dot.pdf", height = 20, width = 10)

dotTree(
  tree,
  fmode,
  colors = setNames(c("red", "blue"), c("diapause", "no_diapause")),
  ftype = "i", fsize = 1
)

# --- Add continuous-value dots ---
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_x <- max(obj$xx) * 1.01
points(rep(tip_x, length(tree$tip.label)), obj$yy[1:length(tree$tip.label)],
       col = cont.col, pch = 19, cex = 1)

# --- Draw compact colorbar legend at top-right ---
# fixed coordinates relative to plot
usr <- par("usr")  # get figure coordinates: c(x1, x2, y1, y2)
legend_xleft <- usr[2] - 0.1 * (usr[2]-usr[1])
legend_xright <- usr[2] - 0.02 * (usr[2]-usr[1])
legend_ybottom <- usr[4] - 0.15 * (usr[4]-usr[3])
legend_ytop <- usr[4] - 0.02 * (usr[4]-usr[3])

ncolors <- length(colors)
for(i in 1:ncolors){
  yb <- legend_ybottom + (i-1)/ncolors * (legend_ytop - legend_ybottom)
  yt <- legend_ybottom + i/ncolors * (legend_ytop - legend_ybottom)
  rect(legend_xleft, yb, legend_xright, yt, col = colors[i], border = NA)
}

# Add numeric labels
legend_labels <- round(seq(min(cont, na.rm=TRUE), max(cont, na.rm=TRUE), length.out=5), 2)
legend_y <- seq(legend_ybottom, legend_ytop, length.out=5)
text(x = legend_xright + 0.01*(usr[2]-usr[1]), y = legend_y,
     labels = legend_labels, adj = 0)

# Legend title
text(x = legend_xleft, y = legend_ytop + 0.01*(usr[4]-usr[3]),
     labels = "Continuous trait", font=2, adj=0)

dev.off()



x
