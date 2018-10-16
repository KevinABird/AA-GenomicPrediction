library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(gridExtra)

files <- dir(path="~/Documents/Bird_et_al/sigma_out/", pattern="*.txt", full.names=TRUE)
files

# convert files to .csv
for (i in 1:length(files)) {
  file=read.table(file=files[i], header=T)
  write.csv(file, file=paste0(sub(".txt","", files[i]), ".csv"))
}

data_path <- "~/Documents/Bird_et_al/results/sigma_out/"
files <- dir(path=data_path, pattern="*.csv")
files

# combine output into single data frame
sigmas <- data_frame(filename=files) %>%
  mutate(file_contents=map(filename, ~ read_csv(file.path(data_path, .))))

sigmas <- unnest(sigmas)

colnames(sigmas) <- c("AA", "n", "GS", "AAS", "C", "e")

sigmas$AA <- gsub(pattern="sigmas_|.csv", "", sigmas$AA)

sigmas$AA <- as.factor(sigmas$AA)
sigmas$AA <- recode(sigmas$AA, "ala"="Ala", "arg"="Arg", "asp"="Asp", "gln"="Gln", "glu"="Glu",
                                       "his"="His", "ile"="Ile", "leu"="Leu", "lys"="Lys", "met"="Met",
                                       "phe"="Phe", "pro"="Pro", "ser"="Ser", "thr"="Thr", "trp"="Trp",
                                       "try"="Tyr", "val"="Val", "Total"="TFAA")

sigmas <- subset(sigmas, ! AA %in% c("AspFam_Asp", "AspFam", "BCAA", "GluFam", "GluFam_glu", 
                                        "PyrFam", "SerFam", "ShikFam"))
sigmas <- droplevels(sigmas)
sigmas$AASprop <- sigmas$AAS/(sigmas$GS+sigmas$AAS+sigmas$C)
sigmas$CSprop <- sigmas$C/(sigmas$GS+sigmas$AAS+sigmas$C)
sigmas$ratio <- sigmas$AASprop - sigmas$CSprop

sigmas$AA <- reorder(sigmas$AA, reorder(gp_results$AA, -gp_results$Diff))

sigmas$ratio[sigmas$ratio==1] <- NA # set outliers to NA

gp_results <- read.table("~/Documents/Bird_et_al/AAS_GSdf.txt", header=TRUE)

gp_results$n <- rep(1:1000)

merged <- merge(sigmas, gp_results, by=c("AA", "n"))

merged$AA <- factor(merged$AA, levels=c("Ala", "Arg", "Asp", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                       "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "TFAA"))

# define colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myColors <- colorRampPalette(cbPalette)(19)
names(myColors) <- levels(merged$AA)
colScale <- scale_fill_manual(name="", values=myColors)




## create legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

leg <- ggplot(merged, aes(AA, ratio, fill=AA)) + geom_boxplot() + colScale + theme_bw()
leg
my_legend <- g_legend(leg)

p1 <- ggplot(merged, aes(reorder(AA, -Diff), Diff, fill=AA)) +
  geom_hline(yintercept=0, colour="gray", size=1.5) +
  geom_boxplot(outlier.size=0.5)

p1 <- p1 + colScale + theme_bw() +
  labs(x="Amino Acid", y=expression(r[AAS]-r[CS])) + 
  theme(legend.position="none") 

p1 <- ggdraw(p1) + draw_plot_label("A")

p2 <- ggplot(merged, aes(reorder(AA, -Diff), ratio, fill=AA)) +
  geom_hline(yintercept=0, colour="gray", size=1.5) +
  geom_boxplot(outlier.size=0.5) 

p2 <- p2 + colScale + theme_bw() +
  labs(x="Amino Acid", y=expression(paste(h[AAS]^{2}-h[CS]^{2}))) + theme(legend.position="none")

p2 <- ggdraw(p2) + draw_plot_label("B")

# setEPS()
png("~/Documents/Bird_et_al/figures/Fig2.png", width=7.5, height=7.5, units="in", res=300)
grid.arrange(p1, p2, my_legend, ncol=2, nrow=2, layout_matrix = rbind(c(1,3), c(2,3)), widths=c(4,0.5))
dev.off()
