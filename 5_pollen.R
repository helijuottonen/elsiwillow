# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
# https://github.com/helijuottonen/elsiwillow

# 5: pollen remaining by site  / honeybees (Fig. 4)

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

library(ggplot2)
library(cowplot)

# read in pollen data
pol <- read.csv2("data/pollen_fig4.csv", header=T)

# rename sites to match manuscript names (I checked that the outcome is correct)
pol$Site <- as.factor(pol$Site)
levels(pol$Site)
levels(pol$Site) <- c("A", "B", "C", "D", "E", "F")

# plotting
pol_plot <- ggplot(pol, aes(x=reorder(Site, trt_ratio_casy, FUN=median), y=trt_ratio_casy, fill=honeybees)) +
  geom_boxplot(color="grey60", outlier.shape=NA, width=0.5) +
  geom_point(position=position_jitterdodge(jitter.width=0.15), shape=21, color="grey60", size=2, aes(fill=honeybees)) +
  labs(y="pollen remaining (%) \n", x="study site") +
  scale_fill_manual(values=c("yellow", "black")) +
  #scale_x_discrete(labels=c("open", "bagged")) +
  guides(fill="none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_cowplot(14)

pol_plot

pol_plot
pdf("plots/Fig4_pol_plot.pdf")
print(pol_plot)
dev.off()

