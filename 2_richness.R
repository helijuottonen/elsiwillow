# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
# https://github.com/helijuottonen/elsiwillow

# 2: OTU richness of bacteria and fungi

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# RUN FIRST: 0_basic_processing.R

#library(vegan)
library(ggplot2) #v.3.3.5
library(cowplot) #v.1.1.1

# richness (OTU number) for bacteria
# note: this is a simple calculation of OTU numbers based on median subsampled data
# I have also tested repeated rarefaction based on the lowest read number across samples and results don't change

# diversity indices for bacteria (note: only richness used in the manuscript)
shannonb <- diversity(elotub.r.all, index="shannon")
richnessb <- specnumber(elotub.r.all)
bdiv.df <- as.data.frame(shannonb)
colnames(bdiv.df) <- "shannon"
bdiv.df$richness <- richnessb

# diversity indices for fungi
shannonf <- diversity(elotuf.r.all, index="shannon")
richnessf <- specnumber(elotuf.r.all)

fdiv.df <- as.data.frame(shannonf)
colnames(fdiv.df) <- "shannon"
fdiv.df$richness <- richnessf

bdiv.df$treat <- elmetab$treat
fdiv.df$treat <- elmetaf$treat
bdiv.df$honeybees <- elmetab$honeybees
fdiv.df$honeybees <- elmetaf$honeybees
bdiv.df$site <- elmetab$site
fdiv.df$site <- elmetaf$site
bdiv.df$pollen_remain <- elmetab$pollen_remain
fdiv.df$pollen_remain <- elmetaf$pollen_remain
bdiv.df$pollen_removed <- elmetab$pollen_removed
fdiv.df$pollen_removed <- elmetaf$pollen_removed

# removing nas
bdiv.df.na <- na.omit(bdiv.df)
fdiv.df.na <- na.omit(fdiv.df)

### richness vs. bagging: Fig. 3

# defining colours
treat_col <- c("#b5e61d", "#ffc90e")

# defining treatment names
treat.labels <- c(ctrl="open", ps="bagged")

# bacteria
brb <- ggplot(bdiv.df, aes(x=treat, y=richness, fill=treat)) +
  geom_boxplot(color="black", outlier.shape=NA, width=0.5) +
  geom_point(position=position_jitterdodge(jitter.width=0.15), shape=21, color="black", size=2, aes(fill=treat)) +
  labs(y="bacterial richness", x=NULL) +
  scale_fill_manual(values=treat_col) +
  scale_x_discrete(labels=c("open", "bagged")) +
  guides(fill="none") +
  ylim(0, 420) +
  theme_cowplot(14)
 
brb

# fungi

frb <- ggplot(fdiv.df, aes(x=treat, y=richness, fill=treat)) +
  geom_boxplot(color="black", outlier.shape=NA, width=0.5) +
  geom_point(position=position_jitterdodge(jitter.width=0.15), shape=21, color="black", size=2, aes(fill=treat)) +
  labs(y="fungal richness", x=NULL) +
  scale_fill_manual(values=treat_col) +
  scale_x_discrete(labels=c("open", "bagged")) +
  guides(fill="none") +
  ylim(0, 420) +
  theme_cowplot(14)

frb

# both plots together
rich.plot <- plot_grid(brb, frb, labels=c("a)", "b)"), label_fontface="plain", ncol=2, hjust=c(0.1, 0.1)) + theme(plot.margin = unit(c(0,0,0,0.2), "cm"))
rich.plot <- plot_grid(brb, frb, labels=c("a)", "b)"), label_fontface="plain")

rich.plot
pdf("plots/Fig3_richplot.pdf")
print(rich.plot)
dev.off()

### regressions for Fig. 5

# basic formula in regressions is lm(y ~ x)

#### regression richness vs. pollen remaining

# bacteria
# 'pollen remain' log transformed
brp <- lm(richness~log(pollen_remain), data=bdiv.df.na)
# plots to check if the model looks ok
plot(brp)
# displaying the model parameters
summary(brp)

#Call:
#  lm(formula = richness ~ log(pollen_remain), data = bdiv.df.na)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-141.599  -55.537    6.899   50.427  168.621 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          330.19      61.58   5.362 2.58e-06 ***
#  log(pollen_remain)   -36.39      15.25  -2.387   0.0212 *  
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 75.72 on 46 degrees of freedom
#Multiple R-squared:  0.1102,	Adjusted R-squared:  0.09083 
#F-statistic: 5.696 on 1 and 46 DF,  p-value: 0.02118

# -> regression equation y=330-36.4x, p=0.021

# fungi
# 'pollen remain' log transformed
frp <- lm(richness~log(pollen_remain), data=fdiv.df.na)
# plots to check if the model looks ok
plot(frp)
# displaying the model parameters
summary(frp)

#Call:
#  lm(formula = richness ~ log(pollen_remain), data = fdiv.df.na)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-79.058 -37.569  -6.856  34.559 110.126 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          324.95      40.84   7.957 3.43e-10 ***
#  log(pollen_remain)   -23.18      10.11  -2.292   0.0265 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 50.22 on 46 degrees of freedom
#Multiple R-squared:  0.1025,	Adjusted R-squared:  0.08299 
#F-statistic: 5.253 on 1 and 46 DF,  p-value: 0.02653

# regression equation: y=328-23.8x, p=0.023

### plots for richness vs. pollen: Fig. 5

# bacteria
bdp <- ggplot(bdiv.df.na, aes(x=pollen_remain, y=richness)) + 
  geom_smooth(method = "lm", se = FALSE, color="grey") +
  geom_point(pch=21, color="black", fill="yellow", size=2)+
  labs(y="bacterial richness", x="pollen remaining (%) \n") +
  ylim(0, 400) +
  scale_x_log10() + 
  theme_cowplot(11) +
  annotate("text", x= 105, y=5, label="y=330-36.4x, P=0.021")

bdp

#fungi
fdp <- ggplot(fdiv.df.na, aes(x=pollen_remain, y=richness)) + 
  geom_smooth(method = "lm", se = FALSE, color="grey") +
  geom_point(pch=21, color="black", fill="yellow", size=2)+
  labs(y="fungal richness", x="pollen remaining (%) \n") +
  ylim(0, 400) +
  scale_x_log10() + 
  theme_cowplot(11) +
  annotate("text", x=105, y= 5, label="y=328-23.8x, P=0.023")

fdp

# both plots together
rich_poll.plot <- plot_grid(bdp, fdp, labels=c("a)", "b)"), label_fontface="plain", ncol=1, hjust=c(0.1, 0.1)) + theme(plot.margin = unit(c(0,0,0,0.2), "cm"))

rich_poll.plot
pdf("plots/Fig5_rich_poll_plot_uusi.pdf")
print(rich_poll.plot)
dev.off()

#### richness values in Table S3:

bdiv_s3 <- aggregate(richness ~site + treat, data=bdiv.df, mean) 
bdiv_s3_sd <- aggregate(richness ~site + treat, data=bdiv.df, sd)
bdiv_s3$sd <- bdiv_s3_sd$richness

fdiv_s3 <- aggregate(richness ~site + treat, data=fdiv.df, mean) 
fdiv_s3_sd <- aggregate(richness ~site + treat, data=fdiv.df, sd)
fdiv_s3$sd <- fdiv_s3_sd$richness
  





