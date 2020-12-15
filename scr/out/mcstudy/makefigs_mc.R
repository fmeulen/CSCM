# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(transport)
library(scales)
library(ggthemes)

theme_set(theme_bw(base_size = 11))  

wd = getwd()

mc <- read_csv("mcstudy.csv", col_types = cols(samplesize = col_character()))

## plots for MC-study

pmc <- mc %>% gather(key="prior", value=x, D, LNGL) %>% ggplot(aes(x=x)) + geom_boxplot(fill="lightblue") +
  facet_wrap(~ prior+samplesize, nrow=1) + xlab("Wasserstein distance") + 
  coord_flip() +   scale_y_continuous(breaks = NULL)
pmc
pdf("wasserstein1.pdf",width=8, height=3)
pmc
dev.off()


pmc2 <- mc %>% ggplot(aes(x=D, y = LNGL, colour=samplesize)) + geom_point(size=0.8) + 
  geom_abline(aes(slope=1,intercept=0))  +
  theme(aspect.ratio=1) +  scale_colour_colorblind() 
pmc2


pdf("wasserstein2.pdf",width=4, height=4)
pmc2
dev.off()


