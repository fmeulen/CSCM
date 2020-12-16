# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(transport)
library(scales)
library(ggthemes)

theme_set(theme_bw(base_size = 12))  

wd = getwd()

dirs <- c("./mnsmall", "./mnmedium", "./mnlarge")
app <- c("_small", "_medium", "_large")
for (i in 1:3)
{
  outdir = dirs[i]
  ### make figs for binerror probabilities
  dd <- read_csv(paste0(outdir,"/binprobs.csv"))
  obs <- read_csv(paste0(outdir,"/observations.csv"))
  
  p <-dd %>%  ggplot(aes(x, y, fill=binprob)) + geom_raster(hjust=0,vjust=0) + facet_wrap(~method, scales="free")+
    scale_fill_gradient(low="lightblue", high="black") + xlab("") + ylab("") + theme(aspect.ratio=1/1) #+ ggtitle("Error")
  p
  pdf(paste0(outdir,"/postmean",app[i],".pdf"),width=8,height=4)
  show(p)
  dev.off()
  
  ddtrue <- dd %>% filter(method=='true')
  p4 <- ggplot() +geom_raster(data=ddtrue,aes(x=x,y=y,fill=binprob),hjust=0,vjust=0)  +
    geom_point(data=obs, aes(x=t, y=y, colour=yobserved), size=1.1) +
    scale_fill_gradient(low="lightblue", high="black") +
    #ggtitle("true bin probabilities and observations")+
    xlab("") + ylab("") + theme(aspect.ratio=1)
  p4
  pdf(paste0(outdir,"/obs.pdf"),width=5,height=5)
  show(p4)   
  dev.off()
}  
  