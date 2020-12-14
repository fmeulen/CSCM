# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(transport)
theme_set(theme_minimal(base_size = 12))  

wd = getwd()

### make figs for binerror probabilities
dd <- read_csv(paste0(wd,"/out/binprobs.csv")) 
obs <- read_csv("out/observations.csv")
trace <- read.csv("out/tracepcn.csv") %>% gather(key="parameter", value="y", theta1, theta10,theta11,logtau)

p <-dd %>%  ggplot(aes(x, y, fill=value)) + geom_raster(hjust=0,vjust=0) + facet_wrap(~method, scales="free")+
  scale_fill_gradient(low="lightblue", high="black") + xlab("") + ylab("") + theme(aspect.ratio=2/1) #+ ggtitle("Error")
p

p_loss <- dd %>% filter(method %in% c("D","LNGL")) %>%  ggplot(aes(x, y, fill=loss)) + geom_raster(hjust=0,vjust=0) + facet_wrap(~method)+
  scale_fill_gradient2() + xlab("") + ylab("") + theme(aspect.ratio=1) #+ ggtitle("Error")
p_loss



pdf("./out/postmean.pdf",width=8,height=4)
  show(p)
dev.off()

pdf("./out/losspostmean.pdf",width=8,height=4)
show(p_loss)
dev.off()

ddtrue <- dd %>% filter(method=='true')
p4 <- ggplot() +geom_raster(data=ddtrue,aes(x=x,y=y,fill=value),hjust=0,vjust=0)  +
  geom_point(data=obs, aes(x=t, y=y, colour=yobserved), size=1.1) +
  scale_fill_gradient2() +
  ggtitle("true bin probabilities and observations")+
  xlab("") + ylab("") + theme(aspect.ratio=1)
p4
pdf("./out/obs.pdf",width=5,height=5)
p4   
dev.off()


p5 <- trace %>% ggplot() + geom_path(aes(x=iterate,y=y)) + facet_wrap(~parameter, scales="free") + ylab("")
pdf("./out/traceplotspcn.pdf",width=8, height=6)
p5  
dev.off()

