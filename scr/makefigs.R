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


d <- dd %>%     gather(key="estimate",value="probability",Dirichlet, graphLaplacian) %>% mutate(error = ptrue-probability)
mincol_lim <- min(d$error); maxcol_lim <- max(d$error)
p <-d %>%  ggplot(aes(x, y, fill=error)) + geom_tile() + facet_wrap(~estimate)+
  scale_fill_gradient2(limits=c(mincol_lim, maxcol_lim)) + xlab("") + ylab("") + theme(aspect.ratio=1)
p  
pdf("./out/binproberror.pdf",width=8,height=4)
  show(p)
dev.off()



p2a <-d %>%  filter(estimate=="Dirichlet") %>% ggplot(aes(x, y, fill=error)) + geom_tile() +
  scale_fill_gradient2() + ggtitle("Dirichlet") + xlab("") + ylab("") + theme(aspect.ratio=1)

p2b <- d %>%  filter(estimate=="graphLaplacian") %>%ggplot(aes(x, y, fill=error)) + geom_tile()+
       scale_fill_gradient2() + ggtitle("graphLaplacian")+ xlab("") + ylab("") + theme(aspect.ratio=1)
p2 <- grid.arrange(p2a,p2b, ncol=2)



pdf("./out/binproberror2.pdf",width=8,height=4)
grid.arrange(p2a,p2b, ncol=2)
dev.off()

