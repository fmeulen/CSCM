# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(transport)
theme_set(theme_bw(base_size = 12))

wd = getwd()

### make figs for binerror probabilities
dd <- read_csv(paste0(wd,"/out/binprobs.csv")) 



d <- dd %>%     gather(key="estimate",value="probability",Dirichlet, graphLaplacian) %>% mutate(error = ptrue-probability)
mincol_lim <- min(d$error); maxcol_lim <- max(d$error)
p <-d %>%  ggplot(aes(x, y, fill=error)) + geom_tile() + facet_wrap(~estimate, ncol=1)+
  scale_fill_gradient2(limits=c(mincol_lim, maxcol_lim)) + xlab("") + ylab("")
p  
pdf("./out/binproberror.pdf",width=8,height=4)
  show(p)
dev.off()


mincol_lim <- min(d$val); maxcol_lim <- max(d$val)
p2 <-d %>%  filter(estimate=="Dirichlet") %>%ggplot(aes(x, y, fill=error)) + geom_tile() +
  scale_fill_gradient2() + ggtitle("Dirichlet") + xlab("") + ylab("")

p3 <- d %>%  filter(estimate=="graphLaplacian") %>%ggplot(aes(x, y, fill=error)) + geom_tile()+
       scale_fill_gradient2() + ggtitle("graphLaplacian")+ xlab("") + ylab("")
grid.arrange(p2,p3)


  scale_fill_gradient2(limits=c(mincol_lim, maxcol_lim)) + xlab("") + ylab("")
p2  

pdf("./out/binprobest.pdf",width=8,height=4)
show(p2)
dev.off()

