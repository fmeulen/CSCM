# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(tidyverse)
theme_set(theme_bw(base_size = 12))

wd = getwd()

### make figs for binerror probabilities
dir_binprob <- read_csv(paste0(wd,"/out/Dirichletbinprob.csv")) %>% mutate(Dirichlet=pest-ptrue)
gl_binprob <- read_csv(paste0(wd,"/out/graphLaplacianbinprob.csv")) %>% mutate(GraphLaplacian=pest-ptrue)

d <- dir_binprob %>% dplyr::select(x,y,Dirichlet) %>%
    mutate(GraphLaplacian =gl_binprob$GraphLaplacian) %>% 
    gather(key="estimator",value="error",Dirichlet, GraphLaplacian)

d

mincol_lim <- min(d$error); maxcol_lim <- max(d$error)
p <- ggplot(data=d,aes(x, y, fill=error)) + geom_tile() + facet_wrap(~estimator)+
  scale_fill_gradient2(limits=c(mincol_lim, maxcol_lim))
p  

pdf("./out/binproberror.pdf",width=8,height=4)
  show(p)
dev.off()


##### trace plots for graphlap
gloutdf <- read_csv("out/gloutdf.csv")
pgl <- gloutdf %>% gather(key="variable",value="value",-iterate) %>%
  ggplot() + geom_line(aes(x=iterate,y=value)) + facet_wrap(~variable,scales="free")


#### trace plots for dirichlet
diroutdf <- read_csv("out/diroutdf.csv")
pdir <- diroutdf %>%  filter(binID %in% c("[1,1]","[1,2]","[1,3]","[2,3]","[3,3]","[5,5]"))%>%
  ggplot() + geom_line(aes(x=iterate,y=w)) + facet_wrap(~binID,scales="free")



pdf("./out/trace_gl.pdf",width=8,height=4)
show(pgl)
dev.off()
pdf("./out/trace_dir.pdf",width=8,height=4)
show(pdir)
dev.off()
