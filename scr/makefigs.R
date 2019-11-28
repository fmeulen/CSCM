# make figs for CSCM

# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(tidyverse)


dir_binprob <- read_csv("~/.julia/dev/CSCM/scr/out/Dirichletbinprob.csv") %>% mutate(Dirichlet=pest-ptrue)
gl_binprob <- read_csv("~/.julia/dev/CSCM/scr/out/graphLaplacianbinprob.csv") %>% mutate(GraphLaplacian=pest-ptrue)

binprobs <- dir_binprob %>% dplyr::select(x,y,Dirichlet) %>%
    mutate(GraphLaplacian =gl_binprob$GraphLaplacian) %>% 
    gather(key="estimator",value="value",Dirichlet, GraphLaplacian)

binprobs
