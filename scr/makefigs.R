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


d <- dd %>%     gather(key="estimate",value="probability",Dirichlet, graphLaplacian) %>%
  mutate(error = ptrue-probability, relerror=error/ptrue)
p <-d %>%  ggplot(aes(x, y, fill=error)) + geom_raster(hjust=0,vjust=0) + facet_wrap(~estimate)+
   scale_fill_gradient2() + xlab("") + ylab("") + theme(aspect.ratio=1) + ggtitle("Error")
p  
pdf("./out/binproberror.pdf",width=8,height=4)
  show(p)
dev.off()

d %>%  ggplot(aes(x, y, fill=relerror)) + geom_raster(hjust=0,vjust=0) + 
  facet_wrap(~estimate)+
  scale_fill_gradient2() + xlab("") + ylab("") + theme(aspect.ratio=1) + ggtitle("Relative error")

d %>% group_by(estimate) %>% summarise(mrelerr=mean(relerror))

dd %>% mutate(relerr = abs((Dirichlet-ptrue)/(graphLaplacian-ptrue))-1) %>%
  ggplot(aes(x, y, fill=relerr)) + geom_raster(hjust=0,vjust=0) +scale_fill_gradient2() + xlab("") +
  ylab("") + theme(aspect.ratio=1) + ggtitle("Ratio of relative errors")


p2a <-d %>%  filter(estimate=="Dirichlet") %>% ggplot(aes(x, y, fill=error)) + geom_raster(hjust=0,vjust=0)+
  scale_fill_gradient2() + ggtitle("Dirichlet") + xlab("") + ylab("") + theme(aspect.ratio=1)

p2b <- d %>%  filter(estimate=="graphLaplacian") %>%ggplot(aes(x, y, fill=error)) + geom_raster(hjust=0,vjust=0)+
       scale_fill_gradient2() + ggtitle("graphLaplacian")+ xlab("") + ylab("") + theme(aspect.ratio=1)
grid.arrange(p2a,p2b, ncol=2)



pdf("./out/binproberror2.pdf",width=8,height=4)
grid.arrange(p2a,p2b, ncol=2)
dev.off()


p3a <-d %>%  filter(estimate=="Dirichlet") %>% ggplot(aes(x, y, fill=probability)) + geom_raster(hjust=0,vjust=0) +
  scale_fill_gradient2() + ggtitle("Dirichlet") + xlab("") + ylab("") + theme(aspect.ratio=1)

p3b <- d %>%  filter(estimate=="graphLaplacian") %>%ggplot(aes(x, y, fill=probability)) + geom_raster(hjust=0,vjust=0)+
  scale_fill_gradient2() + ggtitle("graphLaplacian")+ xlab("") + ylab("") + theme(aspect.ratio=1)
grid.arrange(p3a,p3b, ncol=2)


# to superimpose points to geom_tile, asthetics need to have the same name, so to plot t, it needs to be renamed x (confusing, so first drop x)
# obs2 <- obs %>% select(-x) %>% mutate(x=t)
p4 <- ggplot() +geom_raster(data=d,aes(x=x,y=y,fill=ptrue),hjust=0,vjust=0)  +
  geom_point(data=obs, aes(x=t, y=y, colour=yobserved), size=1.1) +
  scale_fill_gradient2() +
  ggtitle("true bin probabilities and observations")+
  xlab("") + ylab("") + theme(aspect.ratio=1)

pdf("./out/obs.pdf",width=5,height=5)
p4   
dev.off()


p5 <- trace %>% ggplot() + geom_path(aes(x=iterate,y=y)) + facet_wrap(~parameter, scales="free") + ylab("")
pdf("./out/traceplotspcn.pdf",width=8, height=6)
p5  
dev.off()


D <- dd %>%     gather(key="estimate",value="probability",Dirichlet, graphLaplacian, ptrue) 
  
p <- D %>%  ggplot(aes(x, y, fill=probability)) + geom_raster(hjust=0,vjust=0) + facet_wrap(~estimate)+
  scale_fill_gradient2() + xlab("") + ylab("") + theme(aspect.ratio=1) 
p
