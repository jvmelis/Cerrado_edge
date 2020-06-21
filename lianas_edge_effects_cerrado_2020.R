rm(list=ls())
version #
# platform       x86_64-w64-mingw32          
# svn rev        77875                       
# version.string R version 3.6.3 (2020-02-29)
# nickname       Holding the Windsock 
## Setting workspace and reading data
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(lme4)){install.packages('lme4')}           # (g)lmer
if(!require(MASS)){install.packages("MASS")}           # glm.nb
if(!require(car)){install.packages("car")}             # vif
if(!require(emmeans)){install.packages('emmeans')}     # lsmeans
if(!require(glmmTMB)){install.packages("glmmTMB")}     # Zero-Inflated & Hurdle Mixed Models
if(!require(bbmle)){install.packages('bbmle')}         # AICtab
if(!require(vegan)){install.packages("vegan")}         # vegdist, adonis
if(!require(iNEXT)){install.packages("iNEXT")}         # rarefaction curves
if(!require(betapart)){install.packages("betapart")}   # betapart
if(!require(cowplot)){install.packages("cowplot")}     # plot_grid
lrtest<- function(mod_null, mod_full){
  df_null<-attr(logLik(mod_null), "df")
  df_full<-attr(logLik(mod_full), "df")
  df<-abs(df_null - df_full)
  null<-logLik(mod_null)[1]
  full<-logLik(mod_full)[1]
  chisq<-2*abs(full - null)
  pval<-1-pchisq(df=df,q=chisq)
  resulta<-data.frame(df=c(df_null, df_full,df),
                      LL= c(null,full,NA),
                      Chisq=c(NA,NA,chisq),
                      pval=c(NA,NA,pval))
  return(resulta)
}

# setwd(choose.dir())
# source('clean_data.r')
dados<-read.csv("clean_data.csv")
########################################
### 1. Effect of edge on abundance and basal area of lianas and trees
## A. Lianas Abundance
mod1 <- glm(lianas_ab~geo*dist, dados, family = poisson)
mod2 <- MASS::glm.nb(lianas_ab~geo*dist, dados)
par(mfrow=c(2,2))
plot(mod1)
summary(mod2) # overdispersion (159 / 31)
plot(mod2)
AIC(mod1,mod2)
anova(mod1, mod2, test='Chisq')
lrtest(mod1, mod2) # mod2 is better
mod2_null <- MASS::glm.nb(lianas_ab~1, dados)
lrtest(mod2_null, mod2) # mod2 is better
mod2_geo <- MASS::glm.nb(lianas_ab~geo, dados)
lrtest(mod2_geo, mod2) # mod2 is better
mod2_dist <- MASS::glm.nb(lianas_ab~dist, dados)
lrtest(mod2_dist, mod2) # mod2 is better
mod2_int <- MASS::glm.nb(lianas_ab~geo:dist, dados)
(ls_lianas<-lsmeans(mod2, ~geo*dist))
exp(as.data.frame(ls_lianas)$lsmean)
exp(as.data.frame(ls_lianas)$SE)
emmeans::lsmeans(mod2, specs='geo')
emmeans::lsmeans(mod2, specs='dist')
car::Anova(mod2)

# B. Trees Abundance
mod3 <- glm(trees_ab~geo*dist, dados, family=poisson)
summary(mod3) # overdispersion
mod4 <- MASS::glm.nb(trees_ab~geo*dist, dados)
lrtest(mod3,mod4)
AIC(mod3,mod4)
mod4_null <- MASS::glm.nb(trees_ab~1, dados)
mod4_geo <- MASS::glm.nb(trees_ab~geo, dados)
mod4_dist <- MASS::glm.nb(trees_ab~dist, dados)
lrtest(mod4_null,mod4)
lrtest(mod4_geo,mod4_null)
lrtest(mod4_dist, mod4) # ns - only distance is better (edge/interior)
(ls_trees<-lsmeans::lsmeans(mod4, ~ geo*dist))
plot(ls_trees) # to visualize
exp(as.data.frame(ls_trees)$lsmean)
exp(as.data.frame(ls_trees)$SE)
emmeans::lsmeans(mod4_dist, specs='dist')
plot(mod4_dist)
car::Anova(mod4)

## Synthesis
source('Cerrado_syno.r')
# A) Differences in lianas abundace (south more) /50m2
set.seed(42)
group_by(lianas, border,Plot) %>% 
  summarise(n=n()) %>%
  mutate(n=ifelse(is.na(n),0,n))%>%
  ungroup() %>%   group_by(border) %>%
  summarise(mean=mean_cl_boot(n)[[1]],
            min=mean_cl_boot(n)[[2]],
            max=mean_cl_boot(n)[[3]]) %>% 
  mutate(ratio=mean/max(min))
# B) Differences in trees abundace (south more) /50m2
trees %>% mutate(border = paste(geo,dist)) %>%
  group_by(border,Plot) %>%
  summarise(n=n()) %>% 
  mutate(n=ifelse(is.na(n),0,n))%>%
  ungroup() %>% group_by(border) %>%
  summarise(mean=mean_cl_boot(n)[[1]],
            min=mean_cl_boot(n)[[2]],
            max=mean_cl_boot(n)[[3]]) %>% 
  mutate(ratio=mean/min(mean))

# C) Lianas BA
mod5 <- lmer(log(liana_BA)~geo*dist+(1|Plot), lianas)
summary(mod5) # ns
plot(mod5)
mod5_null <- lmer(log(liana_BA)~1+(1|Plot), lianas)
mod5_geo <- lmer(log(liana_BA)~geo+(1|Plot), lianas)
mod5_dist <- lmer(log(liana_BA)~dist+(1|Plot), lianas)
lrtest(mod5, mod5_null) # ns
lrtest(mod5_geo, mod5_null) # ns
lrtest(mod5_dist, mod5_null) # ns
(liana_intrageo<-lsmeans(mod5, ~ dist|geo)) # to visualize
plot(liana_intrageo) # ns
car::Anova(mod5)

# D) Tree BA
mod6 <- lmer(log(tree_BA)~geo*dist+(1|Plot), trees)
summary(mod6) # ok
plot(mod6)
mod6_null <- lmer(log(tree_BA)~1+(1|Plot), trees)
mod6_geo <- lmer(log(tree_BA)~geo+(1|Plot), trees)
mod6_dist <- lmer(log(tree_BA)~dist+(1|Plot), trees)
lrtest(mod6, mod6_null) # <0.01
lrtest(mod6_geo, mod6_null) # ns
lrtest(mod6_dist, mod6_null) # <0.0001 - edge-interior is important
(tree_intrageo<-lsmeans(mod6, ~ dist|geo)) # to visualize
plot(tree_intrageo) 
car::Anova(mod6_dist)

## Synthesis: 
# C) NO differences in lianas BA
group_by(lianas, border,Plot) %>%
  mutate(BA=ifelse(is.na(liana_BA),0,liana_BA))%>%
  ungroup() %>% group_by(border) %>% 
  summarise(mean=mean_cl_boot(BA)[[1]],
          min=mean_cl_boot(BA)[[2]],
          max=mean_cl_boot(BA)[[3]]) %>% 
  mutate(ratio=mean/min(mean))

# D) Differences in trees BA (edge specially in the south has higher BA).
trees %>% group_by(border,Plot) %>%
  summarise(sum_BA = sum(tree_BA, na.rm=T)) %>% 
  mutate(BA.ha = ifelse(is.na(sum_BA),0,.2*sum_BA))%>%
  ungroup() %>% group_by(border) %>%
  summarise(mean=mean_cl_boot(BA.ha)[[1]],
            min=mean_cl_boot(BA.ha)[[2]],
            max=mean_cl_boot(BA.ha)[[3]]) %>% 
  mutate(ratio=mean/min(mean))
trees<-trees %>% mutate(geo = ifelse(geo=="east","East","South"))
FigS1a<-ggplot(lianas, aes(y=log10(liana_BA), x=dist))+
  geom_violin(fill="grey", alpha=0.4)+facet_grid(.~geo)+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "red", width=0.2)+
  scale_y_continuous( labels = scales::math_format(10^.x)) +
  xlab("")+theme_bw()+
  ylab(expression(paste("Basal area (",cm^2, ")")))+
  annotation_logticks(sides="lr")
FigS1b<-ggplot(trees, aes(y=log10(tree_BA/10000), x=dist))+
  geom_violin(fill="grey", alpha=0.4)+facet_grid(.~geo)+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "red", width=0.2)+
  scale_y_continuous( labels = scales::math_format(10^.x)) +
  xlab("")+theme_bw()+
  ylab(expression(paste("Basal area (",m^2, ")")))+
  annotation_logticks(sides="lr")+
  scale_x_discrete(labels=c("edge" = "Edge", "interior" = "Interior")) 
FigS1<-cowplot::plot_grid(FigS1a,FigS1b, ncol = 1, labels = c("A)", "B)"))
FigS1
# ggsave("FigS1_violinplot_BA.jpg")

lista<-ls()
nlista<-lista %in% c("dados", "FigS1","lianas","lrtest","trees")
rm(list=c(lista[!nlista],"lista","nlista"))
########################################################################
#### 2. Edge effect on species richness and composition of lianas and trees
### A) Rarefactions
## LIANAS  
lianas %>% filter(!is.na(Family)) %>% 
  mutate(morpho=paste(Family, Species),
         local = paste(geo,dist)) %>% with(table(morpho, local))

df_lianas <- lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="BL"|site=="IL", 'East','South'),
         dist=ifelse(site=="BL"|site=="BS", 'Edge','Interior'))

df_lianas %>% filter(method=="observed")%>% mutate(dif=y.upr-y) %>% dplyr::select(site,y,dif)

lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance",q=0, size = c(10:20))%>% fortify(type=1) %>% 
  filter(x==14) %>% mutate(dif=y.upr-y) %>% dplyr::select(site,y,dif)

## TREES
df_trees <- trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="East_edge"|site=="East_interior", 'East','South'),
         dist=ifelse(site=="South_edge"|site=="East_edge", 'Edge','Interior'))

df_trees%>% filter(method=="observed") %>% mutate(dif=y.upr-y) %>% dplyr::select(site,y,dif)

trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance", size = 520:550)%>% fortify(type=1) %>% 
  filter(x==548) %>% mutate(dif=y.upr-y) %>% dplyr::select(site,y,dif)

## Figures
p_liana<-ggplot(df_lianas, aes(x=x, y=y, group=geo))+theme_bw()+
  geom_point(shape=16,data=df_lianas[which(df_lianas$method=="observed"),]) +
  geom_line(data=df_lianas[which(df_lianas$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,150,20))+ facet_grid(.~dist)+
  xlab(expression("Liana abundance (ind.plot"^-1*")"))+ 
  ylab("Liana species richness \n (Species per plot)")+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_line(color = 'white'),
        panel.grid.major.x = element_line(color = 'white'),
        panel.grid.major.y = element_line(linetype = 2, colour = 'grey'))

p_tree<-ggplot(df_trees, aes(x=x, y=y, group=geo))+theme_bw()+
  geom_point(shape=16, data=df_trees[which(df_trees$method=="observed"),]) +
  geom_line(data=df_trees[which(df_trees$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0, end = 0.5)+
  scale_y_continuous(expand=c(0,0), breaks =seq(0,90,15))+
  scale_x_continuous(breaks = seq(0,900,150))+
  ggtitle(" ")+  facet_grid(.~dist)+
  xlab(expression("Tree abundance (ind.plot"^-1*")"))+ 
  ylab("Tree species richness \n (Species per plot)")+
  theme(legend.position = 'none',
        panel.grid.minor = element_line(color = 'white'),
        panel.grid.major.x = element_line(color = 'white'),
        panel.grid.major.y = element_line(linetype = 2, colour = 'grey'))
Fig1<-cowplot::plot_grid(p_tree,p_liana,nrow = 2)
Fig1
# ggsave('Fig1_rarefactioncurves.jpg', width = 6, height = 8) # Edited by Illus

### B) PERMANOVA
locality <- c(rep("EE",9), rep("ES",10),
              rep("IE",5), rep("IS",8))
edges <- c(rep("Edge",19), rep("Interior",13))
gps <-c(rep("East",9), rep("South",10),
        rep("East",5), rep("South",8))

## Lianas
lianas.t <- lianas %>% dplyr::select(Species, Plot) %>%
  with(table(Species,Plot)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0, Freq, NA)) %>% filter(!is.na(Freq)) %>%
  with(table(Species, Plot)) %>% t()

lianas_data <- lianas.t %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>%  filter(!is.na(Freq)) %>% 
  with(table(Plot, Species)) %>% as.data.frame.matrix() 

remove <- as.vector(which(apply(lianas_data,1, FUN=sum)==0)) # remove empty samples
lianas_dist <- vegan::vegdist(lianas_data[-remove,])
set.seed(42)
vegan::adonis(lianas_dist ~ geo*dist, 
              data = data.frame(geo = gps, dist = edges),
              permutations = 999)  # results in text

dispersion <- betadisper(lianas_dist, group=locality) # Permdist test
permutest(dispersion) # it is the same as: anova(dispersion)
par(mfrow=c(1,1))
plot(dispersion)

set.seed(42)
(mod_lianas <- metaMDS(lianas_data[-remove,]))
envfit(mod_lianas, dados[,c('geo','dist', 'local')])

data.scores <- as.data.frame(scores(mod_lianas)) %>% rownames_to_column() %>%  
  mutate(site = rowname,grp = locality, geo = gps, dist = edges)
species.scores <- as.data.frame(scores(mod_lianas, "species")) %>% 
  rownames_to_column("species")
Fig2a <- ggplot() + 
  stat_ellipse(data = data.scores, aes(x=NMDS1,y=NMDS2, group=grp), linetype=3)+
  geom_point(data=data.scores,aes(x=NMDS1, y=NMDS2, shape=geo, fill=dist), size=2) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  theme_classic() + theme(legend.position = 'none')+
  geom_text(data = data.frame(local=c("EI", "SE","EE","SI"),
                              NMDS1=c( 2.2,-0.2,-2.5,-1.0),
                              NMDS2=c( 1.1, 1.8, 0.2, 1.1)),
            aes(NMDS1, NMDS2, label = local), size=4)
Fig2a
data.frame(envfit(mod_lianas, lianas_data[-remove,])$vectors$arrows,
           r2=as.vector(envfit(mod_lianas, lianas_data[-remove,])$vectors$r),
           p = as.vector(envfit(mod_lianas, lianas_data[-remove,])$vectors$pvals)) %>% 
  rownames_to_column() %>% filter(p<0.05) # Table S3

## Trees
local_trees <- c(rep("EE",10), rep("ES",10),
                 rep("IE",8), rep("IS",8))
edges_trees <- c(rep("Edge",20), rep("Interior",16))
gps_trees <-c(rep("East",10), rep("South",10),
              rep("East",8), rep("South",8))
trees.t<- trees %>% dplyr::select(Species, Plot) %>%
  with(table(Species,Plot)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,Freq,NA)) %>% 
  filter(!is.na(Freq)) %>% with(table(Species, Plot)) %>% t()
trees_dist<-beta.pair(trees.t, index.family="jaccard")

trees_data <- trees.t %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% filter(!is.na(Freq)) %>% 
  with(table(Plot, Species)) %>% as.data.frame.matrix() 

set.seed(42)
(mod_trees <- metaMDS(trees_data)) # stress>0.2
vegan::adonis(trees_data ~ edges_trees*gps_trees, permutations = 999) # average
trees_dist <- vegan::vegdist(trees_data)
permutest(betadisper(trees_dist, group=gps_trees))   # Permdist test (Variance for east/south)
permutest(betadisper(trees_dist, group=edges_trees)) # Permdist test (Variance for edge/interior)
dispersion<-betadisper(trees_dist, group=local_trees)
permutest(dispersion) # it is the same as: anova(dispersion)
plot(dispersion)

data.frame(envfit(mod_trees, trees_data)$vectors$arrows,
           r2=as.vector(envfit(mod_trees, trees_data)$vectors$r),
           p = as.vector(envfit(mod_trees, trees_data)$vectors$pvals)) %>% 
  rownames_to_column() %>% filter(p<0.05) # Table S4

data.scores <- as.data.frame(scores(mod_trees)) %>% rownames_to_column("site") %>%  
  mutate(grp = local_trees, geo = gps_trees, dist = edges_trees)
species.scores <- as.data.frame(scores(mod_trees, "species")) %>% 
  rownames_to_column("species")

Fig2b <- ggplot() + stat_ellipse(data = data.scores, aes(x=NMDS1,y=NMDS2, 
                                                         group=grp),linetype=3)+
  geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2,
                                   shape=geo, fill=dist),size=2) +
  scale_shape_manual(values=c(21,24))+ scale_fill_grey(start = 1, end = .5)+
  theme_classic() + theme(legend.position = 'none')+
  guides(shape = guide_legend(
    override.aes = list(fill="black")),
    fill = guide_legend(override.aes =list(shape=22, size=3))) +
  geom_text(data = data.frame(local=c("EI", "SE","EE","SI"),
                              NMDS1=c( 0.0, 1.0, -0.5, 0.8),
                              NMDS2=c( 0.3, 0.5,  0.45, 0.3)),
            aes(NMDS1,NMDS2, label=local), size=4)
Fig2b

#### C) Beta partiotining
## Lianas
lianas_table <- lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>% with(table(Species, border)) %>% t()

pair.lianas <- beta.pair(lianas_table)
mean(pair.lianas$beta.sor) # beta sor
mean(pair.lianas$beta.sim) # beta SIM
mean(pair.lianas$beta.sne) # beta SNE

pair.lianas$beta.sor # = SNE + SIM
pair.lianas$beta.sne
pair.lianas$beta.sim

lianas.core <- betapart.core(lianas_table)
lianas.core$shared
lianas.core$not.shared

t(lianas_table) %>% as.data.frame()%>% spread(border,Freq) %>%
  #filter(BL!=0&BS!=0&IL!=0&IS!=0) # All areas
  filter(BL!=0&BS==0&IL==0&IS==0) # exclusive to East Edge
#filter(BL==0&BS!=0&IL==0&IS==0) # exclusive to South Edge
#filter(BL==0&BS==0&IL!=0&IS==0) # exclusive to East Interior = 0
#filter(BL==0&BS==0&IL==0&IS!=0) # exclusive to South East = 0

## Trees
trees_table <- trees %>% 
  mutate(border = paste0(geo,"_",dist)) %>% dplyr::select(Species, border)  %>%
  with(table(Species, border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>%  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>% t()

trees.core <- betapart.core(trees_table)
diag(trees.core$shared)
trees.core$shared
trees.core$not.shared

pair.trees <- beta.pair(trees_table)
mean(pair.trees$beta.sor)
mean(pair.trees$beta.sim)
mean(pair.trees$beta.sne)

pair.trees$beta.sor
pair.trees$beta.sne
pair.trees$beta.sim

source('VennDiagramCerrado.r')
Fig3
# ggsave('Fig3_VennDiagram.jpg')
lista<-ls()
nlista<-lista %in% c("dados","Fig1","Fig2a","Fig2b","Fig3","FigS1","lianas","trees","lrtest")
rm(list = c(lista[!nlista],"lista","nlista"))

###################################################################
### 3. Edge effects, Microenvironmental factors and the abundance and basal area of lianas and trees
## A. betadisp & adonis (PERMANOVA): environmental differences?
env_data<- dados %>% 
  dplyr::select(c(Tree_cov, Regen, Palm_cov, Bare_soil,
                  native, exotic, brom, PAR,SOM, Al, Mn)) %>% 
  decostand(method = 'range')
psych::pairs.panels(env_data) # Tree_cov, Bare_soi, exotic, brom show no variation

env_data<- dados %>% 
  dplyr::select(c(Regen, Palm_cov,native,  PAR,SOM, Al, Mn)) %>% 
  decostand(method = 'range')
psych::pairs.panels(env_data)

env_dist <- vegan::vegdist(env_data)
set.seed(42)
vegan::adonis(env_dist ~ geo*dist, data=dados,permutations = 999)

# 1) There is a difference in the location of the samples (i.e. the average community composition), 
# 2) There is a difference in the dispersion of the samples (i.e. the variability in the community composition), or 
# 3) There is a difference in both the location and the dispersion.
# https://www.researchgate.net/post/How_should_I_correctly_manage_PERMANOVA_for_factors_with_interactions

dispersion<-betadisper(env_dist, group=dados$local) # Permdist test
permutest(dispersion) # it is the same as: anova(dispersion)
par(mfrow=c(1,1))
plot(dispersion)
# Conclusion: there is a difference in dispersion AND location

env_mds<-metaMDS(env_data)
df <- cbind(as.data.frame(scores(env_mds)), 
      dplyr::select(dados, c(geo,dist, local))) 

data.frame(envfit(env_mds, env_data, perm = 999)$vectors$arrows,
           r2 = as.vector(envfit(env_mds, env_data, perm = 999)$vectors$r),
           p = as.vector(envfit(env_mds, env_data, perm = 999)$vectors$pvals)) %>% 
  rownames_to_column()  # Table S5

df_biofit<-scores(env_fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)
levels(df$local)<-c("EE","IE", "ES","IS")
rownames(df_biofit)<-c("Reg","Pal","Nat","PAR","SOM", "Al", "Mn")

Fig2c<-ggplot()+
  stat_ellipse(data = df, aes(x=NMDS1, y=NMDS2, group=local),linetype=2)+
  geom_point(data=df, aes(NMDS1,NMDS2,shape=geo, fill = dist), size=3)+
  scale_shape_manual(values=c(21,24))+ scale_fill_grey(start = 1, end = .5)+
  theme_classic() + theme(legend.position = 'none')+
  geom_text(data = data.frame(local=c("EE","SE","EI","SI"),
                              NMDS1=c(-.7, -.4, -.2, 0.4),
                              NMDS2=c(0.3, 1.0, .5, -.6)),
            aes(NMDS1,NMDS2, label=local), size=4)
Fig2c
Fig2<-cowplot::plot_grid(Fig2a, Fig2b, Fig2c, labels = c('a)','b)','c)'), ncol = 1)
Fig2
# ggsave('Fig2_NMDS.jpg')

lista<-ls()
nlista<-lista %in% c("dados","lianas","lrtest","trees")
rm(list=lista[!nlista] )
rm(list=c("lista","nlista"))
### based on env_fit, we chose: PAR, Regen, native,  SOM, Al, Mn
# A. Lianas Abundance
mod1 <- glm(lianas_ab~PAR+Regen+native+SOM+Al+Mn, 
            dados, family = poisson)
mod2 <- MASS::glm.nb(lianas_ab~PAR+Regen+native+SOM+Al+Mn, dados)

lrtest(mod1,mod2)
summary(mod1) # overdispersion (159 / 31)
AIC(mod1,mod2)
par(mfrow=c(2,2))
plot(mod2)
summary(mod2) # ns
sqrt(car::vif(mod2)) # no significant vif but `Al`
mod2_null <- MASS::glm.nb(lianas_ab~1,dados)
lrtest(mod2,mod2_null)
drop1(mod2) # ns
add1(mod2_null, scope = ~PAR+Regen+native+SOM+Al+Mn)
mod2_a<- MASS::glm.nb(lianas_ab~Mn, dados)
add1(mod2_a, scope = ~PAR+Regen+native+SOM+Al+Mn)
lrtest(mod2_a, mod2_null)
AIC(mod2_a)-AIC(mod2_null)
summary(mod2_a)

# B. Trees Abundance
mod3 <- glm(trees_ab~PAR+Regen+native+SOM+Al+Mn, dados, family=poisson)
summary(mod3) # overdispersion
mod4 <- MASS::glm.nb(trees_ab ~ PAR+Regen+native+SOM+Al+Mn, dados)
sqrt(car::vif(mod4)) # no significant vif
mod4_null <- MASS::glm.nb(trees_ab ~ 1, dados)
lrtest(mod4,mod4_null) # ns
drop1(mod4) # ns
add1(mod4_null,scope = ~PAR+Regen+native+SOM+Al+Mn )
mod4_a <- MASS::glm.nb(trees_ab ~ Al, dados)
lrtest(mod4_a,mod4_null)
AIC(mod4_a)-AIC(mod4_null)
add1(mod4_a,scope = ~PAR+Regen+native+SOM+Al+Mn ) # ns

# C. Lianas BA - it was not significant
geral<-inner_join(lianas, dados, by="Plot") %>% 
  dplyr::select(Plot, liana_BA, PAR,Regen,native,SOM,Al,Mn)
mod5 <- lmer(log(liana_BA) ~ PAR+Regen+native+SOM+Al+Mn+(1|Plot), geral)
isSingular(mod5)
mod5_null<-lmer(log(liana_BA) ~ 1+(1|Plot), geral)
add1(mod5_null, scope = ~PAR+Regen+native+SOM+Al+Mn )
mod5_a <- lmer(log(liana_BA) ~ PAR+(1|Plot), geral)
lrtest(mod5_null, mod5_a)
add1(mod5_a, scope = ~PAR+Regen+native+SOM+Al+Mn )
AIC(mod5_a)-AIC(mod5_null)
plot(mod5_a)

# D. Trees BA
geral<-inner_join(trees, dados, by="Plot") %>% 
  dplyr::select(Plot, tree_BA, PAR,Regen,native,SOM,Al,Mn)
mod6 <- lmer(log(tree_BA) ~ PAR+Regen+native+SOM+Al+Mn+(1|Plot), geral)
mod6_null <- lmer(log(tree_BA) ~ 1 + (1|Plot), geral)
lrtest(mod6_null, mod6) # ns
sqrt(car::vif(mod6)) # Al has high vif value (1.9) & Mn - Regen are redundant in NMDS.
drop1(mod6)
add1(mod6_null, scope = ~PAR+Regen+native+SOM+Al+Mn) # ns

#######################################################################
### 3. N_lianas x trees_BA * local
hosts<-read_xlsx('cerrado_dados.xlsx', sheet = 'host_tree')
h1<-dplyr::select(hosts,c(id_1,n_1))
h2<-dplyr::select(hosts,c(id_2,n_2))
h3<-dplyr::select(hosts,c(id_3,n_3))
names(h1)<-names(h2)<-names(h3)<-c("ind","n")
hosts<-left_join(h1,h2, by="ind") %>% 
  left_join(h3,by='ind') %>%
  group_by(ind) %>% 
  dplyr::summarise(n_lianas=sum(sum(n.x,n.y, 
                                    na.rm = T), n, 
                                na.rm = T))

trees <- left_join(trees, hosts, by='ind')
lista <- ls()
nlista <- lista %in% c("dados","geral","lianas","lrtest","trees")
rm(list = lista[!nlista])
  
trees %>% group_by(geo, dist, Plot) %>%
  summarise(N_trees = n(), N_lianas = sum(n_lianas, na.rm=T)) %>%
  ggplot(aes(y=N_lianas, x=N_trees, group=dist, color=dist))+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se=F)+
  facet_grid(geo~.)+
  geom_point()+theme_bw()

trees<-trees %>% 
  mutate(n_lianas = ifelse(is.na(n_lianas),0,n_lianas)) %>%
  # mutate(AGB = exp(-3.3369+2.7635*log(Diameter)+0.4059*log(Height)+1.2439*log(WD))) %>% # Ribeiro et al. 2011 - m1
  # mutate(AGB = exp(-3.352+2.9853*LN(Diameter)+1.1855*LN(WD))) %>% # Ribeiro et al. 2011 - m4
  # mutate(AGB = exp(-3.3369+2.7635*LN(Diameter)+0.4059*LN(Height)+1.2439*LN(WD))) %>%  # 
  # mutate(AGB = exp(-2.6504+0.8713*LN((Diameter^2)*Height))) %>% # Melo et al. Unpublished
  # mutate(AGB = (28.77*(Diameter^2)*Height)/1000 ) %>% # Delliti et al. 
  mutate(zero = ifelse(n_lianas==0,"No Lianas","With Lianas "),
         tree_BA = pi*((Diameter/2)^2),
         local=paste0(geo,"_",dist))

trees %>%  ggplot(aes(n_lianas, fill=zero)) +
  facet_grid(dist~geo) + geom_hline(yintercept = 0)+
  scale_y_continuous(expand=c(0,0))+
  geom_histogram(color='black')+theme_classic()

# A. Zero Inflated models Poisson
zip_null <-glmmTMB(n_lianas ~ 
                     1 + (1|Plot), data = trees, family=poisson,
                   ziformula=~ 1)
zip_tree0 <- glmmTMB(n_lianas ~ 
                        log(tree_BA) + (1|Plot), data = trees, family=poisson,
                      ziformula=~ 1)
zip_tree1 <- glmmTMB(n_lianas ~ 
                        log(tree_BA) + (1|Plot), data = trees, family=poisson,
                      ziformula=~log(tree_BA))
zip_local0 <- glmmTMB(n_lianas ~ 
                 local + (1|Plot), data = trees, family=poisson,
                ziformula=~ 1)
zip_local1 <- glmmTMB(n_lianas ~ 
                  local + (1|Plot), data = trees, family=poisson,
               ziformula=~local)
zip_ad0 <- glmmTMB(n_lianas ~ 
                  local + log(tree_BA) + (1|Plot), data = trees, family=poisson,
                ziformula=~ 1)
zip_ad1 <- glmmTMB(n_lianas ~ 
                  local + log(tree_BA) + (1|Plot), data = trees, family=poisson,
                ziformula=~local)
zip_ad2 <- glmmTMB(n_lianas ~ 
                     local + log(tree_BA) + (1|Plot), data = trees, family=poisson,
                   ziformula=~log(tree_BA))
zip_int0 <- glmmTMB(n_lianas ~ 
                  local * log(tree_BA) + (1|Plot), data = trees, family=poisson,
                ziformula=~1)
zip_int1 <- glmmTMB(n_lianas ~ 
                      local * log(tree_BA) + (1|Plot), data = trees, family=poisson,
                ziformula=~local)
zip_int2 <- glmmTMB(n_lianas ~ 
                      local * log(tree_BA) + (1|Plot), data = trees, family=poisson,
                    ziformula=~log(tree_BA))
bbmle::AICtab(zip_null,
              zip_tree0, zip_tree1,
              zip_local0,zip_local1,
              zip_ad0,zip_ad1,zip_ad2,
              zip_int0, zip_int1, zip_int2) # zip_int2

# B. Zero Inflated models NBinomial
zinb_null <- glmmTMB(n_lianas ~ 
                       1 + (1|Plot), data = trees, family=nbinom1,
                 ziformula=~1)
zinb_tree0 <- glmmTMB(n_lianas ~ 
                       log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                 ziformula= ~ 1)
zinb_tree1 <- glmmTMB(n_lianas ~ 
                       log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                     ziformula= ~ log(tree_BA))
zinb_local0 <- glmmTMB(n_lianas ~ 
                       local + (1|Plot), data = trees, family=nbinom1,
                     ziformula= ~ 1)
zinb_local1 <- glmmTMB(n_lianas ~ 
                         local + (1|Plot), data = trees, family=nbinom1,
                       ziformula= ~ local)
zinb_ad0 <- glmmTMB(n_lianas ~ 
                         local + log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                       ziformula= ~ 1)
zinb_ad1 <- glmmTMB(n_lianas ~ 
                      local + log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                    ziformula= ~ local)
zinb_ad2 <- glmmTMB(n_lianas ~ 
                      local + log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                    ziformula= ~ log(tree_BA))
zinb_int0 <- glmmTMB(n_lianas ~ 
                      local*log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                    ziformula= ~ 1)
zinb_int1 <- glmmTMB(n_lianas ~ 
                      local*log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                    ziformula= ~ local)
zinb_int2 <- glmmTMB(n_lianas ~ 
                      local*log(tree_BA) + (1|Plot), data = trees, family=nbinom1,
                    ziformula= ~ log(tree_BA))
bbmle::AICtab(zinb_null, 
              zinb_tree0, zinb_tree1, 
              zinb_local0, zinb_local1,
              zinb_ad0, zinb_ad1, zinb_ad2,
              zinb_int0, zinb_int1, zinb_int2) # zinb_int2

# C. Hurdle Poisson
hurdleP_null <- glmmTMB(n_lianas ~ 
                      1  + (1|Plot), data = trees, family=truncated_poisson, 
               ziformula=~1)
hurdleP_tree0 <- glmmTMB(n_lianas ~ 
                          log(tree_BA)  + (1|Plot), data = trees, family=truncated_poisson, 
                        ziformula=~1)
hurdleP_tree1 <- glmmTMB(n_lianas ~ 
                           log(tree_BA)  + (1|Plot), data = trees, family=truncated_poisson, 
                         ziformula=~log(tree_BA))
hurdleP_local0 <- glmmTMB(n_lianas ~ 
                           local  + (1|Plot), data = trees, family=truncated_poisson, 
                         ziformula=~1)
hurdleP_local1 <- glmmTMB(n_lianas ~ 
                           local  + (1|Plot), data = trees, family=truncated_poisson, 
                         ziformula=~local)
hurdleP_ad0 <- glmmTMB(n_lianas ~ 
                            local +log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                          ziformula=~1)
hurdleP_ad1 <- glmmTMB(n_lianas ~ 
                            local  +log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                          ziformula=~local)
hurdleP_ad2 <- glmmTMB(n_lianas ~ 
                            local  +log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                          ziformula=~log(tree_BA))
hurdleP_in0 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                       ziformula=~1)
hurdleP_in1 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                       ziformula=~local)
hurdleP_in2 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_poisson, 
                       ziformula=~log(tree_BA))

bbmle::AICtab(hurdleP_null, 
              hurdleP_ad0, hurdleP_ad1, hurdleP_ad2,
              hurdleP_in0, hurdleP_in1, hurdleP_in2,
              hurdleP_tree0, hurdleP_tree1,
              hurdleP_local0, hurdleP_local1) # hurdleP_ad1

# D. Hurdle Negative Binomial
hurdleNB_null <- glmmTMB(n_lianas ~ 
                       local  + (1 | Plot), data = trees,family=truncated_nbinom1,
                     ziformula=~1)
hurdleNB_tree0 <- glmmTMB(n_lianas ~ 
                           log(tree_BA)  + (1|Plot), data = trees, family=truncated_nbinom1, 
                         ziformula=~1)
hurdleNB_tree1 <- glmmTMB(n_lianas ~ 
                           log(tree_BA)  + (1|Plot), data = trees, family=truncated_nbinom1, 
                         ziformula=~log(tree_BA))
hurdleNB_local0 <- glmmTMB(n_lianas ~ 
                            local  + (1|Plot), data = trees, family=truncated_nbinom1, 
                          ziformula=~1)
hurdleNB_local1 <- glmmTMB(n_lianas ~ 
                            local  + (1|Plot), data = trees, family=truncated_nbinom1, 
                          ziformula=~local)
hurdleNB_ad0 <- glmmTMB(n_lianas ~ 
                         local +log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~1)
hurdleNB_ad1 <- glmmTMB(n_lianas ~ 
                         local  +log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~local)
hurdleNB_ad2 <- glmmTMB(n_lianas ~ 
                         local  +log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~log(tree_BA))
hurdleNB_in0 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~1)
hurdleNB_in1 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~local)
hurdleNB_in2 <- glmmTMB(n_lianas ~ 
                         local*log(tree_BA) + (1|Plot), data = trees, family=truncated_nbinom1, 
                       ziformula=~log(tree_BA))

bbmle::AICtab(hurdleNB_null, 
              hurdleNB_ad0, hurdleNB_ad1, hurdleNB_ad2,
              hurdleNB_in0, hurdleNB_in1, hurdleNB_in2,
              hurdleNB_tree0, hurdleNB_tree1,
              hurdleNB_local0, hurdleNB_local1) # hurdleNB_local1


# E. Final
bbmle::AICtab(zip_int2, zinb_int2, hurdleP_ad1, hurdleNB_local1) # zinb_int2

summary(zinb_int2)
(trend_rel<-contrast(lstrends(zinb_int2, var="tree_BA", "local")))
(mean_rel<-contrast(lsmeans(zinb_int2, ~local|tree_BA)))
plot(trend_rel) # east_edge show less, south_interior show more (mean, but variance is large)
plot(mean_rel) # trees w/ same BA in interior show < lianas

trees %>% filter(zero=='With Lianas ') %>% group_by(local) %>%
  summarise(sum_lianas=sum(n_lianas), N = n()) %>%
  mutate(Mean_rel = sum_lianas/N)
data.frame(Local = as.data.frame(mean_rel)$contrast,
           Mean = exp(as.data.frame(mean_rel)$estimate),
           SE = exp(as.data.frame(mean_rel)$SE))

# Simulation (Brooks et al. 2016)
sims<-simulate(zinb_int2, seed = 1, nsim = 1000)
simdatlist <- lapply(sims, function(n_lianas){
  cbind(n_lianas, trees[,c('tree_BA', 'local', 'Plot')])
})
simdatsums <- lapply(simdatlist, function(x){
  plyr::ddply(x, ~log(tree_BA)+local, summarize,
        absence=mean(n_lianas==0),
        mu=mean(n_lianas))
})
ssd<-do.call(rbind, simdatsums)

ssd<-ssd %>% 
  mutate(geo = ifelse(local=='east_interior'|local=='east_edge', 
                      'East face', 'South face'),
         dist = ifelse(local=='south_edge'|local=='east_edge', 
                       'Edge', 'Interior'),
         tree_BA = exp(`log(tree_BA)`))
head(ssd)

fig4<-ggplot(ssd, aes(y=1-absence, x=tree_BA))+theme_bw()+
  geom_point()+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", 
                                                   function(x) 10^x),
                     labels = scales::trans_format("log10", 
                                                   scales::math_format(10^.x)))+
  facet_grid(geo~dist)+
  ylab("Probability that lianas are climbing")+
  xlab(expression("Basal area of host tree(cm"^2*")"))
fig4
# ggsave('Fig4_TreeBAxLianaAb.jpg')

### END