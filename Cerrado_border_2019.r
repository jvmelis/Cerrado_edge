## Setting workspace and reading data
if(!require(readxl)){install.packages("readxl")}       # read_excel
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(lme4)){install.packages('lme4')}           # (g)lmer
if(!require(MASS)){install.packages("MASS")}           # glm.nb
if(!require(car)){install.packages("car")}             # vif
if(!require(lsmeans)){install.packages('lsmeans')}     # lsmeans
if(!require(glmmTMB)){install.packages("glmmTMB")}     # Zero-Inflated & Hurdle Mixed Models
if(!require(bbmle)){install.packages('bbmle')}         # AICtab
if(!require(vegan)){install.packages("vegan")}         # vegdist, adonis
if(!require(iNEXT)){install.packages("iNEXT")}         # rarefaction curves
if(!require(betapart)){install.packages("betapart")}   # betapart


setwd("C:/Users/User/Documents/ARTIGOS - Meus/2 MS - Cerrado edge effect on lianas/2019 Austral Ecology")
# source('clean_data.r')
dados<-read.csv("clean_data.csv")

#######################################################################
## 1. (G)LM(M): BA & abund ~ dist*geo
#######################################################################
# A. Lianas Abundance
mod1 <- glm(lianas_ab~geo*dist, dados, family = poisson)
mod2 <- MASS::glm.nb(lianas_ab~geo*dist, dados)
par(mfrow=c(2,2))
plot(mod1)
summary(mod2) # overdispersion (159 / 31)
plot(mod2)
AIC(mod1,mod2)
anova(mod2, test='Chisq')
car::Anova(mod2)
(ls_lianas<-lsmeans(mod2, ~geo*dist))
exp(as.data.frame(ls_lianas)$lsmean)
exp(as.data.frame(ls_lianas)$SE)
lsmeans::lsmeans(mod2, specs='geo')
lsmeans::lsmeans(mod2, specs='dist')

# B. Trees Abundance
mod3 <- glm(trees_ab~geo*dist, dados, family=poisson)
plot(mod3)
summary(mod3) # ns
anova(mod3, test='Chisq')
car::Anova(mod3)
(ls_trees<-lsmeans::lsmeans(mod3, ~geo*dist))
plot(ls_trees)
exp(as.data.frame(ls_trees)$lsmean)
exp(as.data.frame(ls_trees)$SE)
lsmeans::lsmeans(mod3, specs='geo')

# C. Lianas BA
lianas<-read_xlsx('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior'))
mod4 <- lmer(log(liana_BA)~geo*dist+(1|Plot), lianas)
plot(mod4)
summary(mod4) # ns
anova(mod4, 'Chisq')
car::Anova(mod4)
(liana_intrageo<-lsmeans(mod4, ~ dist|geo)) # ns between `dist`
plot(liana_intrageo) # ns

# D. Trees BA
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
trees<-read_xlsx('cerrado_dados.xlsx', sheet = 'trees')
trees<-left_join(trees,hosts, by='ind')

mod5 <- trees %>% mutate(tree_BA = pi*((Diameter/2)^2)) %>%
  with(lmer(log(tree_BA)~geo*dist+(1|Plot)))
plot(mod5)
summary(mod5) # ns
anova(mod5, 'Chisq')
car::Anova(mod5)
(tree_intra<-lsmeans(mod5, ~ dist*geo))
plot(tree_intra) # in the south was ***, east not so much
exp(as.data.frame(tree_intra)$lsmean)
exp(as.data.frame(tree_intra)$SE)

# Synthesis: 
# a) Differences in lianas abundace (south more)
# b) Differences in trees abundace (south more)
# c) NO differences in lianas BA
# d) Differences in trees BA (edge specially in the south has higher BA).

#######################################################################
## 2. Environmental differences between edges and interiors
#######################################################################
# A. betadisp e adonis (PERMANOVA): environmental differences?
env_data<- dplyr::select(dados, c(Tree_cov, Regen, Palm_cov, Bare_soil, 
                       native, exotic, brom, PAR,SOM, Al, Mn)) %>% 
  decostand(method = 'range')
psych::pairs.panels(env_data) # Tree_cov, Bare_soi, exotic, brom show no variation
env_data<- dplyr::select(dados, c(Regen, Palm_cov,  
                                  native,  PAR,SOM, Al, Mn)) %>% 
  decostand(method = 'range')
psych::pairs.panels(env_data)

env_dist <- vegdist(env_data)
adonis(env_dist ~ geo*dist, data=dados,
       permutations = 999)

# 1) There is a difference in the location of the samples (i.e. the average community composition), 
# 2) There is a difference in the dispersion of the samples (i.e. the variability in the community composition), or 
# 3) There is a difference in both the location and the dispersion.
# https://www.researchgate.net/post/How_should_I_correctly_manage_PERMANOVA_for_factors_with_interactions

dispersion<-betadisper(env_dist, group=dados$local) # Permdist test
permutest(dispersion) # it is the same as: anova(dispersion)
par(mfrow=c(1,1))
plot(dispersion)
boxplot(dispersion)
# Conclusion: there is a difference in dispersion AND location

env_mds<-metaMDS(env_data)
stressplot(env_mds)
df <- cbind(as.data.frame(scores(env_mds)), 
      dplyr::select(dados, c(geo,dist, local))) 

(env_fit <- envfit(env_mds, env_data, perm = 999))


df_biofit<-scores(env_fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

levels(df$local)<-c("EE","IE", "ES","IS")
rownames(df_biofit)<-c("Reg","Pal","Nat","PAR","SOM", "Al", "Mn")

ggplot()+
  stat_ellipse(data = df, aes(x=NMDS1,y=NMDS2, group=local),
               linetype=2)+
  geom_segment(data=df_biofit, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                  arrow = arrow(length = unit(0.2, "cm")))+
  geom_point(data=df, aes(NMDS1,NMDS2,shape=geo, 
                          fill = dist), size=3)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  theme_classic()+
  theme(legend.position = 'none')+
  geom_text(data = as.data.frame(
    rbind(df_biofit[-7,]*1.45,df_biofit[7,]*1.2)),
    aes(NMDS1, NMDS2, 
        label = rownames(df_biofit)), size=4) +
  geom_text(data = data.frame(local=c("EE","ES","IE","IS"),
                              NMDS1=c(-.7,-.4,-.2,.3),
                              NMDS2=c(-.4,-.6,-.5,-.6)),
            aes(NMDS1,NMDS2, label=local), size=7)
# ggsave('Fig2_EnvironmentalNMDS.png')

# based on env_fit, we chose: PAR, Regen, native,  SOM, Al, Mn

# A. Lianas Abundance
mod1 <- glm(lianas_ab~PAR+Regen+native+SOM+Al+Mn, 
            dados, family = poisson)
mod2 <- MASS::glm.nb(lianas_ab~PAR+Regen+native+SOM+Al+Mn, 
                     dados)
par(mfrow=c(2,2))
plot(mod1)
summary(mod1) # overdispersion (159 / 31)
plot(mod2)
AIC(mod1,mod2)
summary(mod2) # ns
sqrt(car::vif(mod2)) # no significant vif but `Al`

# B. Trees Abundance
mod3 <- glm(trees_ab~PAR+Regen+native+SOM+Al+Mn, 
            dados, family=poisson)
plot(mod3)
summary(mod3) # overdispersion
mod4 <- MASS::glm.nb(trees_ab ~ PAR+Regen+native+SOM+Al+Mn, 
                     dados)
plot(mod4)
AIC(mod3,mod4)
summary(mod4) # ns
sqrt(car::vif(mod4)) # no significant vif

# C. Lianas BA - it was not significant

# D. Trees BA
geral<-inner_join(trees,dados[-1,], by="Plot") %>% 
  mutate(tree_BA = pi*(Diameter/2)) %>%
  dplyr::select(Plot, tree_BA,PAR,Regen,native,SOM,Al,Mn)
mod5 <- lmer(log(tree_BA) ~ PAR+Regen+native+SOM+Al+Mn+
               (1|Plot), 
           geral)
plot(mod5)
summary(mod5) # ns 
anova(mod5, 'Chisq')

sqrt(car::vif(mod5)) # Al has high vif value (1.9) & Mn - Regen are redundant in NMDS.

mod6 <- lmer(log(tree_BA) ~ PAR+Regen+native+SOM+
               (1|Plot), 
             geral)
plot(mod6)
summary(mod6) # ns 
anova(mod6, 'Chisq') # ns

#######################################################################
## 3. N_lianas x trees_BA * local
#######################################################################
trees %>% group_by(geo,dist, Plot) %>%
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
zip <- glmmTMB(n_lianas ~ local + (1 | Plot), 
                data = trees, family=poisson,
                ziformula=~1)
zip0 <- glmmTMB(n_lianas ~ local + (1 | Plot), 
               data = trees, family=poisson,
               ziformula=~local)
zip1 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                data = trees, family=poisson,
                ziformula=~1)
zip2 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                data = trees, family=poisson,
                ziformula=~local)
zip3 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                data = trees, family=poisson,
                ziformula=~1)
zip4 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                data = trees, family=poisson,
                ziformula=~trees)
bbmle::AICtab(zip,zip0,zip1,zip2,zip3, zip4) # zip4

# B. Zero Inflated models NBinomial
zinb <- glmmTMB(n_lianas ~ local + (1 | Plot), 
                 data = trees, family=nbinom1,
                 ziformula=~1)
zinb0 <- glmmTMB(n_lianas ~ local + (1 | Plot), 
                 data = trees, family=nbinom1,
                 ziformula=~local)
zinb1 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                data = trees, family=nbinom1,
                ziformula=~1)
zinb2 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                data = trees, family=nbinom1,
                ziformula=~local)
zinb3 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                data = trees, family=nbinom1,
                ziformula=~1)
zinb4 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                data = trees, family=nbinom1,
                ziformula=~local)
bbmle::AICtab(zinb, zinb0,zinb1,zinb2,zinb3,zinb4) # zinb3

# C. Hurdle Poisson
hurdleP0 <- glmmTMB(n_lianas ~ local  + (1 | Plot), 
               data = trees, ziformula=~1,
               family=truncated_poisson)
hurdleP1 <- glmmTMB(n_lianas ~ local  + (1 | Plot), 
                    data = trees, ziformula=~local,
                    family=truncated_poisson)
hurdleP2 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~1,
                    family=truncated_poisson)
hurdleP3 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~local ,
                    family=truncated_poisson)
hurdleP4 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~1,
                    family=truncated_poisson)
hurdleP5 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~local,
                    family=truncated_poisson)
bbmle::AICtab(hurdleP0,hurdleP1,hurdleP2,hurdleP3,hurdleP4,hurdleP5) # hurdleP3

# D. Hurdle Negative Binomial
hurdleNB0 <- glmmTMB(n_lianas ~ local  + (1 | Plot), 
                    data = trees, ziformula=~1,
                    family=truncated_nbinom1)
hurdleNB1 <- glmmTMB(n_lianas ~ local  + (1 | Plot), 
                    data = trees, ziformula=~local,
                    family=truncated_nbinom1)
hurdleNB2 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~1,
                    family=truncated_nbinom1)
hurdleNB3 <- glmmTMB(n_lianas ~ local + log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~local ,
                    family=truncated_nbinom1)
hurdleNB4 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~1,
                    family=truncated_poisson)
hurdleNB5 <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                    data = trees, ziformula=~local,
                    family=truncated_nbinom1)
bbmle::AICtab(hurdleNB0,hurdleNB1,hurdleNB2,hurdleNB3,hurdleNB4,hurdleNB5) # hurdleNB1

# E. Final
bbmle::AICtab(zip4, zinb3,hurdleP3,hurdleNB1) # zinb3

summary(zinb3)
contrast(lstrends(zinb3, var="tree_BA", "local"))
contrast(lsmeans(zinb3, ~local|tree_BA))
plot(lstrends(zinb3, var="tree_BA", "local"))
plot(lsmeans(zinb3, ~local|tree_BA)) # trees w/ same BA in interior show < lianas

ggplot()+ # View Hurdle Model (Red = LOGIT, Blue = NB)
  geom_point(data=mutate(trees,zero=ifelse(n_lianas>0,"non_zero", "zero")), 
               aes(y=n_lianas, x=log(tree_BA), color=zero))+
  scale_color_manual(values=c("black", "red"))+
  geom_smooth(data=filter(trees,n_lianas>0),
              aes(y=n_lianas, x=log(tree_BA),  group=local),
              method = "glm.nb", se=F)+
  geom_smooth(data=mutate(trees,zero=ifelse(n_lianas>0, 1, 0)),
              aes(y=zero, x=log(tree_BA),  group=local),
              method = "glm", se=F, method.args=list(family='binomial'), 
              color="blue", linetype=2)+
  facet_grid(geo~dist)+theme_bw()

ggplot()+ # View ZINB Model (Red = LOGIT, Blue = NB)
  geom_point(data=mutate(trees,zero=ifelse(n_lianas>0,"non_zero", "zero")), 
             aes(y=n_lianas, x=log(tree_BA), color=zero))+
  scale_color_manual(values=c("black", "red"))+
  geom_smooth(data=filter(trees,n_lianas>0),
              aes(y=n_lianas, x=log(tree_BA),  group=local),
              method = "glm.nb", se=T)+
  geom_hline(data = group_by(trees,geo,dist) %>%
               filter(n_lianas>0) %>%
               summarise(n_lianas=mean(n_lianas)),
             aes(yintercept=n_lianas) )+
  facet_grid(geo~dist)+theme_bw()

# Simulation (Brooks et al. 2016)
sims<-simulate(zinb3, seed = 1, nsim = 1000)
simdatlist <- lapply(sims, function(n_lianas){
  cbind(n_lianas, trees[,c('tree_BA', 'local', 'Plot')])
})
simdatsums <- lapply(simdatlist, function(x){
  plyr::ddply(x, ~log(tree_BA)+local, summarize,
        absence=mean(n_lianas==0),
        mu=mean(n_lianas))
})
ssd<-do.call(rbind, simdatsums)

ggplot(ssd, aes(y=absence, x=`log(tree_BA)`))+theme_bw()+
  geom_point()+
  facet_wrap(~local)+
  ylab("Probability that lianas are not climbing")+
  xlab("log(treeBA)")

########################################################################
# 4. Edge effect on species richness and composition of lianas and trees
########################################################################
##### Rarefaction curves
######################
## LIANAS  
df_lianas <- lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="BL"|site=="IL", 'East','South'),
         dist=ifelse(site=="BL"|site=="BS", 'Edge','Interior'))

ggplot(df_lianas, aes(x=x, y=y))+theme_bw()+
  geom_point(shape=16, data=df_lianas[which(df_lianas$method=="observed"),]) +
  geom_line(data=df_lianas[which(df_lianas$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_color_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,150,20))+
  ggtitle(" ")+facet_grid(dist~.)+
  xlab("Number of Individuals of lianas (ind.plot-1)")+ 
  ylab("Liana Species Richness")+
  guides(color= F, shape = F, 
           fill=guide_legend(override.aes = 
                               list(alpha=.5, color='black'), 
                             title = ""))
## TREES
df_trees <- trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="east_edge"|site=="east_interior", 'East','South'),
         dist=ifelse(site=="south_edge"|site=="east_edge", 'Edge','Interior'))

ggplot(df_trees, aes(x=x, y=y))+theme_bw()+
  geom_point(shape=16, data=df_trees[which(df_trees$method=="observed"),]) +
  geom_line(data=df_trees[which(df_trees$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_color_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,900,100))+
  ggtitle(" ")+facet_grid(dist~.)+
  xlab("Number of Individuals of trees (ind.plot-1)")+ 
  ylab("Tree Species Richness")+
  guides(color= F, shape = F, 
         fill=guide_legend(override.aes = 
                             list(alpha=.5, color='black'), 
                           title = ""))


######################
##### PCoA
#####################
## Lianas
lianas.t<-lianas %>% dplyr::select(Species, Plot) %>%
  with(table(Species,Plot)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,Freq,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, Plot)) %>% t()

locality <- c(rep("EE",9), rep("ES",10),
              rep("IE",5), rep("IS",8))
edges <- c(rep("Edge",19), rep("Interior",13))
gps <-c(rep("East",9), rep("South",10),
       rep("East",5), rep("South",8))
lianas_dist<-beta.pair(lianas.t, index.family="jaccard")

bd_local<-betadisper(lianas_dist[[3]],locality)
bd_edges<-betadisper(lianas_dist[[3]],edges)
bd_gps<-betadisper(lianas_dist[[3]],gps)
par(mfrow=c(1,3))
plot(bd_local); plot(bd_edges); plot(bd_gps)
anova(bd_local) # * IS

## TREES
trees.t<-trees %>% dplyr::select(Species, Plot) %>%
  with(table(Species,Plot)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,Freq,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, Plot)) %>% t()
local_trees <- c(rep("EE",10), rep("ES",10),
              rep("IE",8), rep("IS",8))
edges_trees <- c(rep("Edge",20), rep("Interior",16))
gps_trees <-c(rep("East",10), rep("South",10),
        rep("East",8), rep("South",8))
trees_dist<-beta.pair(trees.t, index.family="jaccard")

bd_local<-betadisper(trees_dist[[3]],local_trees)
bd_edges<-betadisper(trees_dist[[3]],edges_trees)
bd_gps<-betadisper(trees_dist[[3]],gps_trees)
par(mfrow=c(1,3))
plot(bd_local); plot(bd_edges); plot(bd_gps)
anova(bd_local) # * ES
anova(bd_gps)
anova(bd_edges)

#####################
##### PERMANOVA
#####################
# see: https://chrischizinski.github.io/rstats/vegan-ggplot2/
## Lianas
lianas_data <- lianas.t %>% 
  as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Plot, Species)) %>% as.data.frame.matrix() 
(mod <- metaMDS(lianas_data))
par(mfrow=c(1,1))
plot(mod, type="t")
plot(mod); ordihull(mod,groups=dados$local,draw="polygon",col="grey90",label=T)

adonis(lianas_data ~ geo*dist, data=dados, permutations = 999)
adonis(lianas_data ~ local, data=dados, permutations = 999)
envfit(mod, dados[,c('geo','dist', 'local')])

data.scores <- as.data.frame(scores(mod)) 
data.scores$site <- rownames(data.scores)
data.scores$grp <- locality
data.scores$geo <- gps
data.scores$dist <- edges

species.scores <- as.data.frame(scores(mod, "species"))
species.scores$species <- rownames(species.scores)

ggplot() + 
  stat_ellipse(data = data.scores, aes(x=NMDS1,y=NMDS2, 
                                       group=grp),linetype=3)+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,
                                  shape=geo, fill=dist),size=2) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  coord_equal()+ 
  theme_classic() + 
  guides(shape = guide_legend(
    override.aes = list(fill="black")),
    fill = guide_legend(override.aes =list(shape=22, size=3))) +
  geom_text(data = data.frame(local=c("IE", "ES","EE","IS"),
                              NMDS1=c( 4.2,-1.0, 2.4,-1.6),
                              NMDS2=c(-1.0, 1.5, 1.5, 0.8)),
            aes(NMDS1,NMDS2, label=local), size=4)

## Trees
trees_data <- trees.t %>% 
  as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Plot, Species)) %>% as.data.frame.matrix() 

(mod_trees <- metaMDS(trees_data))
par(mfrow=c(1,1))
plot(mod_trees, type="t")
plot(mod_trees); ordihull(mod_trees,
                          groups=dados$local,
                          draw="polygon",col="grey90",label=T)

adonis(trees_data ~ edges_trees*gps_trees, permutations = 999) # stres > 0.2?
adonis(trees_data ~ local_trees, data=dados, permutations = 999)
envfit(mod_trees, data.frame(local_trees,edges_trees,gps_trees))

data.scores <- as.data.frame(scores(mod_trees)) 
data.scores$site <- rownames(data.scores)
data.scores$grp <- local_trees
data.scores$geo <- gps_trees
data.scores$dist <- edges_trees

species.scores <- as.data.frame(scores(mod_trees, "species"))
species.scores$species <- rownames(species.scores)

ggplot() + 
  stat_ellipse(data = data.scores, aes(x=NMDS1,y=NMDS2, 
                                       group=grp),linetype=3)+
  geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2,
                                   shape=geo, fill=dist),size=2) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  coord_equal()+ 
  theme_classic() + 
  guides(shape = guide_legend(
    override.aes = list(fill="black")),
    fill = guide_legend(override.aes =list(shape=22, size=3))) +
  geom_text(data = data.frame(local=c("IE", "ES","EE","IS"),
                              NMDS1=c( 0.2, 1.0, -0.5, 0.8),
                              NMDS2=c( 0.2, 0.6,  0.5, 0.3)),
            aes(NMDS1,NMDS2, label=local), size=4)


#####################
##### Beta partiotining
#####################
## Lianas
lianas_table <- lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>% t()

summary(lianas.core <- betapart.core(lianas_table))
lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  group_by(border) %>% dplyr::summarise(sum(Freq))
lianas.core$shared
lianas.core$not.shared

pair.lianas <- beta.pair(lianas_table)
par(mfrow=c(1,3))
plot(hclust(pair.lianas$beta.sim, method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sim]), line=0.3)
plot(hclust(pair.lianas$beta.sne, 
            method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sne]), line=0.3)
plot(hclust(pair.lianas$beta.sor, method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sor]), line=0.3)

pair.lianas$beta.sne/pair.lianas$beta.sor
pair.lianas$beta.sim

## Trees
trees_table <- trees %>% 
  mutate(border = paste0(geo,"_",dist)) %>%
  dplyr::select(Species, border)  %>%
  with(table(Species, border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>% t()

summary(trees.core <- betapart.core(trees_table))
trees %>% mutate(border = paste0(geo,"_",dist)) %>%
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  group_by(border) %>% dplyr::summarise(sum(Freq))
trees.core$shared
trees.core$not.shared

pair.trees <- beta.pair(trees_table)
par(mfrow=c(1,3))
plot(hclust(pair.trees$beta.sim, method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sim]), line=0.3)
plot(hclust(pair.trees$beta.sne, 
            method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sne]), line=0.3)
plot(hclust(pair.trees$beta.sor, method="average"),
     hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sor]), line=0.3)

pair.trees$beta.sne/pair.lianas$beta.sor
pair.trees$beta.sim