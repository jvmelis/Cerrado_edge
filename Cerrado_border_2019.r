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
set.seed(42)
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
df <- cbind(as.data.frame(scores(env_mds)), 
      dplyr::select(dados, c(geo,dist, local))) 

(env_fit <- envfit(env_mds, env_data, perm = 999))


df_biofit<-scores(env_fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

levels(df$local)<-c("EE","IE", "ES","IS")
rownames(df_biofit)<-c("Reg","Pal","Nat","PAR","SOM", "Al", "Mn")

Fig2a<-ggplot()+
  stat_ellipse(data = df, aes(x=NMDS1,y=NMDS2, group=local),
               linetype=2)+
  geom_point(data=df, aes(NMDS1,NMDS2,shape=geo, 
                          fill = dist), size=3)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  theme_classic()+
  theme(legend.position = 'none')+
  geom_text(data = data.frame(local=c("EE","SE","EI","SI"),
                              NMDS1=c(-.7, -.4, -.2, 0.4),
                              NMDS2=c(0.3, 1.0, .5, -.6)),
            aes(NMDS1,NMDS2, label=local), size=4)
Fig2a

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
(trend_rel<-contrast(lstrends(zinb3, var="tree_BA", "local")))
(mean_rel<-contrast(lsmeans(zinb3, ~local|tree_BA)))
plot(trend_rel)
plot(mean_rel) # trees w/ same BA in interior show < lianas

trees %>% filter(zero=='With Lianas ') %>% group_by(local) %>%
  summarise(sum_lianas=sum(n_lianas), N = n()) %>%
  mutate(Mean_rel = sum_lianas/N)
data.frame(Local = as.data.frame(mean_rel)$contrast,
           Mean = exp(as.data.frame(mean_rel)$estimate),
           SE = exp(as.data.frame(mean_rel)$SE))

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

ssd<-ssd %>% 
  mutate(geo = ifelse(local=='east_interior'|local=='east_edge', 
                      'East face', 'South face'),
         dist = ifelse(local=='south_edge'|local=='east_edge', 
                       'Edge', 'Interior'),
         tree_BA = exp(`log(tree_BA)`))
head(ssd)

fig3<-ggplot(ssd, aes(y=absence, x=tree_BA))+theme_bw()+
  geom_point()+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", 
                                                   function(x) 10^x),
                     labels = scales::trans_format("log10", 
                                                   scales::math_format(10^.x)))+
  facet_grid(dist~geo)+
  ylab("Probability that lianas are not climbing")+
  xlab(expression("Basal area of host tree(cm"^2*")"))
fig3
# ggsave('Fig3_TreeBAxLianaAb.jpg')

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

## TREES
df_trees <- trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="east_edge"|site=="east_interior", 'East','South'),
         dist=ifelse(site=="south_edge"|site=="east_edge", 'Edge','Interior'))

# source('figura_rarefacao.r')

df_lianas[which(df_lianas$method=="observed"),]
lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>% 
  # t()%>%  rarefy(14)
  iNEXT(q=0,datatype = "abundance", size=c(14,14))

df_trees[which(df_trees$method=="observed"),]
trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance", size=c(548,548))
  # t()%>%  rarefy(548)

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
set.seed(42)
(mod_lianas <- metaMDS(lianas_data))
par(mfrow=c(1,1))
plot(mod_lianas, type="t")
plot(mod_lianas); ordihull(mod_lianas,
  groups=dados$local,draw="polygon",col="grey90",label=T)

envfit(mod_lianas, dados[,c('geo','dist', 'local')])
lianas_dist <- vegdist(lianas_data)
set.seed(42)
adonis(lianas_dist ~ geo*dist, 
       data=data.frame(geo=gps,dist=edges),
       permutations = 999)

dispersion<-betadisper(lianas_dist, group=locality) # Permdist test
permutest(dispersion) # it is the same as: anova(dispersion)
par(mfrow=c(1,1))
plot(dispersion)

data.scores <- as.data.frame(scores(mod_lianas)) 
data.scores$site <- rownames(data.scores)
data.scores$grp <- locality
data.scores$geo <- gps
data.scores$dist <- edges

species.scores <- as.data.frame(scores(mod_lianas, "species"))
species.scores$species <- rownames(species.scores)

Fig2b<-ggplot() + 
  stat_ellipse(data = data.scores, 
               aes(x=NMDS1,y=NMDS2, group=grp),
               linetype=3)+
  geom_point(data=data.scores,
             aes(x=NMDS1,y=NMDS2,shape=geo, fill=dist),
             size=2) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  #coord_equal()+ 
  theme_classic()+
  theme(legend.position = 'none')+
  geom_text(data = data.frame(local=c("EI", "SE","EE","SI"),
                              NMDS1=c( 2.2,-0.2,-2.5,-1.0),
                              NMDS2=c( 1.1, 1.8, 0.2, 1.1)),
            aes(NMDS1,NMDS2, label=local), size=4)

Fig2b
resulta <- data.frame(envfit(mod_lianas, lianas_data)$vectors$arrows,
           p = as.vector(envfit(mod_lianas, 
                                lianas_data)$vectors$pvals)) 
resulta[which(resulta$p<0.05),]
## Trees
trees_data <- trees.t %>% 
  as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Plot, Species)) %>% as.data.frame.matrix() 

(mod_trees <- metaMDS(trees_data)) # stress>0.2
par(mfrow=c(1,1))
plot(mod_trees, type="t")
plot(mod_trees); ordihull(mod_trees,
                          groups=local_trees,
                          draw="polygon",col="grey90",label=T)

set.seed(42)
adonis(trees_data ~ edges_trees*gps_trees, permutations = 999)
trees_dist <- vegdist(trees_data)
dispersion<-betadisper(trees_dist, group=gps_trees) # Permdist test
permutest(dispersion) # it is the same as: anova(dispersion)
plot(dispersion)

resulta<-envfit(mod_trees, trees_data)

resulta<-data.frame(resulta$vectors$arrows,
                    r2=resulta$vectors$r,
                    p=resulta$vectors$pvals)

# write.csv(round(resulta[which(resulta$p<0.05),],3),'resulta.csv')

data.scores <- as.data.frame(scores(mod_trees)) 
data.scores$site <- rownames(data.scores)
data.scores$grp <- local_trees
data.scores$geo <- gps_trees
data.scores$dist <- edges_trees

species.scores <- as.data.frame(scores(mod_trees, "species"))
species.scores$species <- rownames(species.scores)

Fig2c<-ggplot() + 
  stat_ellipse(data = data.scores, aes(x=NMDS1,y=NMDS2, 
                                       group=grp),linetype=3)+
  geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2,
                                   shape=geo, fill=dist),size=2) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_grey(start = 1, end = .5)+
  #coord_equal()+ 
  theme_classic()+
  theme(legend.position = 'none')+
  guides(shape = guide_legend(
    override.aes = list(fill="black")),
    fill = guide_legend(override.aes =list(shape=22, size=3))) +
  geom_text(data = data.frame(local=c("EI", "SE","EE","SI"),
                              NMDS1=c( 0.0, 1.0, -0.5, 0.8),
                              NMDS2=c( 0.3, 0.5,  0.45, 0.3)),
            aes(NMDS1,NMDS2, label=local), size=4)
Fig2c
cowplot::plot_grid(Fig2a,Fig2b,Fig2c, labels = c('a)','b)','c)'), ncol = 1)
# ggsave('Fig2_NMDS.jpg')

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
diag(lianas.core$shared) # total species
lianas.core$shared
lianas.core$not.shared

pair.lianas <- beta.pair(lianas_table)

pair.lianas$beta.sor # = SNE + SIM
pair.lianas$beta.sne
pair.lianas$beta.sim

mean(pair.lianas$beta.sor)
mean(pair.lianas$beta.sim)
mean(pair.lianas$beta.sne)

t(lianas_table) %>% as.data.frame()%>% spread(border,Freq) %>%
  filter(BL!=0&BS!=0&IL!=0&IS!=0) # All areas
  filter(BL!=0&BS==0&IL==0&IS==0) # exclusive to East Edge
  #filter(BL==0&BS!=0&IL==0&IS==0) # exclusive to South Edge
  #filter(BL==0&BS==0&IL!=0&IS==0) # exclusive to East Interior = 0
  #filter(BL==0&BS==0&IL==0&IS!=0) # exclusive to South East = 0

## Trees
trees_table <- trees %>% 
  mutate(border = paste0(geo,"_",dist)) %>%
  dplyr::select(Species, border)  %>%
  with(table(Species, border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>% t()

summary(trees.core <- betapart.core(trees_table))
diag(trees.core$shared)
trees.core$shared
trees.core$not.shared

pair.trees <- beta.pair(trees_table)
pair.trees$beta.sor
pair.trees$beta.sne
# pair.trees$beta.sne/pair.lianas$beta.sor
pair.trees$beta.sim

mean(pair.trees$beta.sor)
mean(pair.trees$beta.sim)
mean(pair.trees$beta.sne)

t(trees_table) %>% as.data.frame() %>% spread(border,Freq) %>%
  filter(east_edge!=0&south_edge!=0&east_interior!=0&south_interior!=0) # All = 43 spp
 # filter(east_edge!=0&south_edge==0&east_interior==0&south_interior==0) # exclusive to East Edge = 7 spp
 # filter(east_edge==0&south_edge!=0&east_interior==0&south_interior==0) # exclusive to South Edge = 12 spp
 # filter(east_edge==0&south_edge==0&east_interior!=0&south_interior==0) # exclusive to East Interior = 8 spp
 # filter(east_edge==0&south_edge==0&east_interior==0&south_interior!=0) # exclusive to South East = 12 spp
