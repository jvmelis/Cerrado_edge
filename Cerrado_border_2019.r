## Setting workspace and reading data
if(!require(readxl)){install.packages("readxl")}       # read_excel
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(cowplot)){install.packages('cowplot')}     # graphics
# if(!require(ggalt)){install.packages("ggalt")}         # graphics
# if(!require(ggExtra)){install.packages('ggExtra')}     # graphics (ggMarginal)
# if(!require(lmerTest)){install.packages('lmerTest')}   # 
if(!require(lme4)){install.packages('lme4')}           # (g)lmer
if(!require(MASS)){install.packages("MASS")}           # glm.nb
if(!require(car)){install.packages("car")}             # vif
if(!require(lsmeans)){install.packages('lsmeans')}     # lsmeans
if(!require(glmmTMB)){install.packages("glmmTMB")}     # Zero-Inflated Mixed Models
if(!require(vegan)){install.packages("vegan")}         # vegdist, adonis
if(!require(iNEXT)){install.packages("iNEXT")}         # rarefaction curves
if(!require(betapart)){install.packages("betapart")}   # betapart


setwd("C:/Users/User/Documents/ARTIGOS - Meus/2 MS - Cerrado edge effect on lianas/2019 Austral Ecology")
# source('clean_data.r')
dados<-read.csv("clean_data.csv")

### Liana & Tree Abund and Basal Area
p_lianaBA<-ggplot(dados, aes(y=log(lianas_BA), x=dist))+
  facet_grid(.~geo)+geom_boxplot()+theme_bw()
p_treeBA<-ggplot(dados, aes(y=log(trees_BA), x=dist))+
  facet_grid(.~geo)+geom_boxplot()+theme_bw()
p_lianaAb<-ggplot(dados, aes(y=lianas_ab, x=dist))+
  facet_grid(.~geo)+geom_boxplot()+theme_bw()
p_treeAb<-ggplot(dados, aes(y=trees_ab, x=dist))+
  facet_grid(.~geo)+geom_boxplot()+theme_bw()
cowplot::plot_grid(p_lianaAb,p_lianaBA, p_treeAb, p_treeBA)

p<-ggplot(dados, aes(x=trees_ab, y=lianas_ab, 
                     shape=dist,color=geo, group=geo)) +
  geom_point()+ geom_smooth(se=F, method=glm, 
                            method.args=list(family='poisson'))+
  theme_bw() + theme(legend.position="bottom")
# ggMarginal(p, type="histogram", fill="white")

p+facet_grid(dist~.)

## 1. (G)LM(M): BA & abund ~ dist*geo
#################################### 
# A. Lianas Abundance
mod1 <- glm(lianas_ab~geo*dist, dados, family = poisson)
mod2 <- MASS::glm.nb(lianas_ab~geo*dist, dados)
par(mfrow=c(2,2))
plot(mod1)
summary(mod2) # overdispersion (159 / 31)
plot(mod2)
AIC(mod1,mod2)
summary(mod2) # *** `geo`
lsmeans::lsmeans(mod2, specs='geo')

# B. Trees Abundance
mod3 <- glm(trees_ab~geo*dist, dados, family=poisson)
plot(mod3)
summary(mod3) # ns
lsmeans::lsmeans(mod3, specs='geo')
TukeyHSD(aov(log(lianas_BA)~geo, dados))

# C. Lianas BA
lianas<-read_xlsx('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior'))
lianas %>% dplyr::select(border,geo) %>% table()
lianas %>% dplyr::select(border,dist) %>% table()
mod4 <- lmer(log(liana_BA)~geo*dist+(1|Plot), lianas)
plot(mod4)
summary(mod4) # ns
anova(mod4, 'Chisq')
(liana_intrageo<-lsmeans(mod4, ~ dist|geo)) # ns between `dist`
plot(liana_intrageo) # ns
(liana_intradist<-lsmeans(mod4, ~ geo|dist)) # ns between `geo`
plot(liana_intradist)

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
trees<-left_join(trees,hosts)

mod5 <- trees %>% mutate(tree_BA = pi*((Diameter/2)^2)) %>%
  with(lmer(log(tree_BA)~geo*dist+(1|Plot)))
plot(mod5)
summary(mod5) # ns
anova(mod5, 'Chisq')
(tree_intrageo<-lsmeans(mod5, ~ dist|geo))
plot(tree_intrageo) # in the south was ***, east not so much
(tree_intradist<-lsmeans(mod5, ~ geo|dist)) # ns between `geo`
plot(tree_intradist)


# Synthesis: 
# a) Differences in lianas abundace (south more)
# b) Differences in trees abundace (south more)
# c) NO differences in lianas BA
# d) Differences in trees BA (edge specially in the south has higher BA)


## 2. Environmental differences between edges and interiors


#################################
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
  guides(shape = guide_legend(
    override.aes = list(fill="black")),
         fill = guide_legend(
    override.aes =list(shape=22, size=3)))+
  geom_text(data = as.data.frame(
    rbind(df_biofit[-7,]*1.45,df_biofit[7,]*1.2)),
    aes(NMDS1, NMDS2, 
        label = rownames(df_biofit)), size=4) +
  geom_text(data = data.frame(local=c("EE","ES","IE","IS"),
                              NMDS1=c(-.7,-.4,-.2,.3),
                              NMDS2=c(-.4,-.6,-.5,-.6)),
            aes(NMDS1,NMDS2, label=local), size=7)

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

#######################
## 3. N_trees_w_lianas x N_trees_wo_lianas * local

trees %>% group_by(geo,dist, Plot) %>%
  summarise(N = n(), N_lianas = sum(!is.na(n_lianas))) %>%
  ggplot(aes(y=N_lianas, x=N, group=dist, color=dist))+
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
  mutate(zero = ifelse(n_lianas==0,"No Lianas","With Lianas "))

trees %>%  ggplot(aes(n_lianas, fill=zero)) +
  facet_grid(dist~geo) + geom_hline(yintercept = 0)+
  scale_y_continuous(expand=c(0,0))+
  geom_histogram(color='black')+theme_classic()

# Zero Inflated models (Poisson and NBinomial)

ZIP <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
               data = trees,
               family = poisson)

ZINB <- glmmTMB(n_lianas ~ local * log(tree_BA) + (1 | Plot), 
                data = trees,
                family = nbinom2)
anova(ZINB,ZIP)

summary(ZINB)

contrast(lstrends(ZINB, var="tree_BA", "local"))
plot(lstrends(ZINB, var="tree_BA", "local"))
plot(lsmeans(ZINB, ~tree_BA|local)) # trees in interior show -3 lianas/tree

ggplot()+ # View ZINB Model (Red = LOGIT, Blue = NB)
  geom_point(data=mutate(trees,zero=ifelse(n_lianas>0,"non_zero", "zero")), 
               aes(y=n_lianas, x=log(tree_BA), color=zero))+
  scale_color_manual(values=c("black", "red"))+
  geom_smooth(data=filter(trees,n_lianas>0),
              aes(y=n_lianas, x=log(tree_BA),  group=local),
              method = "glm.nb", se=F)+
  geom_smooth(data=mutate(trees,zero=ifelse(n_lianas>0, 1, 0)),
              aes(y=zero, x=log(tree_BA),  group=local),
              method = "glm", se=F, method.args=list(family='binomial'), 
              color="blue", linetype=1)+
  facet_grid(geo~dist)+theme_bw()

## LIANAS  
# Rarefaction curves
df <- lianas %>% dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="BL"|site=="IL", 'East','South'),
         dist=ifelse(site=="BL"|site=="BS", 'Edge','Interior'))

ggplot(df, aes(x=x, y=y))+theme_bw()+
  geom_point(shape=16, data=df[which(df$method=="observed"),]) +
  geom_line(data=df[which(df$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_color_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,150,20))+
  ggtitle(" ")+facet_grid(dist~.)+
    xlab("Number of Individuals")+ ylab("Species Richness")+
    guides(color= F, shape = F, 
           fill=guide_legend(override.aes = 
                               list(alpha=.5, color='black'), 
                             title = ""))

# PCoA
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

# PERMANOVA
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




# https://chrischizinski.github.io/rstats/vegan-ggplot2/
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

# beta partiotining
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
  group_by(border) %>% summarise(sum(Freq))
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


## TREES
# Rarefaction curves
df <- trees %>%  mutate(border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border, Freq)  %>% dplyr::select(-Species) %>%
  iNEXT(datatype="abundance") %>% fortify(type=1) %>% 
  filter(method!="extrapolated") %>%
  mutate(geo=ifelse(site=="east_edge"|site=="east_interior", 'East','South'),
         dist=ifelse(site=="south_edge"|site=="east_edge", 'Edge','Interior'))

ggplot(df, aes(x=x, y=y))+theme_bw()+
  geom_point(shape=16, data=df[which(df$method=="observed"),]) +
  geom_line(data=df[which(df$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_color_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,900,100))+
  ggtitle(" ")+facet_grid(dist~.)+
  xlab("Number of Individuals")+ ylab("Species Richness")+
  guides(color= F, shape = F, 
         fill=guide_legend(override.aes = 
                             list(alpha=.5, color='black'), 
                           title = ""))

# PCoA
trees.t<-trees %>% dplyr::select(Species, Plot) %>%
  with(table(Species,Plot)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,Freq,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, Plot)) %>% t()
locality <- c(rep("EE",10), rep("ES",10),
              rep("IE",8), rep("IS",8))
edges <- c(rep("Edge",20), rep("Interior",16))
gps <-c(rep("East",10), rep("South",10),
        rep("East",8), rep("South",8))
trees_dist<-beta.pair(trees.t, index.family="jaccard")

bd_local<-betadisper(trees_dist[[3]],locality)
bd_edges<-betadisper(trees_dist[[3]],edges)
bd_gps<-betadisper(trees_dist[[3]],gps)
par(mfrow=c(1,3))
plot(bd_local); plot(bd_edges); plot(bd_gps)
anova(bd_local) # * ES
anova(bd_gps)
anova(bd_edges)

# PERMANOVA
trees_data <- trees.t %>% 
  as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Plot, Species)) %>% as.data.frame.matrix() 

(mod <- metaMDS(trees_data))
par(mfrow=c(1,1))
plot(mod, type="t")
plot(mod); ordihull(mod,groups=dados$local,draw="polygon",col="grey90",label=T)

adonis(trees_data ~ edges*gps, permutations = 999)
adonis(trees_data ~ locality, data=dados, permutations = 999)
envfit(mod, data.frame(locality,edges,gps))


# https://chrischizinski.github.io/rstats/vegan-ggplot2/
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

# beta partiotining - OK
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
  group_by(border) %>% summarise(sum(Freq))
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
