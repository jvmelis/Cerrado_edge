## Setting workspace and reading data
if(!require(readxl)){install.packages("readxl")}       # read_excel
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(iNEXT)){install.packages("iNEXT")}         # rarefaction curves
if(!require(cowplot)){install.packages("cowplot")}     # plot_grid

setwd("C:/Users/User/Documents/ARTIGOS - Meus/2 MS - Cerrado edge effect on lianas/2019 Austral Ecology")
lianas<-read_xlsx('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior'))
trees<-read_xlsx('cerrado_dados.xlsx', sheet = 'trees')

##### Rarefaction curves
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

## Plots
p_liana<-ggplot(df_lianas, aes(x=x, y=y, group=geo))+theme_bw()+
  geom_point(shape=16, 
             data=df_lianas[which(df_lianas$method=="observed"),]) +
  geom_line(data=df_lianas[which(df_lianas$method!="observed"),]) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = seq(0,150,20))+
  facet_grid(.~dist)+
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
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=geo, color=NULL), alpha=0.5)+
  scale_fill_grey(start = 0,end = 0.5)+
  scale_y_continuous(expand=c(0,0), breaks =seq(0,90,15))+
  scale_x_continuous(breaks = seq(0,900,150))+
  ggtitle(" ")+
  facet_grid(.~dist)+
  xlab(expression("Tree abundance (ind.plot"^-1*")"))+ 
  ylab("Tree species richness \n (Species per plot)")+
  theme(legend.position = 'none',
        panel.grid.minor = element_line(color = 'white'),
        panel.grid.major.x = element_line(color = 'white'),
        panel.grid.major.y = element_line(linetype = 2, colour = 'grey'))

### All together
cowplot::plot_grid(p_tree,p_liana,nrow = 2)
ggsave('Fig1_rarefactioncurves.jpg', width = 6, height = 8)
