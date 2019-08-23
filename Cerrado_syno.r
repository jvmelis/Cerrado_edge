if(!require(readxl)){install.packages("readxl")}       # read_excel
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(flora)){install.packages("flora")}

setwd("C:/Users/User/Documents/ARTIGOS - Meus/2 MS - Cerrado edge effect on lianas/2019 Austral Ecology")
trees_table<-read_xlsx('cerrado_dados.xlsx', sheet = 'trees') %>% 
  mutate(border = paste0(geo,"_",dist)) %>%
  dplyr::select(Species, border)  %>%
  with(table(Species, border)) %>% as.data.frame() %>%
  spread(border,Freq)
colnames(trees_table)<-c('Species','ee','ei','se','si')

synom<-flora::get.taxa(trees_table$Species)
synom$Species<-synom$original.search
resulta<-inner_join(trees_table,synom, by="Species") %>%
  dplyr::select(family,Species, scientific.name,taxon.status,
                notes, ee,ei,se,si, taxon.rank) %>%
  mutate(total=ee+ei+se+si,
         scientific.name = ifelse(notes=="not found", 
                                  Species,scientific.name),
         family = ifelse(!is.na(family),
                        family, 
                        sub(" sp.","", Species)),
         scientific.name = ifelse(taxon.rank=='genus',
                                  Species,
                                  scientific.name),
         scientific.name = ifelse(!is.na(scientific.name),
                                  scientific.name,
                                  Species))

resulta$family[resulta$family=='Lippia balansae']<-'Asteraceae' # odd
resulta[order(resulta$family,resulta$scientific.name),] %>%
  dplyr::select(family,scientific.name,ee,ei,se,si, total) %>% 
  write.csv("trees_species.csv")
trees<-read_xlsx('cerrado_dados.xlsx', sheet = 'trees') %>% 
  mutate(border = paste0(geo,"_",dist)) 
trees<-resulta %>% dplyr::select(Species, scientific.name) %>%
  right_join(trees,by="Species") %>% dplyr::select(-Species) %>%
  mutate(Species=scientific.name)%>% dplyr::select(-scientific.name) 
rm(resulta,synom,trees_table)


lianas_table <- read_xlsx('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'east', 'south'),
         dist = ifelse(border=='BL'|border=="BS", 'edge', 'interior'),
         border = paste0(geo,"_",dist)) %>% 
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>%
  spread(border,Freq)
colnames(lianas_table)<-c('Species','ee','ei','se','si')
synom<-flora::get.taxa(lianas_table$Species)
synom$Species<-synom$original.search
resulta<-inner_join(lianas_table,synom, by="Species") %>%
  dplyr::select(family,Species, scientific.name,taxon.status,
                notes, ee,ei,se,si, taxon.rank) %>%
  mutate(total=ee+ei+se+si,
         scientific.name = ifelse(notes=="not found", 
                                  Species,scientific.name),
         family = ifelse(!is.na(family),
                         family, 
                         sub(" sp.","", Species)),
         scientific.name = ifelse(taxon.rank=='genus',
                                  Species,
                                  scientific.name),
         scientific.name = ifelse(!is.na(scientific.name),
                                  scientific.name,
                                  Species))
resulta[order(resulta$family,resulta$scientific.name),] %>%
  dplyr::select(family,scientific.name,ee,ei,se,si, total) %>% 
  write.csv("liana_species.csv")

lianas<-read_xlsx('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior'))

lianas<-resulta %>% dplyr::select(Species, scientific.name) %>%
  right_join(lianas, by="Species") %>% dplyr::select(-Species) %>%
  mutate(Species=scientific.name)%>% dplyr::select(-scientific.name) 

rm(resulta,synom,lianas_table)
