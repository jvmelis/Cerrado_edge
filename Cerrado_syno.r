if(!require(flora)){install.packages("flora")}

trees_table<-readxl::read_excel('cerrado_dados.xlsx', sheet = 'trees') %>% 
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

resulta[order(resulta$family,resulta$scientific.name),] %>%
  dplyr::select(family,scientific.name,ee,ei,se,si, total) %>% 
  write.csv("trees_species.csv")
trees<-readxl::read_excel('cerrado_dados.xlsx', sheet = 'trees') %>% 
  mutate(border = paste0(geo,"_",dist)) 
trees <- resulta %>% dplyr::select(Species, scientific.name) %>%
  right_join(trees) %>% dplyr::select(-Species) %>%
  mutate(Species=scientific.name,
         tree_BA=pi*((Diameter/2)^2))%>% dplyr::select(-scientific.name) 
rm(resulta,synom,trees_table)


lianas_table <- readxl::read_excel('cerrado_dados.xlsx', sheet = 'lianas') %>%
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
resulta<-resulta%>% dplyr::select(Species, scientific.name)
lianas<-readxl::read_excel('cerrado_dados.xlsx', sheet = 'lianas') %>%# filter(border=="BL"&Transect==4)
  mutate(liana_BA = pi*((Total_DBS/2)^2),
         geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior'))
lianas <- right_join(resulta, lianas) %>% dplyr::select(-Species) %>%
  mutate(Species=scientific.name)%>% dplyr::select(-scientific.name) 

rm(resulta,synom,lianas_table)
