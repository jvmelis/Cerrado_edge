setwd("C:/Users/User/Documents/ARTIGOS - Meus/2 MS - Cerrado edge effect on lianas/2019 Austral Ecology")
parcelas<-read_excel("cerrado_dados.xlsx",sheet="plots",na="") %>%
  dplyr::select(-c("AGB","LAGB","AGB_total","Tree_abun",
                   "Liana_Species","Liana_abun","canopy","A","B","C"))

lianas<-read_excel("cerrado_dados.xlsx",sheet="lianas",na="")
lianas<-lianas %>% mutate(BasalArea = (((Total_DBS/2)^2)*pi)*0.0001 )
parcelas<-lianas %>% group_by(Plot)%>% 
  summarise(lianas_BA=sum(BasalArea),
            lianas_ab=n()) %>% ungroup() %>% inner_join(parcelas) 
trees<-read_excel("cerrado_dados.xlsx",sheet="trees",na="")%>%
  dplyr::select(Plot, Diameter, Family, Species)
trees<-trees %>% mutate(Tree_BasalArea = (((Diameter/2)^2)*pi)*0.0001 )

parcelas<-trees %>% group_by(Plot)%>% 
  summarise(trees_BA=sum(Tree_BasalArea),
            trees_ab=n())%>% ungroup() %>%inner_join(parcelas) 
parcelas<-mutate(parcelas,local=paste(geo,dist,sep="_"))

write.csv(parcelas, "clean_data.csv")
