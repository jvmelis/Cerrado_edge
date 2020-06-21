parcelas<-readxl::read_excel("cerrado_dados.xlsx",sheet="plots",na="")

lianas<-readxl::read_excel("cerrado_dados.xlsx",sheet="lianas",na="")
lianas<-lianas %>% mutate(BasalArea = (((Total_DBS/2)^2)*pi)*0.0001 )
parcelas<-lianas %>% group_by(Plot)%>% 
  summarise(lianas_BA=sum(BasalArea),
            lianas_ab=n()) %>% ungroup() %>% inner_join(parcelas) 
trees<-readxl::read_excel("cerrado_dados.xlsx",sheet="trees",na="")%>%
  dplyr::select(Plot, Diameter, Family, Species)
trees<-trees %>% mutate(Tree_BasalArea = (((Diameter/2)^2)*pi)*0.0001 )

parcelas<-trees %>% group_by(Plot)%>% 
  summarise(trees_BA=sum(Tree_BasalArea),
            trees_ab=n())%>% ungroup() %>%inner_join(parcelas) 
parcelas<-mutate(parcelas,local=paste(geo,dist,sep="_"))

write.csv(parcelas, "clean_data.csv")
