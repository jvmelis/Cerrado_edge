if(!require(readxl)){install.packages("readxl")}       # read_excel
if(!require(tidyverse)){install.packages('tidyverse')} # data manipulation & graphs
if(!require(cowplot)){install.packages('cowplot')}     # plot_grid
if(!require(VennDiagram)){install.packages('VennDiagram')}

# source: "https://raw.githubusercontent.com/cran/VennDiagram/master/R/draw.quad.venn.R"
draw.quad.venn <- function(
  area1, area2,area3, area4,
  n12, n13, n14, n23, n24, n34,
  n123, n124, n134, n234, n1234,
  category = rep("", 4),
  lwd = rep(2, 4),
  lty = rep("solid", 4),
  col = rep("black", 4),
  fill = NULL, alpha = rep(0.5, 4),
  label.col = rep("black", 15),
  cex = rep(1, 15),
  fontface = rep("plain", 15),
  fontfamily = rep("serif", 15),
  cat.pos = c(-15, 15, 0, 0),
  cat.dist = c(0.22, 0.22, 0.11, 0.11),
  cat.col = rep("black", 4),
  cat.cex = rep(1, 4),
  cat.fontface = rep("plain", 4),
  cat.fontfamily = rep("serif", 4),
  cat.just = rep(list(c(0.5, 0.5)), 4),
  rotation.degree = 0, rotation.centre = c(0.5, 0.5),
  ind = TRUE, cex.prop=NULL, print.mode = "raw", sigdigs=3,
  direct.area = FALSE, area.vector = 0,
  ...
) {
  cat.pos <- cat.pos + rotation.degree;
  if(direct.area){
    areas <- area.vector;
    #create the variables and assign their values from the area vector
    for(i in 1:15)
    {
      assign(paste("a",i,sep=""),area.vector[i]);
    }
  }
  else {
    # generate partial areas from given arguments
    a6  <- n1234;
    a12 <- n123 - a6;
    a11 <- n124 - a6;
    a5  <- n134 - a6;
    a7  <- n234 - a6;
    a15 <- n12 - a6 - a11 - a12;
    a4  <- n13 - a6 - a5 - a12;
    a10 <- n14 - a6 - a5 - a11;
    a13 <- n23 - a6 - a7 - a12;
    a8  <- n24 - a6 - a7 - a11;
    a2  <- n34 - a6 - a5 - a7;
    a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15;
    a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15;
    a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13;
    a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11;
    
    # check plausibility and 0 partial areas
    areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15);
  }
  areas.error <- c(
    "a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13",
    "a2  <- n34 - a6 - a5 - a7",
    "a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11",
    "a4  <- n13 - a6 - a5 - a12",
    "a5  <- n134 - a6",
    "a6  <- n1234",
    "a7  <- n234 - a6",
    "a8  <- n24 - a6 - a7 - a11",
    "a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15",
    "a10 <- n14 - a6 - a5 - a11",
    "a11 <- n124 - a6",
    "a12 <- n123 - a6",
    "a15 <- n12 - a6 - a11 - a12",
    "a13 <- n23 - a6 - a7 - a12",
    "a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15"
  );
  for (i in 1:length(areas)) {
    if (areas[i] < 0) {
      flog.error(paste("Impossible:", areas.error[i], "produces negative area"),name="VennDiagramLogger")
      stop(paste("Impossible:", areas.error[i], "produces negative area"));
    }
  }
  
  ## rescaling area labels to be proportional to area
  if(length(cex.prop) > 0){
    
    if(length(cex.prop) != 1) {
      flog.error("Value passed to cex.prop is not length 1",name="VennDiagramLogger")
      stop("Value passed to cex.prop is not length 1")
    }
    ## figure out what function to use
    func = cex.prop
    if(class(cex.prop) != "function"){
      if(cex.prop == "lin"){
        func = function(x) x
      }
      else if(cex.prop == "log10"){
        func = log10
      }
      else flog.error(paste0("Unknown value passed to cex.prop: ", cex.prop),name="VennDiagramLogger")
      stop(paste0("Unknown value passed to cex.prop: ", cex.prop))
    }
    ## rescale areas
    maxArea = max(areas)            
    for(i in 1:length(areas)){                
      cex[i] = cex[i] * func(areas[i]) / func(maxArea)
      if(cex[i] <= 0) stop(paste0("Error in rescaling of area labels: the label of area ",
                                  i, " is less than or equal to zero"))
    }
  }
  # initialize gList to hold all Grobs generated
  grob.list <- grid::gList();
  # plot the ellipses of the Venn diagram
  ellipse.positions <- matrix(
    nrow = 4,
    ncol = 7
  );
  colnames(ellipse.positions) <- c('x', 'y', 'a', 'b', 'rotation', 'fill.mapping', 'line.mapping');
  ellipse.positions[1,] <- c(0.65, 0.47, 0.35, 0.20,  45, 2, 2);
  ellipse.positions[2,] <- c(0.35, 0.47, 0.35, 0.20, 135, 1, 1);
  ellipse.positions[3,] <- c(0.50, 0.57, 0.33, 0.15,  45, 4, 4);
  ellipse.positions[4,] <- c(0.50, 0.57, 0.35, 0.15, 135, 3, 3);
  # draw the ellipses themselves
  for (i in 1:4) {
    grob.list <- grid::gList(
      grob.list,
      VennDiagram::ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = grid::gpar(
          lty = 0,
          fill = fill[ellipse.positions[i,'fill.mapping']],
          alpha = alpha[ellipse.positions[i,'fill.mapping']]
        )
      )
    );
  }
  # draw the ellipse borders
  for (i in 1:4) {
    grob.list <- grid::gList(
      grob.list,
      ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = gpar(
          lwd = lwd[ellipse.positions[i,'line.mapping']],
          lty = lty[ellipse.positions[i,'line.mapping']],
          col = col[ellipse.positions[i,'line.mapping']],
          fill = 'transparent'
        )
      )
    );
  }
  
  # create the labels
  label.matrix <- matrix(
    nrow = 15,
    ncol = 3
  );
  colnames(label.matrix) <- c('label', 'x', 'y');
  
  label.matrix[ 1,] <- c(a1,  0.350, 0.77);
  label.matrix[ 2,] <- c(a2,  0.500, 0.69);
  label.matrix[ 3,] <- c(a3,  0.650, 0.77);
  label.matrix[ 4,] <- c(a4,  0.310, 0.67);
  label.matrix[ 5,] <- c(a5,  0.400, 0.58);
  label.matrix[ 6,] <- c(a6,  0.500, 0.47);
  label.matrix[ 7,] <- c(a7,  0.600, 0.58);
  label.matrix[ 8,] <- c(a8,  0.690, 0.67);
  label.matrix[ 9,] <- c(a9,  0.180, 0.58);
  label.matrix[10,] <- c(a10, 0.320, 0.42);
  label.matrix[11,] <- c(a11, 0.425, 0.38);
  label.matrix[12,] <- c(a12, 0.575, 0.38);
  label.matrix[13,] <- c(a13, 0.680, 0.42);
  label.matrix[14,] <- c(a14, 0.820, 0.58);
  label.matrix[15,] <- c(a15, 0.500, 0.28);
  
  processedLabels <- rep("",length(label.matrix[,'label']));
  if(print.mode[1] == "percent"){
    processedLabels <- paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),"%",sep="");
    if(isTRUE(print.mode[2] == "raw"))
    {
      processedLabels <- paste(processedLabels,"\n(",label.matrix[,'label'],")",sep="");
    }
  }
  if(print.mode[1] == "raw"){
    processedLabels <- label.matrix[,'label'];
    if(isTRUE(print.mode[2] == "percent"))
    {
      processedLabels <- paste(processedLabels,"\n(",paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),"%)",sep=""),sep="");
    }
  }
  
  
  for (i in 1:nrow(label.matrix)) {
    grob.list <- grid::gList(
      grob.list,
      textGrob(
        label = processedLabels[i],
        x = label.matrix[i,'x'],
        y = label.matrix[i,'y'],
        gp = gpar(
          col = label.col[i],
          cex = cex[i],
          fontface = fontface[i],
          fontfamily = fontfamily[i]
        )
      )
    );
  }
  
  # find the location and plot all the category names
  cat.pos.x <- c(0.18, 0.82, 0.35, 0.65);
  cat.pos.y <- c(0.58, 0.58, 0.77, 0.77);
  
  for (i in 1:4) {
    
    # work out location of the category label
    this.cat.pos <- find.cat.pos(
      x = cat.pos.x[i],
      y = cat.pos.y[i],
      pos = cat.pos[i],
      dist = cat.dist[i]
    );
    
    # then print it
    grob.list <- gList(
      grob.list,
      textGrob(
        label = category[i],
        x = this.cat.pos$x,
        y = this.cat.pos$y,
        just = cat.just[[i]],
        gp = gpar(
          col = cat.col[i],
          cex = cat.cex[i],
          fontface = cat.fontface[i],
          fontfamily = cat.fontfamily[i]
        )
      )
    );
  }
  
  # adjust grob.list to fit and return grob.list
  grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
  
  # draw diagram before returning gList is specified by user
  if (ind) { grid.draw(grob.list); }
  return(grob.list);
}



# Lianas
lianas_venn<-readxl::read_excel('cerrado_dados.xlsx', sheet = 'lianas') %>%
  mutate(geo = ifelse(border=='BL'|border=="IL", 'East', 'South'),
         dist = ifelse(border=='BL'|border=="BS", 'Edge', 'Interior')) %>%
  dplyr::select(Species, border) %>%
  with(table(Species,border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>%
  as.data.frame()%>% spread(border,Freq)

dados_venn<-data.frame(summarise(filter(lianas_venn,BL!=0),area1=n()),
                  summarise(filter(lianas_venn,BS!=0),area2=n()),
                  summarise(filter(lianas_venn,IL!=0),area3=n()),
                  summarise(filter(lianas_venn,IS!=0),area4=n()),
                  summarise(filter(lianas_venn,BL!=0&BS!=0),n12=n()),
                  summarise(filter(lianas_venn,BL!=0&IL!=0),n13=n()),
                  summarise(filter(lianas_venn,BL!=0&IS!=0),n14=n()),
                  summarise(filter(lianas_venn,BS!=0&IL!=0),n23=n()),
                  summarise(filter(lianas_venn,BS!=0&IS!=0),n24=n()),
                  summarise(filter(lianas_venn,IL!=0&IS!=0),n34=n()), 
                  summarise(filter(lianas_venn,BL!=0&BS!=0&IL!=0),n123=n()),
                  summarise(filter(lianas_venn,BL!=0&BS!=0&IS!=0),n124=n()), 
                  summarise(filter(lianas_venn,BL!=0&IL!=0&IS!=0),n134=n()),
                  summarise(filter(lianas_venn,BS!=0&IL!=0&IS!=0),n234=n()), 
                  summarise(filter(lianas_venn,BL!=0&BS!=0&IL!=0&IS!=0), n1234=n())) 

lianas_venn_plot <- draw.quad.venn(
  area1 = dados_venn$area1, area2 = dados_venn$area2, 
  area3 = dados_venn$area3, area4 = dados_venn$area4,
  n12 = dados_venn$n12, n13 = dados_venn$n13, n14 = dados_venn$n14,
  n23 = dados_venn$n23, n24 = dados_venn$n24, n34 = dados_venn$n34,
  n123 = dados_venn$n123, n124 = dados_venn$n124, n134 = dados_venn$n134,
  n234 = dados_venn$n234, n1234 = dados_venn$n1234,
  category = c("East edge", "South edge", 
               "East interior", "South interior"),
  fill = c("white", "white", "grey", "grey")
)

# Trees
source("Cerrado_syno.r")
trees_venn<-trees %>%
  dplyr::select(Species, border)  %>%
  with(table(Species, border)) %>% as.data.frame() %>% 
  mutate(Freq=ifelse(Freq>0,1,NA)) %>% 
  filter(!is.na(Freq)) %>%
  with(table(Species, border)) %>%as.data.frame() %>% 
  spread(border,Freq)

dados_venn<-cbind(summarise(filter(trees_venn,east_edge!=0),area1=n()),
             summarise(filter(trees_venn,south_edge!=0),area2=n()),
             summarise(filter(trees_venn,east_interior!=0),area3=n()),
             summarise(filter(trees_venn,south_interior!=0),area4=n()),
             summarise(filter(trees_venn,east_edge!=0&south_edge!=0),n12=n()),
             summarise(filter(trees_venn,east_edge!=0&east_interior!=0),n13=n()),
             summarise(filter(trees_venn,east_edge!=0&east_interior!=0),n14=n()),
             summarise(filter(trees_venn,south_edge!=0&east_interior!=0),n23=n()),
             summarise(filter(trees_venn,south_edge!=0&east_interior!=0),n24=n()),
             summarise(filter(trees_venn,east_interior!=0&south_interior!=0),n34=n()),
             summarise(filter(trees_venn,east_edge!=0&south_edge!=0&east_interior!=0),n123=n()),
             summarise(filter(trees_venn,east_edge!=0&south_edge!=0&south_interior!=0),n124=n()),
             summarise(filter(trees_venn,east_edge!=0&east_interior!=0&south_interior!=0),n134=n()),
             summarise(filter(trees_venn,south_edge!=0&east_interior!=0&south_interior!=0),n234=n()),
             summarise(filter(trees_venn, east_edge!=0&south_edge!=0&east_interior!=0&south_interior!=0),n1234=n()),
             summarise(filter(trees_venn,east_edge!=0),area1=n()))

trees_venn_plot <- draw.quad.venn(
  area1 = dados_venn$area1, area2 = dados_venn$area2, 
  area3 = dados_venn$area3, area4 = dados_venn$area4,
  n12 = dados_venn$n12, n13 = dados_venn$n13, n14 = dados_venn$n14,
  n23 = dados_venn$n23, n24 = dados_venn$n24, n34 = dados_venn$n34,
  n123 = dados_venn$n123, n124 = dados_venn$n124, n134 = dados_venn$n134,
  n234 = dados_venn$n234, n1234 = dados_venn$n1234,
  category = c("East edge", "South edge", 
               "East interior", "South interior"),
  fill = c("white", "white", "grey", "grey")
  )
dev.off()

Fig3<-cowplot::plot_grid(lianas_venn_plot,trees_venn_plot, 
                   nrow = 2,
                   labels = c("a) Lianas","b) Trees"))
rm(list=c(dados,lianas_venn, lianas_venn_plot, trees_venn, trees_venn_plot, dados_venn))
