# Manuscript "Contrasting edge influence on lianas and trees in a cerrado savanna remnant"
Authors: 
Juliano van Melis

Maria Gabriela Gutierrez Camargo

Paula Guimarães Carvalho

Leonor Patricia Cerdeira Morellato

Maria Tereza Grombone-Guaratini

Submited to: **Austral Ecology**

Submission date:

This repository has the following files: `clean_data.r`,`clean_data.csv`,`Cerrado_edge_2019.r`,  and `cerrado_dados.xlsx.r`.

### clean_data.r

Based on `cerrado_dados.xlsx`, we edited (clean) the data in this script.

### clean_data.csv

The result of `clean_data.r`. It has the following columns:

`Plot`: plot id (25 m2).	

`trees_BA`: sum of tree stand basal area per plot.	

`trees_ab`: stem abundance of trees per plot.	

`lianas_BA`: sum of liana stand basal area per plot.	

`lianas_ab`: stem abundance of lianas per plot.		

`geo`: If the plot is located in the `south` or `east`.	

`dist`: If the plot is located in the `edge` or `interior`.		

`Tree_cov`: Tree cover, calculated by hemispheric photos.	

`Regen`: % of soil covered by regerative stratum.	

`Palm_cov`% of soil covered by Palms stratum.		

`Bare_soil`: % of soil covered by bare soil stratum.	

`native`: % of soil covered by native grass stratum.	

`exotic`: % of soil covered by exotic grass stratum.	

`brom`:% of soil covered by bromeliad stratum.		

`PAR`: Photosynthetically active radiation.

`SOM`: Soil Organic Matter.	

`Al`: Alumnium concentration in the soil.	

`Mn`: Manganese concentration in the soil.

`local`: Concatenate `geo` and `dist` variables.

### Cerrado_edge_2019.r
All the analyses performed and graphs constructed in a R script.

## cerrado_dados.xlsx
This spreadsheet has five sheets:

### sheet: plots
Data from sampled plots:

`Plot`: plot id (25 m2).

`geo`: If the plot is located in the `south` or `east`.	

`dist`: If the plot is located in the `edge` or `interior`.

`Tree_cov`: % of soil covered by Trees (in %, based on mean of observed % for each plot)

`Regen`: % of soil covered by regerative stratum	(in %, based on mean of observed % for each plot).

`Palm_cov`% of soil covered by Palms stratum	(in %, based on mean of observed % for each plot).

`Bare_soil`: % of soil covered by bare soil stratum		(in %, based on mean of observed % for each plot).

`native`: % of soil covered by native grass stratum		(in %, based on mean of observed % for each plot).

`exotic`: % of soil covered by exotic grass stratum		(in %, based on mean of observed % for each plot).

`brom`:% of soil covered by bromeliad stratum		(in %, based on mean of observed % for each plot).

`canopy`: Tree cover, obtained by hemispheric photos (in %).

`PAR`: Photosynthetically active radiation (in ).

`SOM`: Soil Organic Matter	(in mg).

`Al`: Alumnium concentration in the soil (in mg).	

`Mn`: Manganese concentration in the soil (in mg).

### sheet: lianas
`Date`: Sample date (DD/MM/YYYY). 
`border`: four factor variable regarding plot cardinal orientation and edge influence, where `BL`: east edge; `BS`: south edge; `IS`: south interior; and  `IL`: east interior

`Transect`: Transect ID	(numbered).

`Plot`: Plot ID (`border`+`Transect`).	

`id_tree`: First Host tree ID	(same as `ind` in `trees` sheet).

`DBH_tree`: Host #1 tree Diameter at breat height	(in cm).

`id_tree_1`: Second host tree ID.	

`DBH_tree_1`: Host #2 tree Diameter at breat height (in cm).	

`id_tree_2`: Third host tree ID.		

`DBH_tree_2`: Host #3 tree Diameter at breat height (in cm).	

`id`: Liana ID.	

`DBS_30`: Liana Diameter at 30 cm above soil (in cm).	

`DBS_a`: 	Second stem DBS from the same visually individual (in cm).

`DBS_b`: Third stem DBS from the same visually individual (in cm).	

`DBS_c`: Fourth stem DBS from the same visually individual (in cm).	

`Total_DBS`: Sum of all DBS (in cm).

`DBH_130`: Diameter at 130 cm long from the soil (according Gerwing et al. 2006) (in cm).	

`DBH_a`: Second stem DBH from the same visually individual	(in cm).

`DBH_b`: Third stem DBH from the same visually individual	(in cm).	

`DBH_c`: Fourth stem DBH from the same visually individual	(in cm).	

`Family`: Liana botanical family.	

`Species`: Liana species.	

`Morphospecies`: Codified Liana species.

`Harvested`: where the liana	were harvested: Two factors variable with `NA` (`Alta`: harvested at the canopy; `ok`: deposited in the herbarium, `NA`: identified at field).

`obs`: Field observations.	

### sheet: trees
`Plot`: Plot ID (same as in `lianas` sheet)

`ind`: Tree ID	

`Family`: Botanical Family	

`Species`: Tree species 	

`Girth`: Tree girth at breast height	(in cm)

`Diameter`: Diameter at breat height based on `Girth` ($Girth/\pi$)	(in cm)

`Height`: Tree height (in meters)	

`WD`: Wood Density (based on literature)	

`geo`: If the plot is located in the `south` or `east`.	

`dist`: If the plot is located in the `edge` or `interior`.

### host_tree
`id_1`: First Host tree ID	(same as `ind` in `trees` sheet).	

`n_1`: Number of lianas on host tree #1	

`id_2`: Second Host tree ID	(same as `ind` in `trees` sheet).	

`n_2`: Number of lianas on host tree #2	

`id_3`: Third Host tree ID	(same as `ind` in `trees` sheet).	

`n_3`: Number of lianas on host tree #3


### Field_guide
`Previously`: Field identification

`Now`: Data homogeneization.
