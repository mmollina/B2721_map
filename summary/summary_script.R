#####
## Analitical procedures to construct the B2721 potato map - Map Summary
## -----------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
#####
require(mappoly)
load("~/repos/B2721_map/mapping/B2721_analysis.rda")
load("~/repos/B2721_map/meiosis/meiosis_results.rda")
summary_maps(map.object = final.maps)
