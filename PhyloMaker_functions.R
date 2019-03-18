library (phytools)
library (RCurl)
library (ape)
# library (purrr) # to map the list of the sp.list function.



### Sourcing git https://github.com/jinyizju/S.PhyloMaker
### Paper: Qian, H., & Jin, Y. (2016). An updated megaphylogeny of plants, a tool for generating plant phylogenies and 
###        an analysis of phylogenetic community structure. Journal of Plant Ecology, 9(2), 233-239.



nodes <-read.csv(text=getURL("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/nodes"), sep ="\t", header=T)
example.splist <-read.csv(text=getURL("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/example.splist"), sep ="\t", header=T)
PhytoPhylo<-read.tree ("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/PhytoPhylo")
S.PhyloMaker <- getURL("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/R_codes%20for%20S.PhyloMaker", ssl.verifypeer = FALSE)
S.PhyloMaker
eval(parse(text = S.PhyloMaker))



### Function to generate a list of species to be use with PhyloMaker
sp.list<-function (Genus_species, Family){
  lista<-paste(Genus_species, Family, sep= "-")
  splista<-strsplit(lista, split = "-")
  species<-unlist(map(splista, 1))
  genus<-unlist(map (strsplit(species, split = "_"), 1))
  family<-unlist(map(splista,2))
  sp.list <-data.frame (species,genus, family)
}


######## Example

# example<-example.splist # read in the example species list.
# phylo<-PhytoPhylo # read in the megaphylogeny.
# nodes<-nodes # read in the nodes information of the megaphylogeny.
# result<-S.PhyloMaker(spList=example, tree=phylo, nodes=nodes) # run the function S.PhyloMaker.
# str(result) # the structure of the ouput of S.PhyloMaker.
# par(mfrow=c(1,3),mar=c(0,0,1,0)) # show the phylogenies of the three scenarios.
# plot(result$Scenario.1,cex=1.1,main="Scenarion One")
# 
