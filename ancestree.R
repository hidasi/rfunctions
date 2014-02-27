
# for much more info: http://rfunctions.blogspot.com.br/2013/04/ancestree-function-get-phylogeny-that.html

#################################################

library(ape)
library(picante)


ancestree<-function(sample,tree){

	cnames<-colnames(sample) # get sample species' names

	ca<-mrca(tree) # calculate most common ancestor for all pairs of species in the tree
	ca<-data.frame(ca) # transform into a data frame (just to avoid some errors)

	nameca<-ca[cnames,cnames] # get ancestors only for sample's species

	allanc<-unique(unlist(c(nameca[,1:length(nameca)]))) # all possible common 'ancestor nodes'
	unicnodes<-sort(allanc[allanc>length(tree$tip.label)]) # due to 'mrca' function, avoid ancestors of only one species (e.g. ancestor node of two 'sp1')

	nnodefor<-numeric(length(unicnodes)) # for 'for' =] this will be the number of descendant species for each node
	nprunsizefor<-numeric(length(unicnodes)) # for 'for' =] this will be the number of sample's species within the descendant species of each 'ancestor node'

	for(i in 1:length(unicnodes)){
		phy<-extract.clade(tree,unicnodes[i]) # get each phylogeny starting at ancestor nodes
		nnodefor[i]<-length(phy$tip.label) # get number of descendant species for each node
		prunt<-prune.sample(sample,phy) # prune tree using sample's species
		nprunsizefor[i]<-length(prunt$tip.label) # get number of sample's species in the phylogeny resulted from each ancestor node
	}

	resul<-rbind(unicnodes,nprunsizefor,nnodefor) # combine all results
	resul<-resul[,order(-resul[2,],resul[3,])]  # sort data first by highest to lowest nprunsize, then lowest to highest nnodefor
	bestnode<-resul[1,1] # this is the 'ancestor node' with all sample's species and lowest number of terminal nodes (species); most common ancestor
	bestphylo<-extract.clade(tree,bestnode) # get the phylogeny starting at the best ancestor node

	return(bestphylo)
}