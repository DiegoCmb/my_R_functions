library (MCMCglmm)
library(phangorn)
library (dplyr)

weighting_and_subsetting <-
  function (database, weight_approach = 'inverse_sample_size') {
    # Estimate the weights.
    # remove NAs from sample sizes col and estimates weights 
    # to be used in MCMCglmm
    database <-
      subset(database, database$sample_size != 'NA') # because we want to weight by sample sizes
    # we will have to remove obs with out this info
    if (weight_approach == "inverse_sample_size") {
      database$mev <- 1 / (database$sample_size)
      database$weight <- 1 / database$mev
      database$mesd <- sqrt(database$mev)
      return(database)
    } else {
      return(database)
    }
  }

force.ultrametric <- function(tree, method = c("nnls", "extend")) {
  method <- method[1]
  if (method == "nnls")
    tree <- nnls.tree(cophenetic(tree), tree,
                      rooted = TRUE, trace = 0)
  else if (method == "extend") {
    h <- diag(vcv(tree))
    d <- max(h) - h
    ii <- sapply(1:Ntip(tree), function(x, y)
      which(y == x),
      y = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + d
  } else
    cat("method not recognized: returning input tree\n\n")
  tree
}

mcmcglmm_dignostics <-
  function (modelo_fix,
            modelo_random,
            database,
            mev_col = NULL,
            gelman_diag = "True",
            nitt_x = 13000,
            thin_x = 10 ,
            burnin_x = 3000,
            prior_x,
            save_output = "False",
            saving_path,
            file_name) {
    # Function to run one or several models
    # gelman_diag == False will run only one model that will return the standar output of MCMCglmm
    # gelman_diag default (== True) will run the gelman diagnostic
    # nitt: no. of iterations; 13000 iterations are the original default
    # thin: no. of iterations that must pass after saving one iteration; default thin = 10, so 10 iterations will pass
    # before saving 1 iteration.
    # burnin: period of iterations before start saving iterations using the algorith thin. burnin default = 3000
    # path_name_saving_file: for saving output as .RData Include the path or the file name using the extension .RData
    print (mev_col)
    # Running just one model
    if (gelman_diag == "False") {
      m0.1 <-
        MCMCglmm (
          as.formula(modelo_fix),
          random = as.formula (modelo_random),
          family = "gaussian",
          data = database,
          mev = mev_col ,
          verbose = T,
          nitt = nitt_x,
          thin = thin_x,
          burnin = burnin_x,
          pr = T,
          prior = prior_x
        )
      # to save data in .RData
      if (save_output == "True") {
        path_file_to_save <- file.path(saving_path, file_name)
        save(m0.1, file = path_file_to_save)
        return (m0.1)
      } else {
        return (m0.1)
      }
      
      # Running multiple models for gelman diagnostic
    } else {
      m0.1 <-
        MCMCglmm (
          as.formula(modelo_fix),
          random = as.formula (modelo_random),
          family = "gaussian",
          data = database,
          mev = mev_col,
          verbose = T,
          nitt = nitt_x,
          thin = thin_x,
          burnin = burnin_x,
          pr = T,
          prior = prior_x
        )
      
      m0.2 <-
        MCMCglmm (
          as.formula(modelo_fix),
          random = as.formula (modelo_random),
          family = "gaussian",
          data = database,
          mev = mev_col ,
          verbose = T,
          nitt = nitt_x,
          thin = thin_x,
          burnin = burnin_x,
          pr = T,
          prior = prior_x
        )
      
      m0.3 <-
        MCMCglmm (
          as.formula(modelo_fix),
          random = as.formula (modelo_random),
          family = "gaussian",
          data = database,
          mev = mev_col,
          verbose = T,
          nitt = nitt_x,
          thin = thin_x,
          burnin = burnin_x,
          pr = T,
          prior = prior_x
        )
      
      gel.diag_SOL <-
        gelman.diag(list(m0.1$Sol[, 1], m0.2$Sol[, 1], m0.3$Sol[, 1]))
      gel.diag_VCV <-
        gelman.diag(list(m0.1$VCV[, c(1, 3)], m0.2$VCV[, c(1, 3)], m0.3$VCV[, c(1, 3)]))
      gel.diag_DEV <-
        gelman.diag(list(m0.1$Deviance, m0.2$Deviance, m0.3$Deviance))
      
      m0.1_DIC <- m0.1$DIC
      m0.2_DIC <- m0.2$DIC
      m0.3_DIC <- m0.3$DIC
      
      list_of_DIC <- list (m0.1_DIC, m0.2_DIC, m0.3_DIC)
      names(list_of_DIC) <-
        c("DIC_model_1", "DIC_model_2", "DIC_model_3")
      
      model_list_output <-
        list (
          "/////////////1) First_model  /////////////",
          m0.1,
          # the number within parethesis indicates the parsing indexation in the list
          "/////////////3) Second_model /////////////",
          m0.2,
          "/////////////4) Third_model  /////////////",
          m0.3,
          "/////////////6) Gelman convergence diagnostic SOL /////////////",
          gel.diag_SOL,
          "/////////////8) Gelman convergence diagnostic SOL /////////////",
          gel.diag_VCV,
          "/////////////10) Gelman convergence diagnostic DEV /////////////",
          gel.diag_DEV,
          "/////////////12) DIC for each model /////////////",
          list_of_DIC
        )
      
      if (save_output == "True") {
        path_file_to_save <- file.path(saving_path, file_name)
        save(model_list_output, file = path_file_to_save)
        return (model_list_output)
      } else{
        return (model_list_output)
      }
    }
  }


database_adjments_for_mcmc <- function(my_database, animal, my_tree) {
  # This function helps to create database to be run with MCMCglmm
  # Returns a list with 4 objects
  # 1. Database with column for animal, and estimations of sampling error
  # such as mev. Obs with NAs in sample size were removed.
  # 2. A pruned ultrametric phylogenetic tree with no nodes names and singletones rm
  # 3. Inverse matrix IA
  # 4. a vector including species names that were pruned from tree since lack of match with database
  # about database
  # needs other functions such as weighting_and_subsetting
  # library(RCurl)
  # mcmcglmm_fun <- getURL(paste("https://raw.githubusercontent.com/DiegoCmb/",
  #   "my_R_functions/master/mcmcglmm_extra_functions.R", sep = ""),
  #    ssl.verifypeer = FALSE)
  # eval(parse(text = mcmcglmm_fun))
  
  #-------------- Woring database -----------------------------------------------
  ####
  # If the sampling error around the true value is approximately normal, and the
  # variance of the sampling errors known, then random effect meta-analyses can be
  # fitted by passing the sampling variances to the mev argument of MCMCglmm.
  # (https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf)
  # for doing this I use the weight_approach function using the
  # inverse_sample_sizes option.
  # about sample sizes in He
  # https://www.jstor.org/stable/pdf/30244434.pdf?refreqid=excelsior%3A5833a1ce288f061c96631992d5153a9f
  #### extra information
  # mev.all<-1/d.sex.all$sample_size # sample sizes can be used based on Schmidt and Hunter (1977)
  # or 1/sample sizes based on https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
  # Leimu et al Journal of eCology 2006 How general are positive realtionships between plant
  # population size, fintess and genetic variation? weight by the number of sampled populations
  #mev.all<-1/(d.sex.all$sample_size-3)		# <---------- se usa en el analisis: sampling variance of the effect sizes
  #weight.all<-1/mev.all
  #mesd.all<-sqrt(mev.all)
  # more about weigthing that might be useful
  # https://www.r-bloggers.com/sampling-weights-and-multilevel-modeling-in-r/
  
  
  my_database <-
    weighting_and_subsetting(my_database, # gets the sampling error and rm NAs
                             weight_approach = "inverse_sample_size")
  colnames(my_database)[colnames(my_database) == animal] <- "animal"
  
  #-------------- Working phylogenetic tree
  phylo.tree <- force.ultrametric(tree_scenario3) ## default method
  is.ultrametric(phylo.tree)
  is.ultrametric(tree_scenario3)
  dif_sp <-
    dplyr::setdiff(phylo.tree$tip.label, my_database$animal) # species to be pruned
  dplyr::setdiff(phylo.tree$tip.label, my_database$animal) -> drop_this_sp.all
  phylo.tree.all <- drop.tip(phylo.tree, drop_this_sp.all)
  phylo.tree.all <-
    collapse.singles(phylo.tree.all) # to remove singletones
  phylo.tree.all$node.label <- NULL
  
  #-------------- Working Inverse Matrix for MCMCglmm procedure ----------------------
  ##### Getting the Invers matrix A
  # However, in most statistical applications it is not A
  # that is required, but its inverse A. For pedigrees this
  # matrix can be very large and efficient ways of obtaining
  # the inverse made the fitting of these models practical
  # (Henderson, 1976; Quaas, 1976;
  # Meuwissen & Luo, 1992 IN HADFIELD AND NAKAGAWA 2010)
  
  IA.all <- inverseA(phylo.tree.all, nodes = "TIPS")$Ainv
  
  output <-
    list(
      "database" = my_database,
      "tree" = phylo.tree.all,
      "IA" = IA.all,
      "pruned_species" = dif_sp
    )
  return (output)
}


# based on Scripts from Nakagawa and 
# https://github.com/daniel1noble/metaAidR/blob/master/R/I2.R
# estimating 
# mev<-mev.all
# W<-1/mev  # is this correct?????? this is the formula when is based on 
# # measurment of variance, but I dont know if with sample sizes
# # should follow the same logic. 
# s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()

# heterogeneity for focal factor without considering phylogeny
mm_I2.factor<-function (modelo_mm, focal_comp, mev){
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(modelo_mm$VCV[,focal_comp])/
    (modelo_mm$VCV[,focal_comp]+modelo_mm$VCV[,"units"]+s2)
  estimado2<-data.frame(posterior.mode(estimado))
  names (estimado2)<-'heterogenity (I2)'
  rownames(estimado2)<-paste('mm_I2', focal_comp, sep="." )
  return (estimado2)
}
# heterogeneity for the same focal factor but considering phylogeny
pm_I2.focal.phylo<-function (modelo_pm, focal_comp, mev){ 
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(modelo_pm$VCV[,focal_comp]+modelo_pm$VCV[,"animal"])/
    (modelo_pm$VCV[, focal_comp]+modelo_pm$VCV[,"animal"]+
       modelo_pm$VCV[,"units"]+s2)
  estimado2<-data.frame(posterior.mode(estimado))
  names(estimado2)<-"heterogenity (I2)"
  rownames (estimado2)<-paste("pm_I2",focal_comp, sep= ".")
  return (estimado2)
}
# heterogeneity due to phylogeny considering the factor
pm_I2.phylo<-function (modelo_pm, focal_comp, mev){
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(modelo_pm$VCV[,"animal"])/
    (modelo_pm$VCV[, focal_comp]+modelo_pm$VCV[,"animal"]+
       modelo_pm$VCV[,"units"]+s2)
  estimado2<-data.frame(posterior.mode(estimado))
  names(estimado2)<-"heterogenity (I2)"
  rownames (estimado2)<-paste("pm_I2.animal", focal_comp, sep= ".")
  return (estimado2)
}
# Phylogenetic heredability, similar to lambda pagel
pm_H2.phylo<-function (modelo_pm, focal_comp, mev){
  estimado<-100*(modelo_pm$VCV[,"animal"])/
    (modelo_pm$VCV[, focal_comp]+modelo_pm$VCV[,"animal"]+modelo_pm$VCV[,"units"])
  estimado2<-data.frame(posterior.mode (estimado))
  names (estimado2)<-"phylo heredability(h2)"
  rownames (estimado2)<-paste ("pm_H2", focal_comp, sep=".")
  return (estimado2)
}




