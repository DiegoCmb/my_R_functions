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



########## Heterogeneity and phylogenetic heritability #####
# based on Scripts from Nakagawa and 
# https://github.com/daniel1noble/metaAidR/blob/master/R/I2.R
# estimating 
# mev<-mev.all
# W<-1/mev  # is this correct?????? this is the formula when is based on 
# # measurment of variance, but I dont know if with sample sizes
# # should follow the same logic. 
# s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()

# heterogeneity for focal factor without considering phylogeny
mm_I2.factor<-function (mcmc_modelo_output, random_focal_component, mev){
  ##################### About heterogenity #####################################
  #
  # model_mm is an MCMCglmm object
  # random_focal_component: is the random component that is of interest
  # mev: the component of within variation and it is associated to the database
  # used to generate the MCMCglmm object (model)
  # The reliabilestimates ity of a general trend depnds on the degree of 
  # consistency among studies (heterogeneity). 
  # Two parameters Q and I2. Santos and Nakagawa 2012
  # indicates that neither of both can be use in a multilevel meta-analysi
  # but suggest a new I2 implementation that can be implemented in such models
  # I2 = 25, 50, 75% are low, moderate and high heterogeneity respectively
  # estimation I2specie level = variance.sp/ variance.total 
  # or I2study = variance.study/variance.total
  # variance_tot= var.phylo+var.study+var.sp+var.withinstudy+var.error
  
  model_name<- deparse(substitute (mcmc_modelo_output))
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(mcmc_modelo_output$VCV[,random_focal_component])/
    (mcmc_modelo_output$VCV[,random_focal_component]+mcmc_modelo_output$VCV[,"units"]+s2)
  estimado2<-data.frame("model_name" = model_name, 'heterogenity_I2_random_comp' = posterior.mode(estimado))
  rownames(estimado2)<-paste('mm_I2', random_focal_component, sep="_" )
  return (estimado2)
}

# heterogeneity for the same focal factor but considering phylogeny
pm_I2.focal.phylo<-function (mcmc_modelo_output, random_focal_component, mev){ 
  ##################### About heterogenity #####################################
  #
  # model_mm is an MCMCglmm object
  # random_focal_component: is the random component that is of interest
  # mev: the component of within variation and it is associated to the database
  # used to generate the MCMCglmm object (model)
  # The reliabilestimates ity of a general trend depnds on the degree of 
  # consistency among studies (heterogeneity). 
  # Two parameters Q and I2. Santos and Nakagawa 2012
  # indicates that neither of both can be use in a multilevel meta-analysi
  # but suggest a new I2 implementation that can be implemented in such models
  # I2 = 25, 50, 75% are low, moderate and high heterogeneity respectively
  # estimation I2specie level = variance.sp/ variance.total 
  # or I2study = variance.study/variance.total
  # variance_tot= var.phylo+var.study+var.sp+var.withinstudy+var.error
  model_name<- deparse(substitute (mcmc_modelo_output))
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(mcmc_modelo_output$VCV[,random_focal_component]+mcmc_modelo_output$VCV[,"animal"])/
    (mcmc_modelo_output$VCV[, random_focal_component]+mcmc_modelo_output$VCV[,"animal"]+
       mcmc_modelo_output$VCV[,"units"]+s2)
  estimado2<-data.frame("model_name" = model_name, 'heterogenity_I2_random_comp' = posterior.mode(estimado))
  rownames(estimado2)<-paste('pm_I2', random_focal_component, sep="_" )
  return (estimado2)
}

# heterogeneity due to phylogeny considering the factor
pm_I2.phylo<-function (mcmc_modelo_output, random_focal_component, mev){
  ##################### About heterogenity #####################################
  #
  # The reliabilestimates ity of a general trend depnds on the degree of 
  # consistency among studies (heterogeneity). 
  # Two parameters Q and I2. Santos and Nakagawa 2012
  # indicates that neither of both can be use in a multilevel meta-analysi
  # but suggest a new I2 implementation that can be implemented in such models
  # I2 = 25, 50, 75% are low, moderate and high heterogeneity respectively
  # estimation I2specie level = variance.sp/ variance.total 
  # or I2study = variance.study/variance.total
  # variance_tot= var.phylo+var.study+var.sp+var.withinstudy+var.error
  #####
  model_name<- deparse(substitute (mcmc_modelo_output))
  W<- 1/mev
  s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) # measurment of variance ()
  estimado<-100*(mcmc_modelo_output$VCV[,"animal"])/
    (mcmc_modelo_output$VCV[, random_focal_component]+mcmc_modelo_output$VCV[,"animal"]+
       mcmc_modelo_output$VCV[,"units"]+s2)
  estimado2<-data.frame("model_name" = model_name, 'heterogenity_I2_phylo' = posterior.mode(estimado))
  rownames(estimado2)<-paste('pm_I2', random_focal_component, sep="_" )
  return (estimado2)
}


# Phylogenetic heredability, similar to lambda pagel
pm_H2.phylo<-function (mcmc_modelo_output, random_focal_component, mev){
  
  ##################### About phylogenetic heritability ######################
  # From estimating variance_tot it can be get the phylogenetic heritability
  # that is another measurment indicating phylogenetic relatedness,
  
  # H2= var.phylo/variance_total
  
  # H2 = 0 no phylogenetic relatedness among effect sizes or traits, H2 = 1
  # indicates that effect sizes or traits values among species are exactly 
  # proportional to theri phylogenetic relatedness (Lynch 1991).
  model_name<- deparse(substitute (mcmc_modelo_output))
  estimado<-100*(mcmc_modelo_output$VCV[,"animal"])/
    (mcmc_modelo_output$VCV[, random_focal_component]+mcmc_modelo_output$VCV[,"animal"]+mcmc_modelo_output$VCV[,"units"])
  estimado2<-data.frame("model_name" = model_name, 'heredability_H2_phylo' = posterior.mode(estimado))
  rownames(estimado2)<-paste('pm_I2', random_focal_component, sep="_" )
  return (estimado2)
}

#### Final function that handle the estimation of I2 and H2 using the previous 
# functions
mcmc_I2_H2 <-
  function(model_pm, moderator, mev, analysis = "traditional") {
    # Administrate fungtions to get I2 and H2 from MCMCglmm objects
    # with 1 random factor and the phylogenetic component of variation (animal)  
    # if needed. Returns the output as row in a dataframe object
    # it uses several functions  mm_I2_focal, pm_I2.focal.phylo, pm_I2.phylo
    # pm_H2.phylo
    
    model_pm2 <- deparse(substitute(model_pm))
    
    if ("traditional" == analysis) {
      mm_I2_focal <- mm_I2.factor(model_pm, moderator, mev)
      dataframe1 <- data.frame("models_names" = model_pm2,
                               "I2_moderator" = mm_I2_focal[2],
                               "heterogenity_I2_phylo" = NA,
                               "heredability_H2_phylo"= NA)
      rownames(dataframe1) <- c()
      return (dataframe1)
    } else if ("phylogenetic" == analysis) {
      pm_I2_focal_phylo <- pm_I2.focal.phylo (model_pm, moderator, mev)
      pm_I2_phylo <- pm_I2.phylo(pm_2_he_study, "animal", mev)
      pm_H2_phylo <- pm_H2.phylo(pm_2_he_study, moderator, mev)
      
      dataframe1 <- data.frame(
        "models_names" = model_pm2,
        "I2_moderator" = pm_I2_focal_phylo[2],
        "I2_phylo" = pm_I2_phylo[2],
        "H2_phylo" = pm_H2_phylo[2]
      )
      rownames(dataframe1) <- c()
      return (dataframe1)
    } else{
      print ('Error: wrong name of analysis. Choose between traditional or phylogentic')
    }
  }

##################### Parsing functions
library(devtools) # install.packages("devtools")
library (broom.mixed) #install_github("bbolker/broom.mixed")

row_summary_mcmcglmm<-function(mcmc_modelo_output, moderator, extra_info){
  # allow to generate a summary table from MCMCglmm by using broom.mixed
  # needs to be feeded with the MCMCglmm model. 
  model_name <-deparse(substitute(mcmc_modelo_output))
  
  database1<-tidy(mcmc_modelo_output,  effects = "fixed", conf.int = TRUE,
                  conf.level = 0.95, conf.method = "HPDinterval",  ess = TRUE )
  lista_DIC<-list("DIC" = rep (mcmc_modelo_output$DIC, dim (database1)[1]))
  lista_pMCMC<-list("pMCMC"= summary(mcmc_modelo_output)$solutions[,5])
  models_names<-rep(model_name, dim(database1)[1])
  database2<-data.frame("moderator" = moderator, "model_type"= extra_info, 
                        models_names, database1, lista_pMCMC, lista_DIC)
  rownames(database2) <- c()
  return (database2)
}

#### Putting all the parsed information in one output

complete_summary_mcmcglmm2 <-
  function (mcmc_modelo_output,
            moderator,
            mev,
            analysis) {
    model_name <- deparse(substitute(mcmc_modelo_output))
    print (model_name)
    if (analysis == "traditional") {
      trad <-
        mcmc_I2_H2(mcmc_modelo_output, moderator, mev, analysis = "traditional")
      database1 <-
        row_summary_mcmcglmm(mcmc_modelo_output, moderator, analysis)
      database1 <- database1[, -1]
      database2 <- data.frame(database1, trad[-1])
      
      return (database2)
    } else if (analysis == "phylogenetic") {
      phylo <-
        mcmc_I2_H2(mcmc_modelo_output, moderator, mev, analysis = "phylogenetic")
      database1 <-
        row_summary_mcmcglmm(mcmc_modelo_output, moderator, analysis)
      database2 <- data.frame(database1, phylo[-1])
      return (database2)
    } else if (analysis == "null") {
      database1 <-
        row_summary_mcmcglmm(mcmc_modelo_output, moderator, analysis)
      database2 <-
        data.frame(
          database1,
          "heterogenity_I2_random_comp" = NA,
          "heterogenity_I2_phylo" = NA,
          "heredability_H2_phylo" = NA
        )
      return (database2)
    } else{
      print ("Error. choose between traditional or phylogenetic")
    }
  }

# examples... 
a<-complete_summary_mcmcglmm2(mm_3_he_study_marker, "citation", dato.all$database$mev, analysis = "traditional")
a
# b<-comlete_summary_mcmcglmm(pm_4_he_study_marker, "citation", dato.all$database$mev, analysis = "phylogenetic")
# b


##################### About repeatitivity #####################################

# NEED TO WORK ON THIS 

#http://www.wildanimalmodels.org/tiki-download_wiki_attachment.php?attId=4




