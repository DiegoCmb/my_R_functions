library (MCMCglmm)
library(phangorn)
library (dplyr)

weighting_and_subsetting<-function (database, weight_approach='inverse_sample_size'){
  # to estimate the weights.
  # sample_sizes below 3 are removed since the estimation of
  # mev substract 3 from sample_sizes generating Inf outcomes
  database<-subset(database, database$sample_size != 'NA') # because we want to weight by sample sizes
  # we will have to remove obs with out this info
  if (weight_approach == "inverse_sample_size"){
    database$mev<-1/(database$sample_size)
    database$weight<-1/database$mev
    database$mesd<-sqrt(database$mev)
    return(database)
  }else {
    return(database)
  }
}

force.ultrametric<-function(tree, method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

mcmcglmm_dignostics<-function (modelo_fix, modelo_random, database, mev_col = NULL, gelman_diag="True", nitt_x= 13000, thin_x= 10 , 
                               burnin_x= 3000, prior_x, save_output = "False", saving_path, file_name ){
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
  if (gelman_diag == "False"){
    m0.1<-MCMCglmm (as.formula(modelo_fix), random= as.formula (modelo_random), family = "gaussian", 
                    data = database, mev= mev_col , verbose= T,
                    nitt=nitt_x, thin=thin_x, burnin= burnin_x,pr=T,prior=prior_x)
    # to save data in .RData
    if (save_output == "True"){
      path_file_to_save<-file.path(saving_path, file_name)
      save(m0.1, file= path_file_to_save)
      return (m0.1)
    }else {
      return (m0.1)
    }
    
    # Running multiple models for gelman diagnostic
  } else {
    m0.1<-MCMCglmm (as.formula(modelo_fix), random= as.formula (modelo_random), family = "gaussian", 
                    data = database, mev= mev_col, verbose= T,
                    nitt=nitt_x, thin=thin_x, burnin= burnin_x,pr=T,prior=prior_x)
    
    m0.2<-MCMCglmm (as.formula(modelo_fix), random= as.formula (modelo_random), family = "gaussian", 
                    data = database, mev= mev_col , verbose= T,
                    nitt=nitt_x, thin=thin_x, burnin= burnin_x,pr=T,prior=prior_x)
    
    m0.3<-MCMCglmm (as.formula(modelo_fix), random= as.formula (modelo_random), family = "gaussian", 
                    data = database, mev= mev_col, verbose= T,
                    nitt=nitt_x, thin=thin_x, burnin= burnin_x,pr=T,prior=prior_x)
    
    gel.diag_SOL<-gelman.diag(list(m0.1$Sol[,1],m0.2$Sol[,1],m0.3$Sol[,1]))
    gel.diag_VCV<-gelman.diag(list(m0.1$VCV[,c(1,3)],m0.2$VCV[,c(1,3)],m0.3$VCV[,c(1,3)]))
    gel.diag_DEV<-gelman.diag(list(m0.1$Deviance,m0.2$Deviance,m0.3$Deviance))
    
    m0.1_DIC<-m0.1$DIC
    m0.2_DIC<-m0.2$DIC
    m0.3_DIC<-m0.3$DIC
    
    list_of_DIC<-list (m0.1_DIC, m0.2_DIC, m0.3_DIC)
    names(list_of_DIC)<-c("DIC_model_1", "DIC_model_2","DIC_model_3")
    
    model_list_output<-list ("/////////////1) First_model  /////////////", m0.1, # the number within parethesis indicates the parsing indexation in the list
                             "/////////////3) Second_model /////////////", m0.2, 
                             "/////////////4) Third_model  /////////////", m0.3,
                             "/////////////6) Gelman convergence diagnostic SOL /////////////", gel.diag_SOL,
                             "/////////////8) Gelman convergence diagnostic SOL /////////////", gel.diag_VCV,
                             "/////////////10) Gelman convergence diagnostic DEV /////////////", gel.diag_DEV,
                             "/////////////12) DIC for each model /////////////", list_of_DIC)
    
    if (save_output == "True"){
      path_file_to_save<-file.path(saving_path, file_name)
      save(model_list_output, file= path_file_to_save)
      return (model_list_output)
    }else{
      return (model_list_output)
    }
  }
}
