## ---------------------------
##
## Script name: data_processing_occurrence.r
##
## Purpose of script: This script analyses the relationship between the number of 
## interactions a species has and the number of sites in which is present.
##
## Author:Dr Nuria Galiana
##
## Date Created: 11-03-2022
##
## Email: galiana.nuria@gmail.com
##
##
## The raw data used within this script was provided by 
## Professor Tomas Roslin
## Swedish University of Agricultural Sciences
## Sweden
##
## Sections of the code below for reading raw data files were provided by Prof Roslin
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries

if(!require(magrittr)){ 
  install.packages("magrittr")
}

require(magrittr)
require(reshape2)
require(igraph)
require(bipartite)
require(AICcmodavg)
require(ggplot2)
require(broom)
require(data.table)

## The section below reads the original raw data
####### THE FOLLOWING CODE WAS PROVIDED BY TOMAS ROSLIN TO READ DATA FILES #######
unlink("./raw-data/csv", recursive = TRUE) 
unlink("./raw-data/rdata", recursive = TRUE)

# Change accordingly to your directory
setwd('./raw-data/')
source("../lib/format4R.r")
get_formatData("Salix_webs.csv")
setwd('..')

df_site <- readRDS("./raw-data/rdata/df_site.rds") 
str(df_site, strict.width="cut")

df_interact <- readRDS("./raw-data/rdata/df_interact.rds") 
df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR/df_interact$N_GALLS
str(df_interact, strict.width="cut")

site_interact <- merge(df_site, df_interact, by="REARING_NUMBER") 
head(site_interact)

############### SALGAL NETWORK (i.e. Salix-Galler interactions) ##########################

# Number of sites in which each species of the lower trophic level (i.e. Salix) occurs
output_occupancy_rsal <- NULL
for(i in unique(site_interact$RSAL)){
  occupation <- length(unique(site_interact[site_interact$RSAL==i,]$SITE))
  cur_out <- data.frame(trophic_level='RSAL', species=i, sites=occupation)
  if(is.null(output_occupancy_rsal)){
    output_occupancy_rsal <- cur_out
  }else{
    output_occupancy_rsal <- rbind(output_occupancy_rsal, cur_out)
  }
}

# Number of sites in which each species of the higher trophic level (i.e. Gallers) occurs
output_occupancy_rgal <- NULL
for(i in unique(site_interact$RGALLER)){
  occupation <- length(unique(site_interact[site_interact$RGALLER==i,]$SITE))
  cur_out <- data.frame(trophic_level='RGALLER', species=i, sites=occupation)
  if(is.null(output_occupancy_rgal)){
    output_occupancy_rgal <- cur_out
  }else{
    output_occupancy_rgal <- rbind(output_occupancy_rgal, cur_out)
  }
}

# Data frame with species occupancies of the salgal network
output_occupancy_salgal <- rbind(output_occupancy_rgal,output_occupancy_rsal)

# Here we create the metaweb for the salgal network 
mweb_salgal <- df_interact[,c("RSAL","RGALLER")] %>% unique 
igr_salgal <- data.frame(
  from = mweb_salgal$RSAL,
  to = mweb_salgal$RGAL
) %>% graph_from_data_frame(directed=TRUE)


# Calculating species diet breadth in the metaweb (i.e. number of species interactions) 
# for the lower trophic level (i.e. Salix)
output_diet_rsal <- NULL
for(d in unique(mweb_salgal$RSAL)){
  diet <- dim(mweb_salgal[mweb_salgal$RSAL==d,])[1]
  cur_out <- data.frame(trophic_level='RSAL', species=d, interactions=diet)
  if(is.null(output_diet_rsal)){
    output_diet_rsal <- cur_out
  }else{
    output_diet_rsal <- rbind(output_diet_rsal, cur_out)
  }
}

# Calculating species diet breadth in the metaweb (i.e. number of species interactions) 
# for the higher trophic level (i.e. Gallers)
output_diet_rgal_salgal <- NULL
for(d in unique(mweb_salgal$RGALLER)){
  diet <- dim(mweb_salgal[mweb_salgal$RGALLER==d,])[1]
  cur_out <- data.frame(trophic_level='RGALLER', species=d, interactions=diet)
  if(is.null(output_diet_rgal_salgal)){
    output_diet_rgal_salgal <- cur_out
  }else{
    output_diet_rgal_salgal <- rbind(output_diet_rgal_salgal, cur_out)
  }
}

# Data frame with species interactions information for the salgal network
output_diet_salgal <- rbind(output_diet_rsal,output_diet_rgal_salgal)

# Data frame with species interactions and occupancy information for the salgal network
occu_diet_salgal <- merge(output_diet_salgal, output_occupancy_salgal, by= c("species",'trophic_level')) 

############### GALPAR NETWORK (i.e. galler-parasitoid interactions) ##########################
# Number of sites in which each species of the galpar network occurs
output_occupancy_rpar <- NULL
for(i in unique(site_interact$RPAR)){
  occupation <- length(unique(site_interact[site_interact$RPAR==i,]$SITE))
  cur_out <- data.frame(trophic_level='RPAR', species=i, sites=occupation)
  if(is.null(output_occupancy_rpar)){
    output_occupancy_rpar <- cur_out
  }else{
    output_occupancy_rpar <- rbind(output_occupancy_rpar, cur_out)
  }
}

# Data frame with species occupancies of the galpar network
output_occupancy_galpar <- rbind(output_occupancy_rgal,output_occupancy_rpar)

# Here we create the metaweb for the network galpar
id <- df_interact$RPAR!="none"
mweb_galpar <- df_interact[id,c("RPAR","RGALLER","PAR_RATE")] %>% unique 
igr_salpar <- data.frame(
  from = mweb_galpar$RGAL,
  to = mweb_galpar$RPAR
) %>% graph_from_data_frame(directed=TRUE)

# Calculating species diet breadth in the metaweb (i.e. number of species interactions) 
# for the lower trophic level (i.e. Gallers)
output_diet_rgal_galpar <- NULL
for(d in unique(mweb_galpar$RGALLER)){
  diet <- length(unique(mweb_galpar[mweb_galpar$RGALLER==d,]$RPAR))
  mean_par_rate <-  mean(mweb_galpar[mweb_galpar$RGALLER==d,]$PAR_RATE)
  sd_par_rate <-  sd(mweb_galpar[mweb_galpar$RGALLER==d,]$PAR_RATE)
  cur_out <- data.frame(trophic_level='RGALLER', species=d, interactions=diet,mean_parasitism.rate=mean_par_rate, sd_parasitism.rate=sd_par_rate)
  if(is.null(output_diet_rgal_galpar)){
    output_diet_rgal_galpar <- cur_out
  }else{
    output_diet_rgal_galpar <- rbind(output_diet_rgal_galpar, cur_out)
  }
}

# Calculating species diet breadth in the metaweb (i.e. number of species interactions) 
# for the higher trophic level (i.e. Parasitoids)
output_diet_rpar <- NULL
for(d in unique(mweb_galpar$RPAR)){
  diet <- length(unique(mweb_galpar[mweb_galpar$RPAR==d,]$RGALLER))
  mean_par_rate <-  mean(mweb_galpar[mweb_galpar$RPAR==d,]$PAR_RATE)
  sd_par_rate <-  sd(mweb_galpar[mweb_galpar$RPAR==d,]$PAR_RATE)
  cur_out <- data.frame(trophic_level='RPAR', species=d, interactions=diet, mean_parasitism.rate=mean_par_rate, sd_parasitism.rate=sd_par_rate)
  if(is.null(output_diet_rpar)){
    output_diet_rpar <- cur_out
  }else{
    output_diet_rpar <- rbind(output_diet_rpar, cur_out)
  }
}

# Data frame with species interactions information for the galpar network
output_diet_galpar <- rbind(output_diet_rpar,output_diet_rgal_galpar)

# Data frame with species interactions and occupancy information for the galpar network
occu_diet_galpar <- merge(output_diet_galpar, output_occupancy_galpar, by= c("species",'trophic_level')) 

# Calculating the total number of interactions per species

number_interactions_pred <- unique(occu_diet_galpar[occu_diet_galpar$trophic_level=='RPAR',c(1,3)])

number_interactions_prey <- unique(occu_diet_galpar[occu_diet_galpar$trophic_level=='RGALLER',c(1,3)])

number_interactions <- rbind(number_interactions_pred, number_interactions_prey)

# Calculating species occupancy

number_sites_pred <- unique(occu_diet_galpar[occu_diet_galpar$trophic_level=='RPAR',c(1,6)])

number_sites_prey <- unique(occu_diet_galpar[occu_diet_galpar$trophic_level=='RGALLER',c(1,6)])


# We Check the match between parasitoid occupancy and galler occupancy. 

prey_occupancy <- NULL

for(sp in unique(site_interact$RPAR)){
  all_prey <- unique(site_interact[site_interact$RPAR==sp,]$RGALLER)
  sites_pred <- unique(site_interact[site_interact$RPAR==sp,]$SITE)
  for(prey in unique(all_prey)){
    sites_prey <- unique(site_interact[site_interact$RGALLER==prey,]$SITE)
    followed <- length(sites_pred[which(sites_pred %in% sites_prey)])
    sites_prey_alone <- length(setdiff(sites_prey,sites_pred))
    cur_out <- data.frame(species=sp, sp_gall=prey, sites_unfollowed=sites_prey_alone,sites_followed=followed,proportion_followed=followed/length(sites_prey),sites_par=length(sites_pred),sites_gall=length(sites_prey),ratio_occu=length(sites_pred)/length(sites_prey))
    if(is.null(prey_occupancy)){
      prey_occupancy <- cur_out
    }else{
      prey_occupancy <- rbind(prey_occupancy, cur_out)
    }
  }
}

# We merge the information of galler occupancies with parasioid occupancies
occu_prey <- merge(occu_diet_galpar,prey_occupancy, by='species')

########################################################################
# Computing the netwwork of co-occurrences for predators
output_cooc_per_site_pred <- NULL
output_cooccurrence_pred <- NULL

for(sp in unique(site_interact$RPAR)){
  sp_sites <- unique(site_interact[site_interact$RPAR==sp,]$SITE)
  for(sit in unique(sp_sites)){
    local_com <- site_interact[site_interact$SITE==sit,]
    cur_potential_prey <- unique(local_com$RGALLER)
    cur_out <- data.frame(site=sit, species=sp, potential_prey=cur_potential_prey)
    if(is.null(output_cooc_per_site_pred)){
      output_cooc_per_site_pred <- cur_out
    }else{
      output_cooc_per_site_pred <- rbind(output_cooc_per_site_pred, cur_out)
    }
  }
  total_potential <- length(unique(output_cooc_per_site_pred[output_cooc_per_site_pred$species==sp,]$potential_prey))
  potential_prey<- unique(output_cooc_per_site_pred[output_cooc_per_site_pred$species==sp,]$potential_prey)
  
  cur_out_fin <- data.frame(species=sp, total_pot_prey=total_potential, coocurring_species=potential_prey)
  if(is.null(output_cooccurrence_pred)){
    output_cooccurrence_pred <- cur_out_fin
  }else{
    output_cooccurrence_pred <- rbind(output_cooccurrence_pred, cur_out_fin)
  }
}

# Computing the netwwork of co-occurrences for predators
output_cooc_per_site_prey <- NULL
output_cooccurrence_prey <- NULL

for(sp in unique(site_interact$RGALLER)){
  sp_sites <- unique(site_interact[site_interact$RGALLER==sp,]$SITE)
  for(sit in unique(sp_sites)){
    local_com <- site_interact[site_interact$SITE==sit,]
    cur_potential_pred <- unique(local_com$RPAR)
    cur_out <- data.frame(site=sit, species=sp, potential_pred=cur_potential_pred)
    if(is.null(output_cooc_per_site_prey)){
      output_cooc_per_site_prey <- cur_out
    }else{
      output_cooc_per_site_prey <- rbind(output_cooc_per_site_prey, cur_out)
    }
  }
  total_potential <- length(unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred))
  potential_pred<- unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred)
  
  cur_out_fin <- data.frame(species=sp, total_pot_pred=total_potential, coocurring_species=potential_pred)
  if(is.null(output_cooccurrence_prey)){
    output_cooccurrence_prey <- cur_out_fin
  }else{
    output_cooccurrence_prey <- rbind(output_cooccurrence_prey, cur_out_fin)
  }
}


# We calculate the proportion of sites in which two species co-occur
proportion_coocurrence <- NULL

for(sp in unique(output_cooc_per_site_pred$species)){
  total_n_sites <- length(unique(output_cooc_per_site_pred$site))
  cur_sp <- output_cooc_per_site_pred[output_cooc_per_site_pred$species==sp,]
  for(pot_prey in unique(cur_sp$potential_prey)){
    n_cooc <- dim(cur_sp[cur_sp$potential_prey==pot_prey,])[1]
    prop_cooc <- n_cooc#/total_n_sites
    
    cur_out <- data.frame(species=sp, potential_prey=pot_prey, prop.cooc=prop_cooc)
    
    if(is.null(proportion_coocurrence)){
      proportion_coocurrence <- cur_out
    }else{
      proportion_coocurrence <- rbind(proportion_coocurrence, cur_out)
    }
  }
}


# N_alpha: sum of all consumers co-occurrences frequencies

proportion_coocurrence <- proportion_coocurrence %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(N_alpha = sum(prop.cooc),
                relative_N_alpha = sum(prop.cooc)/length(potential_prey),
                number_potential_interactions = length(unique(potential_prey)))

# N_i: sum of all resources co-occurrences frequencies

proportion_coocurrence <- proportion_coocurrence %>%
  dplyr:: group_by(potential_prey) %>%
  dplyr::mutate(N_i = sum(prop.cooc),
                relative_N_i = sum(prop.cooc)/length(species),
                number_potential_interactions_resources = length(unique(species)))


interaction_proportion <- merge(number_interactions, proportion_coocurrence, by="species") 

# We built a model to infer interactions based on frequency of co-occurrence

p <- 0.121

expected_int_model <- NULL

for(pred in unique(interaction_proportion$species)){
  cur_pred <- interaction_proportion[interaction_proportion$species==pred,]
  N_alpha <- unique(cur_pred$N_alpha)
  for(prey in unique(cur_pred$potential_prey)){
    cur_pred_prey <- cur_pred[cur_pred$potential_prey==prey,]
    n_i_alpha <- cur_pred_prey$prop.cooc
    prob_interaction_expected <- (1-(1-p)^n_i_alpha)/(1-(1-p)^N_alpha)
    prob_interaction <- (1-(1-p)^n_i_alpha)
    
    cur_out <- data.frame(species=pred, potential_prey=prey, prob_expected= prob_interaction_expected, prob_int = prob_interaction)
    
    if(is.null(expected_int_model)){
      expected_int_model <- cur_out
    }else{
      expected_int_model <- rbind(expected_int_model, cur_out)
    }
  }
}

sum(expected_int_model$prob_expected)
sum(number_interactions_pred$interactions)


## We now generate 100 realisations using the probability of interaction 

replicates <- 100
realisations_indegree <- NULL
realisations_outdegree <- NULL

for(r in 1:replicates){
  int_mat <- xtabs(prob_int~species+potential_prey, data=expected_int_model)
  random_mat <- matrix(runif(dim(int_mat)[1]*dim(int_mat)[2]), ncol=dim(int_mat)[2])
  random_mat[random_mat>int_mat]<-0 #we define the probability of having links in the matrix by comparing the numbers given by the uniform distribution with the probability of interaction
  random_mat[random_mat!=0]<-1
  colnames(random_mat) <- colnames(int_mat)
  rownames(random_mat) <- rownames(int_mat)
  
  if(any(rowSums(random_mat)==0)){
    for(pred in rownames(random_mat[rowSums(random_mat)==0,])){
      new_interact <- expected_int_model[order(expected_int_model[expected_int_model$species==pred,]$prob_int,decreasing=TRUE),][1,]$potential_prey
      random_mat[pred , new_interact] <- 1
    }
  }
  
  
  cur_out_indegree <- data.frame(species=rownames(random_mat),indegree=rowSums(random_mat), rep=r)
  cur_out_indegree <- merge(number_interactions_pred, cur_out_indegree, by='species')
  
  
  cur_out_outdegree <- data.frame(species=colnames(random_mat),outdegree=colSums(random_mat),rep=r)
  cur_out_outdegree <- merge(number_interactions_prey, cur_out_outdegree, by='species')
  
  if(is.null(realisations_outdegree)){
    realisations_outdegree <- cur_out_outdegree
  }else{
    realisations_outdegree <- rbind(realisations_outdegree, cur_out_outdegree)
  }
  
  if(is.null(realisations_indegree)){
    realisations_indegree <- cur_out_indegree
  }else{
    realisations_indegree <- rbind(realisations_indegree, cur_out_indegree)
  }
}

write.csv(realisations_indegree, file = 'random_realisations_expected_indegree_galpar.csv')
write.csv(realisations_outdegree, file = 'random_realisations_expected_outdegree_galpar.csv')

# Sum all the probabilities of a given predator to obtain the expected number of links

expected_int_model_indegree <- expected_int_model %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(expected_indegree = sum(prob_int))

expected_realised_indegree <- merge(number_interactions_pred, expected_int_model_indegree, by='species')

expected_realised_indegree_freq <- merge(expected_realised_indegree, interaction_proportion[,c('species','relative_N_alpha')], by='species')

colnames(expected_realised_indegree_freq)[7] <- 'mean_cooccurrence'

write.csv(expected_realised_indegree_freq, file = 'expected_realised_indegree_galpar.csv')

# Sum all the probabilities of a given predator to obtain the expected number of links

expected_int_model_outdegree <- expected_int_model[,c('potential_prey','prob_expected')] %>%
  dplyr::group_by(potential_prey) %>%
  dplyr::mutate(expected_outdegree = sum(prob_expected))

colnames(expected_int_model_outdegree)[1] <- 'species'

expected_realised_outdegree <- merge(number_interactions_prey, expected_int_model_outdegree, by='species')

interaction_proportion_outdegree <- interaction_proportion[,c('potential_prey','relative_N_i')]
colnames(interaction_proportion_outdegree)[1] <- 'species'

expected_realised_outdegree_freq <- merge(expected_realised_outdegree, interaction_proportion_outdegree, by='species')

colnames(expected_realised_outdegree_freq)[5] <- 'mean_cooccurrence'

write.csv(expected_realised_outdegree_freq, file = 'expected_realised_outdegree_galpar.csv')

# Here we have the occupancy distributions for both trophic levels.
occupancy_pred <- unique(number_sites_pred)

occupancy_prey <- unique(number_sites_prey)

write.csv(occupancy_pred, file = 'occupancy_pred_galpar.csv')
write.csv(occupancy_prey, file = 'occupancy_prey_galpar.csv')

# Here we have the outputs from both networks the realised and the co-occurrence network.

network_real_pred <- unique(number_interactions_pred)
network_real_prey <- unique(number_interactions_prey)

occupancy_pred <- unique(number_sites_pred)
occupancy_prey <- unique(number_sites_prey)

network_coocurrence_pred <- unique(output_cooccurrence_pred[,1:2])
network_coocurrence_prey <- unique(output_cooccurrence_prey[,1:2])

write.csv(network_real_pred, file = 'network_real_pred_galpar.csv')
write.csv(network_real_prey, file = 'network_real_prey_galpar.csv')

write.csv(occupancy_pred, file = 'occupancy_pred_galpar.csv')
write.csv(occupancy_prey, file = 'occupancy_prey_galpar.csv')

write.csv(network_coocurrence_pred, file = 'network_coocurrence_pred_galpar.csv')
write.csv(network_coocurrence_prey, file = 'network_coocurrence_prey_galpar.csv')


# We now built a model to infer interactions without using frequency of co-occurrence to predict interactions.
# This corresponds to a random pruning of the network of co-occurrences.

p <- 0.29 # p here is the percentage of interactions that where realised from the co-occurrence network to the network of biotic interactions

expected_int_null_model <- NULL

for(pred in unique(interaction_proportion$species)){
  cur_pred <- interaction_proportion[interaction_proportion$species==pred,]
  N_alpha <- unique(cur_pred$N_alpha)
  for(prey in unique(cur_pred$potential_prey)){
    cur_pred_prey <- cur_pred[cur_pred$potential_prey==prey,]
    n_i_alpha <- 1
    prob_interaction_expected <- (1-(1-p)^n_i_alpha)/(1-(1-p)^N_alpha)
    prob_interaction <- (1-(1-p)^n_i_alpha)
    
    cur_out <- data.frame(species=pred, potential_prey=prey, prob_expected= prob_interaction_expected, prob_int = prob_interaction)
    
    if(is.null(expected_int_null_model)){
      expected_int_null_model <- cur_out
    }else{
      expected_int_null_model <- rbind(expected_int_null_model, cur_out)
    }
  }
}

sum(expected_int_null_model$prob_expected)
sum(number_interactions_pred$interactions)

# We now generate 100 realisations using the probability of interaction

replicates <- 100
realisations_indegree_null_model <- NULL
realisations_outdegree_null_model <- NULL

for(r in 1:replicates){
  int_mat <- xtabs(prob_int~species+potential_prey, data=expected_int_null_model)
  random_mat <- matrix(runif(dim(int_mat)[1]*dim(int_mat)[2]), ncol=dim(int_mat)[2])
  random_mat[random_mat>int_mat]<-0 #we define the probability of having links in the matrix by comparing the numbers given by the uniform distribution with the probability of interaction
  random_mat[random_mat!=0]<-1
  colnames(random_mat) <- colnames(int_mat)
  rownames(random_mat) <- rownames(int_mat)
  
  if(any(rowSums(random_mat)==0)){
    for(pred in rownames(random_mat[rowSums(random_mat)==0,])){
      new_interact <- expected_int_null_model[order(expected_int_null_model[expected_int_null_model$species==pred,]$prob_int,decreasing=TRUE),][1,]$potential_prey
      random_mat[pred , new_interact] <- 1
    }
  }
  
  
  cur_out_indegree <- data.frame(species=rownames(random_mat),indegree=rowSums(random_mat), rep=r)
  cur_out_indegree <- merge(number_interactions_pred, cur_out_indegree, by='species')
  
  
  cur_out_outdegree <- data.frame(species=colnames(random_mat),outdegree=colSums(random_mat),rep=r)
  cur_out_outdegree <- merge(number_interactions_prey, cur_out_outdegree, by='species')
  
  if(is.null(realisations_outdegree_null_model)){
    realisations_outdegree_null_model <- cur_out_outdegree
  }else{
    realisations_outdegree_null_model <- rbind(realisations_outdegree_null_model, cur_out_outdegree)
  }
  
  if(is.null(realisations_indegree_null_model)){
    realisations_indegree_null_model <- cur_out_indegree
  }else{
    realisations_indegree_null_model <- rbind(realisations_indegree_null_model, cur_out_indegree)
  }
}


write.csv(realisations_indegree_null_model, file = 'null_random_realisations_expected_indegree_galpar.csv')
write.csv(realisations_outdegree_null_model, file = 'null_random_realisations_expected_outdegree_galpar.csv')

