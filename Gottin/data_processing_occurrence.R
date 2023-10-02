## ---------------------------
##
## Script name: data_processing_newproject.r
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
## Notes:
##
## The raw data used within this script was provided by 
## Dr Ingo Grass
## Georg-August-Universitat
## Gottingen, Germany
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(bipartite)
source("../utils.r")

########################################################################################
## The section below reads the original raw data

## network data
ant_raw <- read.csv("./raw-data/host_para_all_interactions.csv", head=T, sep=";")
mut_raw <- read.csv("./raw-data/plant_poll_all_interactions.csv", head=T, sep=";")

## metadata environmental variables (e.g., size of grassland fragments)
metadata <- read.csv("metadata.csv", head=T, sep=",")

##########
## create individual webs for each study site
ant_nets = frame2webs(ant_raw, varnames=c("Genus.Species","P1.Species.Genus","Site","P1.cells"))
mut_nets = frame2webs(mut_raw, varnames=c("P.Genus.Species","Genus.Species","Site"))


################################## ANTAGONISTIC NETWORK #############################

## Calculating species occupancy ANT_NETS

occupancy <- NULL

for(p in unique(metadata$site)){
  local_com <- as.data.frame(ant_nets[[p]])
  
  species_present <- as.data.frame(names(local_com))
  species_present <- setNames(species_present, c("species"))
  species_present$trophic_level <- 'consumer'
  
  resources_present <- as.data.frame(rownames(local_com))
  resources_present <- setNames(resources_present, c("species"))
  resources_present$trophic_level <- 'resource'
  
  
  cur_data <- data.frame(site=p,rbind(species_present,resources_present))
  
  if(is.null(occupancy)){
    occupancy <- cur_data
  }else{
    occupancy <- rbind(occupancy, cur_data)
  }
}

number_sites <- NULL
for(sp in unique(occupancy$species)){
  cur_sp <- occupancy[occupancy$species==sp,]
  num_sites <- length(cur_sp$site)
  cur_data <- data.frame(species=sp,trophic_level=cur_sp$trophic_level,sites=num_sites)
  
  if(is.null(number_sites)){
    number_sites <- cur_data
  }else{
    number_sites <- rbind(number_sites, cur_data)
  }
}

number_sites <- unique(number_sites[number_sites$sites!=0,])

number_sites_pred <- number_sites[number_sites$trophic_level=='consumer',]
number_sites_prey <- number_sites[number_sites$trophic_level=='resource',]

####### ATTENTION #######
##### it was detected that networks number 7, 18, 20, 29 of ant_nets had a column without a name. 
##### This column was removed from the analyses in each network respectively 
##### because we couldn't know which species that column was representing
ant_nets[[7]] <- ant_nets[[7]][,-1]
ant_nets[[18]] <- ant_nets[[18]][,-1]
ant_nets[[20]] <- ant_nets[[20]][,-1]
ant_nets[[29]] <- ant_nets[[29]][,-1]

## create metaweb for each interaction type
ant_raw$dummy = 1

ant_meta = frame2webs(ant_raw, varnames=c("Genus.Species","P1.Species.Genus","dummy","P1.cells"))[[1]]

sum(ant_meta)

## Calculating total number of interactions per species (of the higher trophic level)
metaweb <- graph_from_incidence_matrix(ant_meta, directed=T, mode='out')

number_interactions_pred <- as.data.frame(igraph::degree(metaweb, V(metaweb)[which(V(metaweb)$type == TRUE)]$name, mode='in'))

number_interactions_pred <- tibble::rownames_to_column(number_interactions_pred, "species")

number_interactions_pred <- setNames(number_interactions_pred, c("species","interactions"))

number_interactions_prey <- as.data.frame(igraph::degree(metaweb, V(metaweb)[which(V(metaweb)$type == FALSE)]$name, mode='out'))

number_interactions_prey <- tibble::rownames_to_column(number_interactions_prey, "species")

number_interactions_prey <- setNames(number_interactions_prey, c("species","interactions"))

number_interactions <- rbind(number_interactions_pred, number_interactions_prey)

real_interactions <- as.data.frame(get.edgelist(metaweb))
real_interactions <- setNames(real_interactions, c("resource","consumer"))

site_interaction <- merge(number_interactions, number_sites, by="species") 


# We calculate the number of times two species interact relative to the number of times they co-occur
output_prop.inter <- NULL
output_site.inter <- NULL

for(con in unique(real_interactions$consumer)){
  cur_sp <- real_interactions[real_interactions$consumer==con,]
  for(res in unique(cur_sp$resource)){
    consumer_sites <- unique(occupancy[occupancy$species==con,]$site)
    resource_sites <- unique(occupancy[occupancy$species==res,]$site)
    number_cooc <- sum(consumer_sites %in% resource_sites, na.rm=TRUE)
    for(p in unique(intersect(consumer_sites,resource_sites))){ 
      n <- as.data.frame(ant_nets[[p]])
      
      local_com <- igraph::graph_from_incidence_matrix(n, directed = T, mode='out')
      local_com <- igraph::delete.vertices(local_com, which(igraph::degree(local_com,  igraph::V(local_com)) == 0))
      
      cur_interactions <- as.data.frame(igraph::get.edgelist(local_com))
      cur_interactions <- setNames(cur_interactions, c("resource","consumer"))
      
      interactions <- dim(cur_interactions[cur_interactions$consumer==con & cur_interactions$resource==res,])[1]
      if(interactions>1){
        interactions <- 1
      } 
      
      cur_out <- data.frame(consumer=con, resource=res, site=p, interact=interactions)
      if(is.null(output_site.inter)){
        output_site.inter <- cur_out
      }else{
        output_site.inter <- rbind(output_site.inter, cur_out)
      }
    }
    
    sites_interact <- dim(output_site.inter[output_site.inter$consumer==con & output_site.inter$resource==res & output_site.inter$interact==1,])[1]
    
    cooccurring_sites <- dim(output_site.inter[output_site.inter$consumer==con & output_site.inter$resource==res,])[1]
    
    cur_proportion <- data.frame(consumer=con, resource=res, coocurring_sites=cooccurring_sites, interacting_sites=sites_interact, proportion_interact=sites_interact/cooccurring_sites)
    
    if(is.null(output_prop.inter)){
      output_prop.inter <- cur_proportion
    }else{
      output_prop.inter <- rbind(output_prop.inter, cur_proportion)
    }
  }
}

# Computing the network of co-occurrences for predators
output_cooc_per_site_pred <- NULL
output_cooccurrence_pred <- NULL

for(sp in unique(occupancy[occupancy$trophic_level=='consumer',]$species)){
  sp_sites <- unique(occupancy[occupancy$species==sp,]$site)
  for(sit in unique(sp_sites)){
    local_com <- occupancy[occupancy$site==sit,]
    cur_potential_prey <- unique(local_com[local_com$trophic_level=='resource',]$species)
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

# We now compute the network of co-occurrences for resources

output_cooc_per_site_prey <- NULL
output_cooccurrence_prey <- NULL

for(sp in unique(occupancy[occupancy$trophic_level=='resource',]$species)){
  sp_sites <- unique(occupancy[occupancy$species==sp,]$site)
  for(sit in unique(sp_sites)){
    local_com <- occupancy[occupancy$site==sit,]
    cur_potential_pred <- unique(local_com[local_com$trophic_level=='consumer',]$species)
    cur_out <- data.frame(site=sit, species=sp, potential_pred=cur_potential_pred)
    if(is.null(output_cooc_per_site_prey)){
      output_cooc_per_site_prey <- cur_out
    }else{
      output_cooc_per_site_prey<- rbind(output_cooc_per_site_prey, cur_out)
    }
  }
  # We calculate the total number of co-occurring resources
  total_potential <- length(unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred))
  potential_pred<- unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred)
  
  cur_out_fin <- data.frame(species=sp, total_pot_pred=total_potential, coocurring_species=potential_pred)
  if(is.null(output_cooccurrence_prey)){
    output_cooccurrence_prey <- cur_out_fin
  }else{
    output_cooccurrence_prey <- rbind(output_cooccurrence_prey, cur_out_fin)
  }
}

# We calculate the proportion of sites in which every pair of species co-occur
proportion_coocurrence <- NULL

for(sp in unique(output_cooc_per_site_pred$species)){
  total_n_sites <- length(unique(output_cooc_per_site_pred$site))
  cur_sp <- output_cooc_per_site_pred[output_cooc_per_site_pred$species==sp,]
  for(pot_prey in unique(cur_sp$potential_prey)){
    n_cooc <- dim(cur_sp[cur_sp$potential_prey==pot_prey,])[1]
    prop_cooc <- n_cooc#/total_n_sites
    
    cur_out <- data.frame(species=sp, potential_prey=pot_prey, coocurring_sites=n_cooc, prop.cooc=prop_cooc)
    
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

p <- 0.074

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

write.csv(realisations_indegree, file = 'random_realisations_expected_indegree_gottin-HP.csv')
write.csv(realisations_outdegree, file = 'random_realisations_expected_outdegree_gottin-HP.csv')

# Sum all the probabilities of a given predator to obtain the expected number of links

expected_int_model_indegree <- expected_int_model %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(expected_indegree = sum(prob_int))

expected_realised_indegree <- merge(number_interactions_pred, expected_int_model_indegree, by='species')

expected_realised_indegree_freq <- merge(expected_realised_indegree, interaction_proportion[,c('species','relative_N_alpha')], by='species')


write.csv(expected_realised_indegree_freq, file = 'expected_realised_indegree_gottin-HP.csv')

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

write.csv(expected_realised_outdegree_freq, file = 'expected_realised_outdegree_gottin-HP.csv')

# Here we have the occupancy distributions for both trophic levels.
occupancy_pred <- unique(number_sites_pred)

occupancy_prey <- unique(number_sites_prey)

write.csv(occupancy_pred, file = 'occupancy_pred_Gottin-HP.csv')
write.csv(occupancy_prey, file = 'occupancy_prey_Gottin-HP.csv')

# Here we have the outputs from both networks the realised and the co-occurrence network.
network_real <- unique(number_interactions)

network_real_pred <- unique(number_interactions_pred)

network_real_prey <- unique(number_interactions_prey)

network_coocurrence_pred <- unique(output_cooccurrence_pred[,1:2])

network_coocurrence_prey <- unique(output_cooccurrence_prey[,1:2])

write.csv(network_real_pred, file = 'network_real_pred_Gottin-HP.csv')
write.csv(network_real_prey, file = 'network_real_prey_Gottin-HP.csv')
write.csv(network_coocurrence_pred, file = 'network_coocurrence_pred_Gottin-HP.csv')
write.csv(network_coocurrence_prey, file = 'network_coocurrence_prey_Gottin-HP.csv')


# We now built a model to infer interactions without using frequency of co-occurrence to predict interactions.
# This corresponds to a random pruning of the network of co-occurrences.

p <- 0.21 # p here is the percentage of interactions that where realised from the co-occurrence network to the network of biotic interactions

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

write.csv(realisations_indegree_null_model, file = 'null_random_realisations_expected_indegree_gottin-HP.csv')
write.csv(realisations_outdegree_null_model, file = 'null_random_realisations_expected_outdegree_gottin-HP.csv')


############################## MUTUALISTIC NETWORK ####################################

## Calculating species occupancy MUT_NETS

occupancy <- NULL

for(p in unique(metadata$site)){
  local_com <- as.data.frame(mut_nets[[p]])
  
  species_present <- as.data.frame(names(local_com))
  species_present <- setNames(species_present, c("species"))
  species_present$trophic_level <- 'consumer'
  
  resources_present <- as.data.frame(rownames(local_com))
  resources_present <- setNames(resources_present, c("species"))
  resources_present$trophic_level <- 'resource'
  
  
  cur_data <- data.frame(site=p,rbind(species_present,resources_present))
  
  if(is.null(occupancy)){
    occupancy <- cur_data
  }else{
    occupancy <- rbind(occupancy, cur_data)
  }
}

number_sites <- NULL
for(sp in unique(occupancy$species)){
  cur_sp <- occupancy[occupancy$species==sp,]
  num_sites <- length(cur_sp$site)
  cur_data <- data.frame(species=sp,trophic_level=cur_sp$trophic_level,sites=num_sites)
  
  if(is.null(number_sites)){
    number_sites <- cur_data
  }else{
    number_sites <- rbind(number_sites, cur_data)
  }
}

number_sites <- unique(number_sites[number_sites$sites!=0,])

number_sites_pred <- number_sites[number_sites$trophic_level=='consumer',]
number_sites_prey <- number_sites[number_sites$trophic_level=='resource',]

## create metaweb for each interaction type
mut_raw$dummy = 1

mut_meta = frame2webs(mut_raw, varnames=c("P.Genus.Species","Genus.Species","dummy"))[[1]]

sum(mut_meta)

## Calculating total number of interactions per species (of the higher trophic level)
metaweb <- graph_from_incidence_matrix(mut_meta, directed=T, mode='out')

number_interactions_pred <- as.data.frame(igraph::degree(metaweb, V(metaweb)[which(V(metaweb)$type == TRUE)]$name, mode='in'))

number_interactions_pred <- tibble::rownames_to_column(number_interactions_pred, "species")

number_interactions_pred <- setNames(number_interactions_pred, c("species","interactions"))

number_interactions_prey <- as.data.frame(igraph::degree(metaweb, V(metaweb)[which(V(metaweb)$type == FALSE)]$name, mode='out'))

number_interactions_prey <- tibble::rownames_to_column(number_interactions_prey, "species")

number_interactions_prey <- setNames(number_interactions_prey, c("species","interactions"))

number_interactions <- rbind(number_interactions_pred, number_interactions_prey)


real_interactions <- as.data.frame(get.edgelist(metaweb))
real_interactions <- setNames(real_interactions, c("resource","consumer"))

site_interaction <- merge(number_interactions, number_sites, by="species") 


# We calculate the number of times two species itneract relative to the number of times they co-occur
output_prop.inter <- NULL
output_site.inter <- NULL

for(con in unique(real_interactions$consumer)){
  cur_sp <- real_interactions[real_interactions$consumer==con,]
  for(res in unique(cur_sp$resource)){
    consumer_sites <- unique(occupancy[occupancy$species==con,]$site)
    resource_sites <- unique(occupancy[occupancy$species==res,]$site)
    number_cooc <- sum(consumer_sites %in% resource_sites, na.rm=TRUE)
    for(p in unique(intersect(consumer_sites,resource_sites))){ 
      n <- as.data.frame(mut_nets[[p]])
      
      local_com <- igraph::graph_from_incidence_matrix(n, directed = T, mode='out')
      local_com <- igraph::delete.vertices(local_com, which(igraph::degree(local_com,  igraph::V(local_com)) == 0))
      
      cur_interactions <- as.data.frame(igraph::get.edgelist(local_com))
      cur_interactions <- setNames(cur_interactions, c("resource","consumer"))
      
      interactions <- dim(cur_interactions[cur_interactions$consumer==con & cur_interactions$resource==res,])[1]
      if(interactions>1){
        interactions <- 1
      } 
      
      cur_out <- data.frame(consumer=con, resource=res, site=p, interact=interactions)
      if(is.null(output_site.inter)){
        output_site.inter <- cur_out
      }else{
        output_site.inter <- rbind(output_site.inter, cur_out)
      }
    }
    
    sites_interact <- dim(output_site.inter[output_site.inter$consumer==con & output_site.inter$resource==res & output_site.inter$interact==1,])[1]
    
    cooccurring_sites <- dim(output_site.inter[output_site.inter$consumer==con & output_site.inter$resource==res,])[1]
    
    cur_proportion <- data.frame(consumer=con, resource=res, coocurring_sites=cooccurring_sites, interacting_sites=sites_interact, proportion_interact=sites_interact/cooccurring_sites)
    
    if(is.null(output_prop.inter)){
      output_prop.inter <- cur_proportion
    }else{
      output_prop.inter <- rbind(output_prop.inter, cur_proportion)
    }
  }
}

# Computing the network of co-occurrences for predators
output_cooc_per_site_pred <- NULL
output_cooccurrence_pred <- NULL

for(sp in unique(occupancy[occupancy$trophic_level=='consumer',]$species)){
  sp_sites <- unique(occupancy[occupancy$species==sp,]$site)
  for(sit in unique(sp_sites)){
    local_com <- occupancy[occupancy$site==sit,]
    cur_potential_prey <- unique(local_com[local_com$trophic_level=='resource',]$species)
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

# We now compute the network of co-occurrences for resources

output_cooc_per_site_prey <- NULL
output_cooccurrence_prey <- NULL

for(sp in unique(occupancy[occupancy$trophic_level=='resource',]$species)){
  sp_sites <- unique(occupancy[occupancy$species==sp,]$site)
  for(sit in unique(sp_sites)){
    local_com <- occupancy[occupancy$site==sit,]
    cur_potential_pred <- unique(local_com[local_com$trophic_level=='consumer',]$species)
    cur_out <- data.frame(site=sit, species=sp, potential_pred=cur_potential_pred)
    if(is.null(output_cooc_per_site_prey)){
      output_cooc_per_site_prey <- cur_out
    }else{
      output_cooc_per_site_prey<- rbind(output_cooc_per_site_prey, cur_out)
    }
  }
  # We calculate the total number of co-occurring resources
  total_potential <- length(unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred))
  potential_pred<- unique(output_cooc_per_site_prey[output_cooc_per_site_prey$species==sp,]$potential_pred)
  
  cur_out_fin <- data.frame(species=sp, total_pot_pred=total_potential, coocurring_species=potential_pred)
  if(is.null(output_cooccurrence_prey)){
    output_cooccurrence_prey <- cur_out_fin
  }else{
    output_cooccurrence_prey <- rbind(output_cooccurrence_prey, cur_out_fin)
  }
}

# We calculate the proportion of sites in which every pair of species co-occur
proportion_coocurrence <- NULL

for(sp in unique(output_cooc_per_site_pred$species)){
  total_n_sites <- length(unique(output_cooc_per_site_pred$site))
  cur_sp <- output_cooc_per_site_pred[output_cooc_per_site_pred$species==sp,]
  for(pot_prey in unique(cur_sp$potential_prey)){
    n_cooc <- dim(cur_sp[cur_sp$potential_prey==pot_prey,])[1]
    prop_cooc <- n_cooc#/total_n_sites
    
    cur_out <- data.frame(species=sp, potential_prey=pot_prey, coocurring_sites = n_cooc, prop.cooc=prop_cooc)
    
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

p <- 0.044

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
  
  while(any(rowSums(random_mat)==0)){
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

write.csv(realisations_indegree, file = 'random_realisations_expected_indegree_gottin-PP.csv')
write.csv(realisations_outdegree, file = 'random_realisations_expected_outdegree_gottin-PP.csv')

# Sum all the probabilities of a given predator to obtain the expected number of links

expected_int_model_indegree <- expected_int_model %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(expected_indegree = sum(prob_int))

expected_realised_indegree <- merge(number_interactions_pred, expected_int_model_indegree, by='species')

expected_realised_indegree_freq <- merge(expected_realised_indegree, interaction_proportion[,c('species','relative_N_alpha')], by='species')

colnames(expected_realised_indegree_freq)[7] <- 'mean_cooccurrence'

write.csv(expected_realised_indegree_freq, file = 'expected_realised_indegree_gottin-PP.csv')

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

write.csv(expected_realised_outdegree_freq, file = 'expected_realised_outdegree_gottin-PP.csv')

# Here we have the occupancy distributions for both trophic levels.
occupancy_pred <- unique(number_sites_pred)

occupancy_prey <- unique(number_sites_prey)

write.csv(occupancy_pred, file = 'occupancy_pred_Gottin-PP.csv')
write.csv(occupancy_prey, file = 'occupancy_prey_Gottin-PP.csv')

# Here we have the outputs from both networks the realised and the co-occurrence network.
network_real <- unique(number_interactions)

network_real_pred <- unique(number_interactions_pred)

network_real_prey <- unique(number_interactions_prey)

network_coocurrence_pred <- unique(output_cooccurrence_pred[,1:2])

network_coocurrence_prey <- unique(output_cooccurrence_prey[,1:2])

write.csv(network_real_pred, file = 'network_real_pred_Gottin-PP.csv')
write.csv(network_real_prey, file = 'network_real_prey_Gottin-PP.csv')
write.csv(network_coocurrence_pred, file = 'network_coocurrence_pred_Gottin-PP.csv')
write.csv(network_coocurrence_prey, file = 'network_coocurrence_prey_Gottin-PP.csv')


# We now built a model to infer interactions without using frequency of co-occurrence to predict interactions.
# This corresponds to a random pruning of the network of co-occurrences.

p <- 0.12 # p here is the percentage of interactions that where realised from the co-occurrence network to the network of biotic interactions

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


write.csv(realisations_indegree_null_model, file = 'null_random_realisations_expected_indegree_gottin-PP.csv')
write.csv(realisations_outdegree_null_model, file = 'null_random_realisations_expected_outdegree_gottin-PP.csv')




