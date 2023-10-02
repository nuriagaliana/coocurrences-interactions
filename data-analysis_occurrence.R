## ---------------------------
##
## Script name: data-analysis_occurrence.r
##
##
## Author: Dr Nuria Galiana
##
## Marie Curie Fellow at the National Museum 
## of Natural Sciences (CSIC)
##
## Date Created: 17-05-2022
##
## Email: galiana.nuria@gmail.com
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Galiana, et al. (2022) 
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(dplyr)
require(purrr)
require(sars)
require(ggplot2)


############################################## DATA UPLOAD ################################################

### DATA FOR THE OCCUPANCY DISTRIBUTION
output_garraf_hp_occupancy_pred <- read.table("./Garraf-Montseny-Olot/occupancy_pred_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_occupancy_prey <- read.table("./Garraf-Montseny-Olot/occupancy_prey_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_occupancy_pred <- read.table("./Garraf-Montseny-Olot/occupancy_pred_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_occupancy_prey <- read.table("./Garraf-Montseny-Olot/occupancy_prey_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_occupancy_pred <- read.table("./Garraf-Montseny-Olot/occupancy_pred_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_occupancy_prey <- read.table("./Garraf-Montseny-Olot/occupancy_prey_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_montseny_occupancy_pred <- read.table("./Garraf-Montseny-Olot/occupancy_pred_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_occupancy_prey <- read.table("./Garraf-Montseny-Olot/occupancy_prey_montseny.csv", sep=",", header=TRUE)[-1]
output_olot_occupancy_pred <- read.table("./Garraf-Montseny-Olot/occupancy_pred_olot.csv", sep=",", header=TRUE)[-1]
output_olot_occupancy_prey <- read.table("./Garraf-Montseny-Olot/occupancy_prey_olot.csv", sep=",", header=TRUE)[-1]
output_nahuel_occupancy_pred <- read.table("./Nahuel/occupancy_pred_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_occupancy_prey <- read.table("./Nahuel/occupancy_prey_nahuel.csv", sep=",", header=TRUE)[-1]
output_quercus_occupancy_pred <- read.table("./Quercus/occupancy_pred_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_occupancy_prey <- read.table("./Quercus/occupancy_prey_quercus.csv", sep=",", header=TRUE)[-1]
output_galpar_occupancy_pred <- read.table("./Salix/occupancy_pred_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_occupancy_prey <- read.table("./Salix/occupancy_prey_galpar.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_occupancy_pred <- read.table("./Gottin/occupancy_pred_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_occupancy_prey <- read.table("./Gottin/occupancy_prey_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_occupancy_pred <- read.table("./Gottin/occupancy_pred_Gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_occupancy_prey <- read.table("./Gottin/occupancy_prey_Gottin-PP.csv", sep=",", header=TRUE)[-1]


### DATA FOR NETWORKS OF BIOTIC INTERACTIONS
output_garraf_hp_real_pred <- read.table("./Garraf-Montseny-Olot/network_real_pred_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_real_prey <- read.table("./Garraf-Montseny-Olot/network_real_prey_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_real_pred <- read.table("./Garraf-Montseny-Olot/network_real_pred_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_real_prey <- read.table("./Garraf-Montseny-Olot/network_real_prey_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_real_pred <- read.table("./Garraf-Montseny-Olot/network_real_pred_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_real_prey <- read.table("./Garraf-Montseny-Olot/network_real_prey_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_montseny_real_pred <- read.table("./Garraf-Montseny-Olot/network_real_pred_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_real_prey <- read.table("./Garraf-Montseny-Olot/network_real_prey_montseny.csv", sep=",", header=TRUE)[-1]
output_olot_real_pred <- read.table("./Garraf-Montseny-Olot/network_real_pred_olot.csv", sep=",", header=TRUE)[-1]
output_olot_real_prey <- read.table("./Garraf-Montseny-Olot/network_real_prey_olot.csv", sep=",", header=TRUE)[-1]
output_nahuel_real_pred <- read.table("./Nahuel/network_real_pred_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_real_prey <- read.table("./Nahuel/network_real_prey_nahuel.csv", sep=",", header=TRUE)[-1]
output_quercus_real_pred <- read.table("./Quercus/network_real_pred_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_real_prey <- read.table("./Quercus/network_real_prey_quercus.csv", sep=",", header=TRUE)[-1]
output_galpar_real_pred <- read.table("./Salix/network_real_pred_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_real_prey <- read.table("./Salix/network_real_prey_galpar.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_real_pred <- read.table("./Gottin/network_real_pred_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_real_prey <- read.table("./Gottin/network_real_prey_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_real_pred <- read.table("./Gottin/network_real_pred_Gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_real_prey <- read.table("./Gottin/network_real_prey_Gottin-PP.csv", sep=",", header=TRUE)[-1]


### DATA FOR CO-OCCURRENCE NETWORKS
output_garraf_hp_cooc_pred <- read.table("./Garraf-Montseny-Olot/network_coocurrence_pred_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_cooc_prey <- read.table("./Garraf-Montseny-Olot/network_coocurrence_prey_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_cooc_pred <- read.table("./Garraf-Montseny-Olot/network_coocurrence_pred_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_cooc_prey <- read.table("./Garraf-Montseny-Olot/network_coocurrence_prey_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_cooc_pred <- read.table("./Garraf-Montseny-Olot/network_coocurrence_pred_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_cooc_prey <- read.table("./Garraf-Montseny-Olot/network_coocurrence_prey_garraf-pp2.csv", sep=",", header=TRUE)[-1]
output_montseny_cooc_pred <- read.table("./Garraf-Montseny-Olot/network_coocurrence_pred_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_cooc_prey <- read.table("./Garraf-Montseny-Olot/network_coocurrence_prey_montseny.csv", sep=",", header=TRUE)[-1]
output_olot_cooc_pred <- read.table("./Garraf-Montseny-Olot/network_coocurrence_pred_olot.csv", sep=",", header=TRUE)[-1]
output_olot_cooc_prey <- read.table("./Garraf-Montseny-Olot/network_coocurrence_prey_olot.csv", sep=",", header=TRUE)[-1]
output_nahuel_cooc_pred <- read.table("./Nahuel/network_coocurrence_pred_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_cooc_prey <- read.table("./Nahuel/network_coocurrence_prey_nahuel.csv", sep=",", header=TRUE)[-1]
output_quercus_cooc_pred <- read.table("./Quercus/network_coocurrence_pred_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_cooc_prey <- read.table("./Quercus/network_coocurrence_prey_quercus.csv", sep=",", header=TRUE)[-1]
output_galpar_cooc_pred <- read.table("./Salix/network_coocurrence_pred_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_cooc_prey <- read.table("./Salix/network_coocurrence_prey_galpar.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_cooc_pred <- read.table("./Gottin/network_coocurrence_pred_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_cooc_prey <- read.table("./Gottin/network_coocurrence_prey_Gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_cooc_pred <- read.table("./Gottin/network_coocurrence_pred_Gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_cooc_prey <- read.table("./Gottin/network_coocurrence_prey_Gottin-PP.csv", sep=",", header=TRUE)[-1]


### DATA FOR EXPECTED INTERACTIONS BASED ON FREQUENCY OF CO-OCCURRENCE
output_garraf_hp_expected_pred <- read.table("./Garraf-Montseny-Olot/expected_realised_indegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_expected_prey <- read.table("./Garraf-Montseny-Olot/expected_realised_outdegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_expected_pred <- read.table("./Garraf-Montseny-Olot/expected_realised_indegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_expected_prey <- read.table("./Garraf-Montseny-Olot/expected_realised_outdegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_expected_pred <- read.table("./Garraf-Montseny-Olot/expected_realised_indegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_expected_prey <- read.table("./Garraf-Montseny-Olot/expected_realised_outdegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_montseny_expected_pred <- read.table("./Garraf-Montseny-Olot/expected_realised_indegree_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_expected_prey <- read.table("./Garraf-Montseny-Olot/expected_realised_outdegree_montseny.csv", sep=",", header=TRUE)[-1]
output_olot_expected_pred <- read.table("./Garraf-Montseny-Olot/expected_realised_indegree_olot.csv", sep=",", header=TRUE)[-1]
output_olot_expected_prey <- read.table("./Garraf-Montseny-Olot/expected_realised_outdegree_olot.csv", sep=",", header=TRUE)[-1]
output_nahuel_expected_pred <- read.table("./Nahuel/expected_realised_indegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_expected_prey <- read.table("./Nahuel/expected_realised_outdegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_quercus_expected_pred <- read.table("./Quercus/expected_realised_indegree_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_expected_prey <- read.table("./Quercus/expected_realised_outdegree_quercus.csv", sep=",", header=TRUE)[-1]
output_galpar_expected_pred <- read.table("./Salix/expected_realised_indegree_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_expected_prey <- read.table("./Salix/expected_realised_outdegree_galpar.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_expected_pred <- read.table("./Gottin/expected_realised_indegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_expected_prey <- read.table("./Gottin/expected_realised_outdegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_expected_pred <- read.table("./Gottin/expected_realised_indegree_gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_expected_prey <- read.table("./Gottin/expected_realised_outdegree_gottin-PP.csv", sep=",", header=TRUE)[-1]


### DATA FOR THE RANDOMIZATION FOR EXPECTED INTERACTIONS BASED ON FREQUENCY OF CO-OCCURRENCE
output_garraf_hp_random_expected_pred <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_indegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_random_expected_prey <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_outdegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_random_expected_pred <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_indegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_random_expected_prey <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_outdegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_random_expected_pred <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_indegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_random_expected_prey <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_outdegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_olot_random_expected_pred <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_indegree_olot.csv", sep=",", header=TRUE)[-1]
output_olot_random_expected_prey <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_outdegree_olot.csv", sep=",", header=TRUE)[-1]
output_montseny_random_expected_pred <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_indegree_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_random_expected_prey <- read.table("./Garraf-Montseny-Olot/random_realisations_expected_outdegree_montseny.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_random_expected_pred <- read.table("./Gottin/random_realisations_expected_indegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_random_expected_prey <- read.table("./Gottin/random_realisations_expected_outdegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_random_expected_pred <- read.table("./Gottin/random_realisations_expected_indegree_gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_random_expected_prey <- read.table("./Gottin/random_realisations_expected_outdegree_gottin-PP.csv", sep=",", header=TRUE)[-1]
output_nahuel_random_expected_pred <- read.table("./Nahuel/random_realisations_expected_indegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_random_expected_prey <- read.table("./Nahuel/random_realisations_expected_outdegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_galpar_random_expected_pred <- read.table("./Salix/random_realisations_expected_indegree_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_random_expected_prey <- read.table("./Salix/random_realisations_expected_outdegree_galpar.csv", sep=",", header=TRUE)[-1]
output_quercus_random_expected_pred <- read.table("./Quercus/random_realisations_expected_indegree_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_random_expected_prey <- read.table("./Quercus/random_realisations_expected_outdegree_quercus.csv", sep=",", header=TRUE)[-1]

### DATA FOR RANDOM PRUNING WHERE FREQUENCY OF CO-OCCURRENCE IS NOT CONSIDERED (NULL MODEL)
output_garraf_hp_null_random_expected_pred <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_indegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_hp_null_random_expected_prey <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_outdegree_garraf-HP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_null_random_expected_pred <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_indegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_null_random_expected_prey <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_outdegree_garraf-PP.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_null_random_expected_pred <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_indegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_null_random_expected_prey <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_outdegree_garraf-PP2.csv", sep=",", header=TRUE)[-1]
output_olot_null_random_expected_pred <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_indegree_olot.csv", sep=",", header=TRUE)[-1]
output_olot_null_random_expected_prey <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_outdegree_olot.csv", sep=",", header=TRUE)[-1]
output_montseny_null_random_expected_pred <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_indegree_montseny.csv", sep=",", header=TRUE)[-1]
output_montseny_null_random_expected_prey <- read.table("./Garraf-Montseny-Olot/null_random_realisations_expected_outdegree_montseny.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_null_random_expected_pred <- read.table("./Gottin/null_random_realisations_expected_indegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_hp_null_random_expected_prey <- read.table("./Gottin/null_random_realisations_expected_outdegree_gottin-HP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_null_random_expected_pred <- read.table("./Gottin/null_random_realisations_expected_indegree_gottin-PP.csv", sep=",", header=TRUE)[-1]
output_gottin_pp_null_random_expected_prey <- read.table("./Gottin/null_random_realisations_expected_outdegree_gottin-PP.csv", sep=",", header=TRUE)[-1]
output_nahuel_null_random_expected_pred <- read.table("./Nahuel/null_random_realisations_expected_indegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_nahuel_null_random_expected_prey <- read.table("./Nahuel/null_random_realisations_expected_outdegree_nahuel.csv", sep=",", header=TRUE)[-1]
output_galpar_null_random_expected_pred <- read.table("./Salix/null_random_realisations_expected_indegree_galpar.csv", sep=",", header=TRUE)[-1]
output_galpar_null_random_expected_prey <- read.table("./Salix/null_random_realisations_expected_outdegree_galpar.csv", sep=",", header=TRUE)[-1]
output_quercus_null_random_expected_pred <- read.table("./Quercus/null_random_realisations_expected_indegree_quercus.csv", sep=",", header=TRUE)[-1]
output_quercus_null_random_expected_prey <- read.table("./Quercus/null_random_realisations_expected_outdegree_quercus.csv", sep=",", header=TRUE)[-1]

###########################################################################################################################



############################################## FIGURE 2 ##################################################################

### Comparing degree distribution of co-occurrence versus realised networks of interactions
### Change the name of the dataset according to the desired plot

df <- unique(output_galpar_real_pred) # Data corresponding to the network of biotic interactions
occur_pot_prey = as.vector(table(df$interactions))
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$interactions)))

temp_real <- data.frame(x, y)

df <- unique(output_galpar_cooc_pred) # Data corresponding to the network of co-occurrences
occur_pot_prey = as.vector(table(df$total_pot_prey))
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$total_pot_prey)))

temp_cooc <- data.frame(x, y)

# Functions that we want to fit and plot to the degree distribution of the network of biotic interactions
exponential <- nls(y ~ (exp(-x/b)), data = temp_real, start = list(b = 2), control=nls.control(maxiter = 1000))
power_law <- nls(y ~ (x^-a), data = temp_real, start = list(a = .01), control=nls.control(maxiter = 1e3))
truncated <- nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp_real, start = list(a = .0001, b = 2), control=nls.control(maxiter = 1e3))

# Plots in Figure 2 (and Figure S1 for resources)
ggplot(temp_real, aes(x=(x), y=(y))) + 
  geom_line(data=temp_real, aes(x=(x), y=(y)), color='darkgoldenrod3', size=1.5) +
  geom_point(color='darkgoldenrod3', size=3) +  
  geom_line(data=temp_real,aes(x=x, y=predict(truncated)), size=1.5,linetype=2) + # Change the function to match the best fit for each dataset
  geom_line(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=1.5) +
  geom_point(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  xlab('Number of links') + ylab('Cum. Probability') + ggtitle('Galpar consumers')

### Calculating the proportion of links kept from the co-occurrence network in the network of biotic interactions.
### Change the name of the dataset accordingly

df_real <- unique(output_quercus_real_pred)
df_cooc <- unique(output_quercus_cooc_pred)

proportion_realised <- sum(df_real$interactions)/sum(df_cooc$total_pot_prey)

############################################ END FIGURE 2 ################################################################



############ FITTING MODELS TO DEGREE DISTRIBUTIONS (Supplementary Table S2) ###############################

require(AICcmodavg)

## NETWORK OF BIOTIC INTERACTIONS

degree_dist_params <- NULL

failed <- tryCatch({
  mod1 <- tryCatch({
    nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp_real, start = list(a = .0001, b = 2), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod2 <- tryCatch({
    nls(y ~ (exp(-x/b)), data = temp_real, start = list(b = 2), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod3 <- tryCatch({
    nls(y ~ (x^-a), data = temp_real, start = list(a = .01), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod4 <- tryCatch({
    nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp_real, start = list(a = .3, b = .3), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  FALSE
}, warning = function(w) {
  FALSE
}, error = function(e) {
  TRUE
}, finally = {
  
})

if(!failed){
  model_list <- list(mod1,mod2,mod3,mod4)
  names_list <- c('mod1','mod2','mod3','mod4')
  if(length(which(is.na(model_list))) != 0){
    names_list <- names_list[-which(is.na(model_list))]
    model_list <- model_list[-which(is.na(model_list))]
  }
  names(model_list) <- names_list
  
  if(length(model_list) == 0) next
  
  if(length(model_list) == 1){
    model_name <- names_list[1]
    model <- summary(eval(as.symbol(model_name)))
  }else{
    aic_comp <- aictab(model_list)
    model_name <- tolower(as.character(aic_comp$Modnames[1]))
    model <- summary(eval(as.symbol(model_name)))
  }
  
  if(model_name == 'mod1' | model_name == 'mod4'){
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
  }else if(model_name == 'mod2'){
    cur_out <- data.frame( model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
  }else{
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
  }
  
  if(is.null(degree_dist_params)){
    degree_dist_params <- cur_out
  }else{
    degree_dist_params <- rbind(degree_dist_params, cur_out)
  }
}


# CO-OCCURRENCE NETWORK

degree_dist_params_coocur <- NULL

failed <- tryCatch({
  mod1 <- tryCatch({
    nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp_cooc, start = list(a = .0001, b = 2), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod2 <- tryCatch({
    nls(y ~ (exp(-x/b)), data = temp_cooc, start = list(b = 2), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod3 <- tryCatch({
    nls(y ~ (x^-a), data = temp_cooc, start = list(a = .01), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod4 <- tryCatch({
    nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp_cooc, start = list(a = .3, b = .3), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  FALSE
}, warning = function(w) {
  FALSE
}, error = function(e) {
  TRUE
}, finally = {
  
})

if(!failed){
  model_list <- list(mod1,mod2,mod3,mod4)
  names_list <- c('mod1','mod2','mod3','mod4')
  if(length(which(is.na(model_list))) != 0){
    names_list <- names_list[-which(is.na(model_list))]
    model_list <- model_list[-which(is.na(model_list))]
  }
  names(model_list) <- names_list
  
  if(length(model_list) == 0) next
  
  if(length(model_list) == 1){
    model_name <- names_list[1]
    model <- summary(eval(as.symbol(model_name)))
  }else{
    aic_comp <- aictab(model_list)
    model_name <- tolower(as.character(aic_comp$Modnames[1]))
    model <- summary(eval(as.symbol(model_name)))
  }
  
  if(model_name == 'mod1' | model_name == 'mod4'){
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
  }else if(model_name == 'mod2'){
    cur_out <- data.frame( model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
  }else{
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
  }
  
  if(is.null(degree_dist_params_coocur)){
    degree_dist_params_coocur <- cur_out
  }else{
    degree_dist_params_coocur <- rbind(degree_dist_params_coocur, cur_out)
  }
}

######################################## END MODEL FITTING #################################################



########################################### FIGURE 3 #######################################################

### Relationship between the number of potential interactions (based on co-occurrence) and the number of 
### realised interactions (biotic interactions) and our model predictions.

### Change the name of the dataset accordingly

df_real <- unique(output_galpar_real_pred) # data on biotic interactions
df_cooc <- unique(output_galpar_cooc_pred) # data on co-occurrences
df_expected <- unique(output_galpar_expected_pred[,c('species','expected_indegree')]) # model predictions. Change 'expected_indegree' to 'expected_outdegree' to analyse resources.

df_all <- merge(df_real, df_cooc, by= c("species")) 
df_all_expectation <- merge(df_all, df_expected, by= c("species")) 

# Plots in Figure 3 (and Figure S2 for resources)
ggplot(df_all_expectation, aes(x=total_pot_prey, y=interactions)) + # Change 'x=total_pot_prey' to 'x=total_pot_pred' to analyse resources throughout the plot
  geom_point(size=3)+ #geom_smooth()+
  geom_point(data=df_all_expectation, aes(x=total_pot_prey, y=expected_indegree), color='cyan4', size=4) + #Change 'expected_indegree' to 'expected_outdegree' to analyse resources throughout the plot
  #stat_smooth(fullrange = F, level = 0.8, color='cyan4',size=1.5)+
  geom_smooth(data=df_all_expectation, aes(x=total_pot_prey, y=expected_indegree), color='cyan4', size=1.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) + 
  theme(plot.margin = margin(t=5.5, r=10, b=5.5, l=5.5, "pt")) +
  xlab("Potential interactions") + ylab("Realised interactions") + ggtitle('Garraf PP consumers')

############################################# END FIGURE 3 ###################################################



############################################# FIGURE 5 #######################################################

### Relationship between number of potential interactions of species and their average frequency of 
### co-occurrence with their potential interacting partners in each dataset.

### Change the name of the dataset accordingly

df_freq_potential <- merge(unique(output_garraf_pp_cooc_pred), unique(output_garraf_pp_expected_pred[,c('species','mean_cooccurrence')]), by= c("species")) 

#Plots in Figure 4 (and Figure S3 for resources)
ggplot(df_freq_potential, aes(x=total_pot_prey, y=mean_cooccurrence)) + # Change'x=total_pot_prey' to 'x=total_pot_pred' to analyse the resource pattern
  geom_point(size=3)+ geom_smooth(method = 'loess')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) + 
  theme(plot.margin = margin(t=5.5, r=10, b=5.5, l=5.5, "pt")) +
  xlab("Potential interactions") + ylab("Mean frequency of co-occurrence") + ggtitle('Garraf PP consumers')

############################################# END FIGURE 5 ###################################################



########################################### FIGURE 6 #######################################################
### Comparing degree distribution of the pruned networks based on frequency of co-occurrences using multiple random 
### realisations (our theoretical model), network of co-occurrences and network of biotic interactions.

### Change the name of the dataset accordingly
df <- unique(output_galpar_real_pred) # Data corresponding to the network of biotic interactions
occur_pot_prey = as.vector(table(df$interactions))
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$interactions)))

temp_real <- data.frame(x, y)

df <- unique(output_galpar_cooc_pred) # Data corresponding to the network of co-occurrences
occur_pot_prey = as.vector(table(df$total_pot_prey)) # Change 'df$total_pot_prey' to 'df$total_pot_pred' to analyse resources
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$total_pot_prey))) # Change 'df$total_pot_prey' to 'df$total_pot_pred' to analyse resources

temp_cooc <- data.frame(x, y)

expected_degree<-NULL

for(r in unique(output_galpar_random_expected_pred$rep)){
  df <- unique(output_galpar_random_expected_pred[output_galpar_random_expected_pred$rep==r,]) 
  occur_pot_prey = as.vector(table(df$indegree)) # Change df$indegree to df$outdegree to analyse resources
  p = occur_pot_prey/sum(occur_pot_prey)
  y = rev(cumsum(rev(p)))
  x = as.numeric(names(table(df$indegree))) # Change df$indegree to df$outdegree to analyse resources
  
  temp_rando <- data.frame(x, y, replicate=r)
  
  if(is.null(expected_degree)){
    expected_degree <- temp_rando
  }else{
    expected_degree <- rbind(expected_degree, temp_rando)
  } 
}

# Plots in Figure 5 (and figure S5 for resources)
ggplot(expected_degree[expected_degree$x>0,], aes(x=(x), y=(y)),color=as.factor(replicate)) + 
  geom_line(data=expected_degree[expected_degree$x>0,], aes(x=(x), y=(y),color=as.factor(replicate)), size=1.5, alpha=0.1) +
  geom_point(data=expected_degree[expected_degree$x>0,], aes(x=(x), y=(y), color=as.factor(replicate)), size=3, alpha=0.1) +
  geom_line(data=temp_real, aes(x=(x), y=(y)), color='darkgoldenrod3', size=2) +
  geom_point(data=temp_real, aes(x=(x), y=(y)), color='darkgoldenrod3', size=4) +  
  geom_line(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=1.5) +
  #stat_smooth(data=temp_cooc, fullrange = F, level = 0.7, color='violetred3',size=1.5) + 
  geom_point(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  xlab('Number of links') + ylab('Cum. Probability') + ggtitle('Galpar resources')

############################################### END FIGURE 6 #############################################################


########################################## FIGURE S4 #######################################################
### Expected (based on our model) versus realised interactions (biotic interactions)

# For consumers
output_garraf_hp_expected_pred$dataset <- 'Garraf HP'
output_garraf_pp_expected_pred$dataset <- 'Garraf PP'
output_garraf_pp2_expected_pred$dataset <- 'Garraf PP2'
output_olot_expected_pred$dataset <- 'Olot'
output_montseny_expected_pred$dataset <- 'Montseny'
output_gottin_pp_expected_pred$dataset <- 'Gottin PP'
output_gottin_hp_expected_pred$dataset <- 'Gottin HP'
output_galpar_expected_pred$dataset <- 'Galpar'
output_nahuel_expected_pred$dataset <- 'Nahuel'
output_quercus_expected_pred$dataset <- 'Quercus'

output_all_expected_pred <- rbind(output_quercus_expected_pred,output_nahuel_expected_pred,output_galpar_expected_pred,
                                  output_gottin_hp_expected_pred,output_gottin_pp_expected_pred,output_montseny_expected_pred,
                                  output_olot_expected_pred,output_garraf_pp_expected_pred,output_garraf_pp2_expected_pred,output_garraf_hp_expected_pred)

# Plot in Figure S4a
ggplot(data=unique(output_all_expected_pred[,-c(3:4)]), aes(x=(expected_indegree), y=(interactions),color=dataset)) + 
  geom_point(size=3, alpha=0.6)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) + 
  geom_abline(slope=1)+
  theme(plot.margin = margin(t=5.5, r=10, b=5.5, l=5.5, "pt")) +
  xlab("Expected indegree") + ylab("Realised indegree")

#For resources
output_garraf_hp_expected_prey$dataset <- 'Garraf HP'
output_garraf_pp_expected_prey$dataset <- 'Garraf PP'
output_garraf_pp2_expected_prey$dataset <- 'Garraf PP2'
output_olot_expected_prey$dataset <- 'Olot'
output_montseny_expected_prey$dataset <- 'Montseny'
output_gottin_pp_expected_prey$dataset <- 'Gottin PP'
output_gottin_hp_expected_prey$dataset <- 'Gottin HP'
output_galpar_expected_prey$dataset <- 'Galpar'
output_nahuel_expected_prey$dataset <- 'Nahuel'
output_quercus_expected_prey$dataset <- 'Quercus'

output_all_expected_prey <- rbind(output_quercus_expected_prey,output_nahuel_expected_prey,output_galpar_expected_prey,
                                  output_gottin_hp_expected_prey,output_gottin_pp_expected_prey,output_montseny_expected_prey,
                                  output_olot_expected_prey,output_garraf_pp_expected_prey,output_garraf_pp2_expected_prey,output_garraf_hp_expected_prey)

# Plot in Figure S4b
ggplot(data=unique(output_all_expected_prey[,-3]), aes(x=(expected_outdegree), y=(interactions),color=dataset)) + 
  geom_point(size=3, alpha=0.6)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) + 
  geom_abline(slope=1)+
  theme(plot.margin = margin(t=5.5, r=10, b=5.5, l=5.5, "pt")) +
  xlab("Expected outdegree") + ylab("Realised outdegree")

##################################################### END FIGURE S4 #######################################################



########################################### FIGURE S6 #######################################################
### Comparing degree distribution of the pruned networks removing links randomly (null model), network of co-occurrences and network of biotic interactions.

### Change the name of the dataset accordingly

df <- unique(output_galpar_real_prey) # Data corresponding to the network of biotic interactions
occur_pot_prey = as.vector(table(df$interactions))
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$interactions)))

temp_real <- data.frame(x, y)

df <- unique(output_galpar_cooc_prey) # Data corresponding to the network of co-occurrences
occur_pot_prey = as.vector(table(df$total_pot_pred)) # Change 'df$total_pot_prey' to 'df$total_pot_pred' to analyse resources
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$total_pot_pred))) # Change df$total_pot_pred for consumers and df$total_pot_prey

temp_cooc <- data.frame(x, y)


expected_degree<-NULL # Data corresponding to the pruned network randomly

for(r in unique(output_galpar_null_random_expected_prey$rep)){ # Change to appropiate dataset
  df <- unique(output_galpar_null_random_expected_prey[output_galpar_null_random_expected_prey$rep==r,]) # Change to appropiate dataset
  occur_pot_prey = as.vector(table(df$outdegree)) # Change df$indegree to df$outdegree to analyse resources
  p = occur_pot_prey/sum(occur_pot_prey)
  y = rev(cumsum(rev(p)))
  x = as.numeric(names(table(df$outdegree))) # Change df$indegree to df$outdegree to analyse resources
  
  temp_null <- data.frame(x, y, replicate=r)
  
  if(is.null(expected_degree)){
    expected_degree <- temp_null
  }else{
    expected_degree <- rbind(expected_degree, temp_null)
  } 
}

# Plots in Figure S6 (null model)
ggplot(expected_degree[expected_degree$x>0,], aes(x=(x), y=(y)),color=as.factor(replicate)) + 
  geom_line(data=expected_degree[expected_degree$x>0,], aes(x=(x), y=(y),color=as.factor(replicate)), size=1.5, alpha=0.1) +
  geom_point(data=expected_degree[expected_degree$x>0,], aes(x=(x), y=(y), color=as.factor(replicate)), size=3, alpha=0.1) +
  geom_line(data=temp_real[temp_real$x>0,], aes(x=(x), y=(y)), color='darkgoldenrod3', size=2) +
  geom_point(data=temp_real[temp_real$x>0,], aes(x=(x), y=(y)), color='darkgoldenrod3', size=4) +  
  geom_line(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=1.5) +
  #stat_smooth(data=temp_cooc, fullrange = F, level = 0.7, color='violetred3',size=1.5) + 
  geom_point(data=temp_cooc, aes(x=(x), y=(y)), color='violetred3', size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 24 ,face="bold"), plot.title = element_text(size=20),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", color = 'black', size = 1.7)) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  xlab('Number of links') + ylab('Cum. Probability') + ggtitle('Galpar resource')

############################################### END FIGURE S6 #############################################################


