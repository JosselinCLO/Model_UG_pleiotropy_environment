## This code has been used to generate the results of the manuscript
## Title : How environment and genetic architecture of unreduced gametes shape the establishment of autopolyploids
## Authors : Chen, Yu; Schmickl, Roswitha; Kolar, Filip; and Clo, Josselin.
## Contacts : josselin.clo@gmail.com; chengyu@natur.cuni.cz

### Functions to perform the simulations:


# A function to introduce mutations within newly formed gametes for the quantitative trait

Mutation<-function(haplotype, U, Locus_trait,Locus_UG, var.add.eff)
{
  
  Nb_mut_trait<-rpois(1,U) ## The number of mutations to be done for trait
  Nb_mut_UG_F<-rpois(1,U) ## The number of mutations to be done for UG female
  Nb_mut_UG_M<-rpois(1,U) ## The number of mutations to be done for UG male
  
  ## The positions of the mutations on the haploid genome
  
  position_mut<-c(sample(1:Locus_trait, Nb_mut_trait, replace = F), sample((Locus_trait+1):(Locus_trait+Locus_UG), Nb_mut_UG_F, replace = F),sample((Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG), Nb_mut_UG_M, replace = F)) 
  
  # Make the mutations
  
  if((Nb_mut_trait+Nb_mut_UG_F+Nb_mut_UG_M) != 0){
    
    ## If you have at least 1 mutation, do the mutations
    
    for(i in 1:length(position_mut)){
      
      haplotype[position_mut[i]]= as.numeric(haplotype[position_mut[i]]) + rnorm(n = 1, mean = 0, sd = var.add.eff^0.5)
      
    }
    
  }
  
  return(haplotype)
  
}

# A function to select the reproducers 

reproduction <- function(Diploid_ind, Tetraploid_ind, Locus_trait,Locus_UG, U, Npop, fitness, step, env_effect, var.add.eff, pleiotropy, pleio_type){
  
  # The probability of being selected for reproduction depends on the fitness of individuals
  
  vector_proba_sampling = fitness
  
  # The probabilities of producing UG for female and male functions
  
  vector_proba_unred_gam_female = proba_unreduced_gametes_female(Locus_trait,Locus_UG,Npop,Diploid_ind, Tetraploid_ind)
  vector_proba_unred_gam_male = proba_unreduced_gametes_male(Locus_trait,Locus_UG,Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type)
  
  # During the adaptation process (step == 2), environment can influence UG production
  
  if(step == 2){
    
    effect_env_UG = c(rep(env_effect, time(length(vector_proba_unred_gam_male))))
    
    vector_proba_unred_gam_female = vector_proba_unred_gam_female + effect_env_UG
    vector_proba_unred_gam_male = vector_proba_unred_gam_male + effect_env_UG
    
  }
  
  ## Reproduction
  
  Nb_offsprings = 0
  
  diplo_ind_temp = c(NULL)
  tetra_ind_temp = c(NULL)
  
  repeat{
    
    # Select the first parent (female by default)
    
    Parent_1 = sample(1:Npop, 1, replace = F, prob = vector_proba_sampling)
    if(Parent_1 <= (nrow(Diploid_ind)/2)){Ploidy_parent1 = "diploid"}
    else{Ploidy_parent1 = "tetraploid"}
    
    # Select the kind of gamete (reduced or unreduced)
    
    if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam_female[Parent_1] <= runif(1, 0, 1)){Gamete_parent1 = "haplo"}
    else if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam_female[Parent_1] >= runif(1, 0, 1)){Gamete_parent1 = "diplo"}
    else if(Ploidy_parent1 == "tetraploid" && vector_proba_unred_gam_female[Parent_1] <= runif(1, 0, 1)){Gamete_parent1 = "diplo"}
    else{Gamete_parent1 = "tetra"}
    
    # We assume fully outcrossing populations
    
    repeat{Parent_2 = sample(1:Npop, 1, replace = F, prob = vector_proba_sampling)
      if(Parent_2 != Parent_1){break}
    }
      
    if(Parent_2 <= (nrow(Diploid_ind)/2)){Ploidy_parent2 = "diploid"}
    else{Ploidy_parent2 = "tetraploid"}
      
    if(Ploidy_parent2 == "diploid" && vector_proba_unred_gam_male[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "haplo"}
    else if(Ploidy_parent2 == "diploid" && vector_proba_unred_gam_male[Parent_2] >= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
    else if(Ploidy_parent2 == "tetraploid" && vector_proba_unred_gam_male[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
    else{Gamete_parent2 = "tetra"}
    
    # We assume no triploids 
    
    if(Ploidy_parent1 == "diploid" && Ploidy_parent2 == "diploid" && Gamete_parent1 == "haplo" && Gamete_parent2 == "haplo"){
      
      # First case, two diploids parents producing reduced gametes (diploid offsprings)
      
      Nb_offsprings = Nb_offsprings + 1
      
      # Gametes 1 & 2
      
      Gamete_1_temp = c(NULL)
      Gamete_2_temp = c(NULL)
      
      for(l in 1:(Locus_trait+Locus_UG+Locus_UG)){
        
        temp_ind_1 = c(Diploid_ind[2*Parent_1,l],Diploid_ind[(2*Parent_1)-1,l])
        temp_ind_2 = c(Diploid_ind[2*Parent_2,l],Diploid_ind[(2*Parent_2)-1,l])
        
        Gamete_1_temp[l] = sample(temp_ind_1, 1, replace = F)
        Gamete_2_temp[l] = sample(temp_ind_2, 1, replace = F)
        
      }
      
      Gamete_1_mut = Mutation(Gamete_1_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_2_mut = Mutation(Gamete_2_temp, U, Locus_trait,Locus_UG, var.add.eff)
      
      ploidy = c(rep("Diploid", times = 2))
      
      Offspirng_temp_fin = rbind(Gamete_1_mut,Gamete_2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(diplo_ind_temp) == 0){diplo_ind_temp = Offspirng_fin}
      else{diplo_ind_temp = rbind(diplo_ind_temp,Offspirng_fin)}
    }
    else if(Ploidy_parent1 == "diploid" && Ploidy_parent2 == "diploid" && Gamete_parent1 == "diplo" && Gamete_parent2 == "diplo"){
      
      # Second case, two diploids parents making unreduced gametes (tetraploid offsprings)
      
      Nb_offsprings = Nb_offsprings + 1
      
      # Gametes 1 & 2
      
      Gamete_1.1_temp = Diploid_ind[((2*Parent_1)-1),1:(Locus_trait+Locus_UG+Locus_UG)]
      Gamete_1.2_temp = Diploid_ind[(2*Parent_1),1:(Locus_trait+Locus_UG+Locus_UG)]
      Gamete_2.1_temp = Diploid_ind[((2*Parent_2)-1),1:(Locus_trait+Locus_UG+Locus_UG)]
      Gamete_2.2_temp = Diploid_ind[(2*Parent_2),1:(Locus_trait+Locus_UG+Locus_UG)]
      
      Gamete_1.1_mut = Mutation(Gamete_1.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_1.2_mut = Mutation(Gamete_1.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_2.1_mut = Mutation(Gamete_2.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_2.2_mut = Mutation(Gamete_2.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
      
      ploidy = c(rep("Tetraploid", times = 4))
      
      Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
      else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
    }
    else if(Ploidy_parent1 != Ploidy_parent2 && Gamete_parent1 == Gamete_parent2){
      
      # Third case, one diploid parent making unreduced gametes and one tetraploid parent making reduced gametes (tetraploid offsprings)
      
      Nb_offsprings = Nb_offsprings + 1
      
      if(Ploidy_parent1 == "diploid"){
        
        Gamete_1.1_temp = Diploid_ind[((2*Parent_1)-1),1:(Locus_trait+Locus_UG+Locus_UG)]
        Gamete_1.2_temp = Diploid_ind[(2*Parent_1),1:(Locus_trait+Locus_UG+Locus_UG)]
        
        Gamete_2.1_temp = c(NULL)
        Gamete_2.2_temp = c(NULL)
        
        for(l in 1:(Locus_trait+Locus_UG+Locus_UG)){
          
          temp_ind_2 = c(Tetraploid_ind[4*(Parent_2-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-1,l])
          
          Gamete_2_temp = c(sample(temp_ind_2, 2, replace = F))
          
          Gamete_2.1_temp[l] = Gamete_2_temp[1]
          Gamete_2.2_temp[l] = Gamete_2_temp[2]
        }
        
        Gamete_1.1_mut = Mutation(Gamete_1.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_1.2_mut = Mutation(Gamete_1.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_2.1_mut = Mutation(Gamete_2.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_2.2_mut = Mutation(Gamete_2.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
        
        ploidy = c(rep("Tetraploid", times = 4))
        
        Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
        Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
        
        if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
        else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
      }
      else{
        Gamete_2.1_temp = Diploid_ind[((2*Parent_2)-1),1:(Locus_trait+Locus_UG+Locus_UG)]
        Gamete_2.2_temp = Diploid_ind[(2*Parent_2),1:(Locus_trait+Locus_UG+Locus_UG)]
        
        Gamete_1.1_temp = c(NULL)
        Gamete_1.2_temp = c(NULL)
        
        for(l in 1:(Locus_trait+Locus_UG+Locus_UG)){
          
          temp_ind_1 = c(Tetraploid_ind[4*(Parent_1-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-1,l])
          
          Gamete_1_temp = c(sample(temp_ind_1, 2, replace = F))
          
          Gamete_1.1_temp[l] = Gamete_1_temp[1]
          Gamete_1.2_temp[l] = Gamete_1_temp[2]
        }
        
        Gamete_1.1_mut = Mutation(Gamete_1.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_1.2_mut = Mutation(Gamete_1.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_2.1_mut = Mutation(Gamete_2.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
        Gamete_2.2_mut = Mutation(Gamete_2.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
        
        ploidy = c(rep("Tetraploid", times = 4))
        
        Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
        Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
        
        if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
        else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}}
    }
    else if(Ploidy_parent1 == "tetraploid" && Ploidy_parent2 == "tetraploid" && Gamete_parent1 == "diplo" && Gamete_parent2 == "diplo"){
      
      # Fourth case, two tetraploids parents making reduced gametes (tetraploid offsprings)
      
      Nb_offsprings = Nb_offsprings + 1
      
      Gamete_1.1_temp = c(NULL)
      Gamete_1.2_temp = c(NULL)
      
      Gamete_2.1_temp = c(NULL)
      Gamete_2.2_temp = c(NULL)
      
      for(l in 1:(Locus_trait+Locus_UG+Locus_UG)){
        
        temp_ind_1 = c(Tetraploid_ind[4*(Parent_1-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-1,l])
        temp_ind_2 = c(Tetraploid_ind[4*(Parent_2-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-1,l])
        
        Gamete_1_temp = c(sample(temp_ind_1, 2, replace = F))
        Gamete_2_temp = c(sample(temp_ind_2, 2, replace = F))
        
        Gamete_1.1_temp[l] = Gamete_1_temp[1]
        Gamete_1.2_temp[l] = Gamete_1_temp[2]
        
        Gamete_2.1_temp[l] = Gamete_2_temp[1]
        Gamete_2.2_temp[l] = Gamete_2_temp[2]
      }
      
      Gamete_1.1_mut = Mutation(Gamete_1.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_1.2_mut = Mutation(Gamete_1.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_2.1_mut = Mutation(Gamete_2.1_temp, U, Locus_trait,Locus_UG, var.add.eff)
      Gamete_2.2_mut = Mutation(Gamete_2.2_temp, U, Locus_trait,Locus_UG, var.add.eff)
      
      ploidy = c(rep("Tetraploid", times = 4))
      
      Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
      else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
    }
    
    if(Nb_offsprings == Npop){break}   
  }
  
  All_ind_temp = rbind(diplo_ind_temp, tetra_ind_temp)  
  
  return(All_ind_temp)
  
}

proba_unreduced_gametes_female <-function( Locus_trait,Locus_UG, Npop, Diploid_ind, Tetraploid_ind){
  
  # Prepare the production of UG for female function, for diploids and tetraploids individuals
  
  list_proba_unreduced = c(NULL)
  
  k = 0
  
  if(nrow(Diploid_ind) !=0){
    
    for(i in 1:(nrow(Diploid_ind)/2)){
      
      k = k + 1
      
      val_temp = sum(abs(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+1):(Locus_trait+Locus_UG)]))) + sum(abs(as.numeric(Diploid_ind[(2*i), (Locus_trait+1):(Locus_trait+Locus_UG)])))

      if(val_temp < 0){val_temp = 0}
      else if(val_temp > 1){val_temp = 1}
      else{val_temp = val_temp}
      
      list_proba_unreduced[k] = val_temp
    }
  }
  
  if(k < Npop && is.null(nrow(Tetraploid_ind)) == F){
    for(i in 1:(nrow(Tetraploid_ind)/4)){
      k = k + 1
      
      val_temp = sum(abs(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+1):(Locus_trait+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+1):(Locus_trait+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+1):(Locus_trait+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+1):(Locus_trait+Locus_UG)])))

      if(val_temp < 0){val_temp = 0}
      else if(val_temp > 1){val_temp = 1}
      else{val_temp = val_temp}
      
      list_proba_unreduced[k] = val_temp
    }
  }
  
  return(list_proba_unreduced)
  
}

proba_unreduced_gametes_male <-function( Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type){
  
  # Prepare the production of UG for male function, for diploids and tetraploids individuals
  
  list_proba_unreduced = c(NULL)
  
  k = 0
  
  if(pleio_type != 0){
    if(nrow(Diploid_ind) !=0){
      
      for(i in 1:(nrow(Diploid_ind)/2)){
        
        k = k + 1
        
        #val_temp = sum(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+1):(Locus_trait+pleiotropy)]))+ sum(as.numeric(Diploid_ind[(2*i), (Locus_trait+1):(Locus_trait+pleiotropy)]))+sum(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Diploid_ind[(2*i), (Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)]))
        val_temp = sum(abs(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+1):(Locus_trait+pleiotropy)])))+ sum(abs(as.numeric(Diploid_ind[(2*i), (Locus_trait+1):(Locus_trait+pleiotropy)])))+sum(abs(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Diploid_ind[(2*i), (Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)])))

        if(val_temp < 0){val_temp = 0}
        else if(val_temp > 1){val_temp = 1}
        else{val_temp = val_temp}
        
        list_proba_unreduced[k] = val_temp
      }
    }
    
    if(k < Npop && is.null(nrow(Tetraploid_ind)) == F){
      for(i in 1:(nrow(Tetraploid_ind)/4)){
        k = k + 1
        
        #val_temp = sum(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+pleiotropy+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))
         val_temp = sum(abs(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+1):(Locus_trait+pleiotropy)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+1):(Locus_trait+pleiotropy)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+1):(Locus_trait+pleiotropy)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+1):(Locus_trait+pleiotropy)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+Locus_UG+pleiotropy+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+pleiotropy+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])))

        if(val_temp < 0){val_temp = 0}
        else if(val_temp > 1){val_temp = 1}
        else{val_temp = val_temp}
        
        list_proba_unreduced[k] = val_temp
      }
    }
  }
  else{
  if(nrow(Diploid_ind) !=0){
    
    for(i in 1:(nrow(Diploid_ind)/2)){
      
      k = k + 1
      
      #val_temp = sum(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Diploid_ind[(2*i), (Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))
       val_temp = sum(abs(as.numeric(Diploid_ind[(2*i)-1, (Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Diploid_ind[(2*i), (Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])))

      if(val_temp < 0){val_temp = 0}
      else if(val_temp > 1){val_temp = 1}
      else{val_temp = val_temp}
      
      list_proba_unreduced[k] = val_temp
    }
  }
  
  if(k < Npop && is.null(nrow(Tetraploid_ind)) == F){
    for(i in 1:(nrow(Tetraploid_ind)/4)){
      k = k + 1
      
      #val_temp = sum(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])) + sum(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))
       val_temp = sum(abs(as.numeric(Tetraploid_ind[(4*i)-3,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-2,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i)-1,(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)]))) + sum(abs(as.numeric(Tetraploid_ind[(4*i),(Locus_trait+Locus_UG+1):(Locus_trait+Locus_UG+Locus_UG)])))

      if(val_temp < 0){val_temp = 0}
      else if(val_temp > 1){val_temp = 1}
      else{val_temp = val_temp}
      
      list_proba_unreduced[k] = val_temp
    }
  }
}
  
  return(list_proba_unreduced)
  
}

# A function to compute the phenotypic values of individuals

Phenotype.offspring <-function(Genotype, Npop)
{
  
  pheno.off = c(NULL)
  
  # An empty vector of phenotypic value
  
  for(i in 1:Npop){
    
    pheno.off[i] = Genotype[i] + rnorm(1, 0, 1) # genotype + environmental effects
    
    # Adding an environmental effect to the genotypic value to have the phenotypic value
    
  }
  
  return(pheno.off)
  
}

# A function to compute the genotypic values of individuals

genotype.offspring <-function(Diploid_ind, Tetraploid_ind, dosage, Npop, Locus_trait, pleiotropy, pleio_type)
{
  
    geno.off = c(NULL)
    
    ## An empty vector of genotypic values
    
    k = 0
    
    if(pleio_type == 2){
      if(nrow(Diploid_ind) !=0){
        for(i in 1:(nrow(Diploid_ind)/2)){
          
          ## With additivity, genotype is just the sum of all values stored in the genome file
          
          k = k + 1
          geno.off[k] = sum(as.numeric(Diploid_ind[2*i,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Diploid_ind[(2*i) - 1,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Diploid_ind[2*i,(pleiotropy+1):Locus_trait])) + sum(as.numeric(Diploid_ind[(2*i) - 1,(pleiotropy+1):Locus_trait]))
          
        }
      }
      
      if(k < Npop && is.null(nrow(Tetraploid_ind)) == F){
        for(i in 1:(nrow(Tetraploid_ind)/4)){
          
          k = k + 1
          geno.off[k] = (1+dosage)*(sum(as.numeric(Tetraploid_ind[4*i,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i) - 1,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i) - 2,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[(4*i) - 3,(Locus_trait+1):(Locus_trait+pleiotropy)])) + sum(as.numeric(Tetraploid_ind[4*i,(pleiotropy+1):Locus_trait])) + sum(as.numeric(Tetraploid_ind[(4*i) - 1,(pleiotropy+1):Locus_trait])) + sum(as.numeric(Tetraploid_ind[(4*i) - 2,(pleiotropy+1):Locus_trait]))+ sum(as.numeric(Tetraploid_ind[(4*i) - 3,(pleiotropy+1):Locus_trait])))
          
        }
      }
    }
    else{
    if(nrow(Diploid_ind) !=0){
    for(i in 1:(nrow(Diploid_ind)/2)){
      
      ## With additivity, genotype is just the sum of all values stored in the genome fill
      
      k = k + 1
      geno.off[k] = sum(as.numeric(Diploid_ind[2*i,1:Locus_trait])) + sum(as.numeric(Diploid_ind[(2*i) - 1,1:Locus_trait]))
      
      }
    }
  
    if(k < Npop && is.null(nrow(Tetraploid_ind)) == F){
    for(i in 1:(nrow(Tetraploid_ind)/4)){
      
      k = k + 1
      geno.off[k] = (1+dosage)*(sum(as.numeric(Tetraploid_ind[4*i,1:Locus_trait])) + sum(as.numeric(Tetraploid_ind[(4*i) - 1,1:Locus_trait]))+ sum(as.numeric(Tetraploid_ind[(4*i) - 2,1:Locus_trait])) + sum(as.numeric(Tetraploid_ind[(4*i) - 3,1:Locus_trait])))
      
      }
    }
  }
      
  return(geno.off)
}

# A function to check if population is at equilibrium

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# The main simulation function

simulation_model <-function(Locus_trait,Locus_UG, Npop, Nb_rep, U, pleiotropy, pleio_type, env_effect, var.add.eff, om_2, dosage){
  
  rep = 0 # initiating the counting of the repetitions
  
  for (r in 1:Nb_rep){
    
    rep = rep + 1 # Add +1 to the repetition counter
    
    ### Initialisation of the genome for the fitness trait
    
    # For simplicity, we merged the phenotypic trait and UG production in a single genome
    
    geno.init = c(rep(0, times=(Locus_trait+Locus_UG+Locus_UG))) ## Generating the null haplotypes with only zero values at all loci
    
    for(i in 1:(2*Npop)){
      
      ## Npop is population size, and merge the null haplotypes to generation N individuals genetically similar at the beginning of simulations
      
      if(i == 1)(genome = geno.init)
      else{genome = rbind(genome, geno.init)}
      
    }
    
    phenotype = c(rnorm(Npop, 0, 1)) ## A first vector of phenotypic values (genotype = 0, and a random environmental effect with rnorm)
    
    ## At t = 0, the population is diploid, and only produce haploid gametes
    
    Ploidy = c(rep("Diploid", times = 2*Npop))
    genome = cbind(genome, Ploidy)
    
    Diploid_ind = subset(genome, genome[,Locus_trait+Locus_UG+Locus_UG+1]=="Diploid")
    Tetraploid_ind = c(NULL)
    
    # Functions
    
    gen = 0 ## Initiating the counting of generations 
    
    # Vectors to store the variable of interest
    
    freq_tetraploids = c(NULL) # Will store frequency of tetraploids at equilibrium and during environmental change
    fitness_population = c(NULL) # Will store mean population fitness at equilibrium and during environmental change
    UG_production = c(NULL) # Will store frequency of Unreduced Gametes at equilibrium and during environmental change
    var_UG_production = c(NULL) # Will store the genetic variance for Unreduced Gametes production at equilibrium and during environmental change
    
    step = 1 # Step 1 is the population under stabilizing selection, before the environmental change
    
    mean.fitness = c(NULL) # Will store the mean fitness at each generation, use to check equilibrium
    
    repeat{
      
      gen = gen + 1 ## Adding one generation to the counter
      
      fitness = exp(-(phenotype^2)/(2*om_2)) ## Computing the fitness of individuals based on phenotype-fitness function
      
      mean.fitness[gen] = mean(fitness) ## Inferring the mean fitness value of the population at generation t
      
      ## To understand the functions, have a look to the annotations in each function 
      
      ## To prepare the new generation we are making 1) Reproduction
      ## 2) Inferring the genotypes of offsprings
      ## 3) Inferring the phenotypes of offsprings
      
      # 1) Reproduction
      
      genome = reproduction(Diploid_ind, Tetraploid_ind, Locus_trait,Locus_UG, U, Npop, fitness, step, env_effect, var.add.eff, pleiotropy, pleio_type)
      
      # Seperation of diploids and tetraploids individuals
      
      Diploid_ind = subset(genome, genome[,Locus_trait+Locus_UG+Locus_UG+1]=="Diploid")
      Tetraploid_ind = subset(genome, genome[,Locus_trait+Locus_UG+Locus_UG+1]=="Tetraploid")
      
      # 2) Genotypes 
      
      genotype = genotype.offspring(Diploid_ind, Tetraploid_ind, dosage, Npop, Locus_trait, pleiotropy, pleio_type)
      
      # 3) Phenotypes
      
      phenotype = Phenotype.offspring(genotype, Npop)
      
      ## Check if you have run enough generation
      ## And check if we reached the equilibirum for genetic diversity
      
      if( is.wholenumber(gen/1000) == TRUE && (gen/1000) > 1 ){
        
        mean.1<-mean(mean.fitness[(gen-999):gen]) ## Variance of the first 1000 generations
        mean.2<-mean(mean.fitness[(gen-1999):(gen-1000)]) ## Variance of the second 1000 generations
        
        if( abs(1 - (mean.1/mean.2)) <= 0.01){break} ## Have a look if it's deviate to much or not
        
      }
      
    }
    
    ## Store quantities at equilibrium
    
    freq_tetraploids[1] = (nrow(Tetraploid_ind)/4)/Npop
    fitness_population[1] = mean(fitness)
    UG_production[1] = (mean(proba_unreduced_gametes_male(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type))+ mean(proba_unreduced_gametes_female(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind)))/2
    var_UG_production[1] = (var(proba_unreduced_gametes_male(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type))+ var(proba_unreduced_gametes_female(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind)))/2 
      
    ## Second step = introduction of the environmental change (on optimum and UG production)
    ## Possibility to remove mutations or not (U = 0 or not)
    
    step = 2
    opt = 5
    
    # Follow the population for g generations
    
    for(g in 1:500){
    
    fitness = exp(-((phenotype-opt)^2)/(2*om_2)) ## Computing the fitness of individuals based on phenotype-fitness function
    mean.fitness[g] = mean(fitness) ## Inferring the mean fitness value of the population at generation t
    
    ## To understand the functions, have a look to the annotations in each function 
    
    ## To prepare the new generation we are making 1) Reproduction
    ## 2) Inferring the genotypes of offsprings
    ## 3) Inferring the phenotypes of offsprings
    
    # 1) Reproduction
    
    genome = reproduction(Diploid_ind, Tetraploid_ind, Locus_trait,Locus_UG, U, Npop, fitness, step, env_effect, var.add.eff, pleiotropy, pleio_type)
    
    # Seperation of diploids and tetraploids individuals
    
    Diploid_ind = subset(genome, genome[,Locus_trait+Locus_UG+Locus_UG+1]=="Diploid")
    Tetraploid_ind = subset(genome, genome[,Locus_trait+Locus_UG+Locus_UG+1]=="Tetraploid")
    
    # 2) Genotypes 
    
    genotype = genotype.offspring(Diploid_ind, Tetraploid_ind, dosage, Npop, Locus_trait, pleiotropy, pleio_type)
    
    # 3) Phenotypes
    
    phenotype = Phenotype.offspring(genotype, Npop)
    
    # Store quantities at each generations during adaptation process
    
    freq_tetraploids[g+1] = (nrow(Tetraploid_ind)/4)/Npop
    fitness_population[g+1] = mean(fitness)
    UG_production[g+1] = (mean(proba_unreduced_gametes_male(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type))+ mean(proba_unreduced_gametes_female(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind)))/2
    var_UG_production[g+1] = (var(proba_unreduced_gametes_male(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind, pleiotropy, pleio_type))+ var(proba_unreduced_gametes_female(Locus_trait,Locus_UG, Npop,Diploid_ind, Tetraploid_ind)))/2 
    
    }

## Choose what you want to store    
 
    if(rep == 1){
      
      final_freq_tetra = freq_tetraploids
      final_fitness = fitness_population
      final_mean_UG = UG_production
      final_var_UG = var_UG_production
    }
    
    else{
      
      final_freq_tetra = rbind(final_freq_tetra, freq_tetraploids)
      final_fitness = rbind(final_fitness, fitness_population)
      final_mean_UG = rbind(final_mean_UG, UG_production)
      final_var_UG = rbind(final_var_UG, var_UG_production)
    }
       
}
  
  ### Write the outputs of the model as .txt files
  
  write.table(final_freq_tetra, file=paste("freq_tetra_L",Locus_trait,"_LUG",Locus_UG,"_N",Npop,"_U",U,"_pleio",pleiotropy,pleio_type,"_env",env_effect,".txt",sep=""),row.names = F, dec = ".")
  write.table(final_fitness, file=paste("fitness_L",Locus_trait,"_LUG",Locus_UG,"_N",Npop,"_U",U,"_pleio",pleiotropy,pleio_type,"_env",env_effect,"_var",var.add.eff,".txt",sep=""),row.names = F, dec = ".")
  write.table(final_mean_UG, file=paste("UG_L",Locus_trait,"_LUG",Locus_UG,"_N",Npop,"_U",U,"_pleio",pleiotropy,pleio_type,"_env",env_effect,".txt",sep=""),row.names = F, dec = ".")
  write.table(final_var_UG, file=paste("VG_UG_L",Locus_trait,"_LUG",Locus_UG,"_N",Npop,"_U",U,"_pleio",pleiotropy,pleio_type,"_env",env_effect,".txt",sep=""),row.names = F, dec = ".")

  
}

# Example of simulation parameters set

# The parameters are:

# Locus_trait = Number of loci for the quantitative trait
# Locus_UG = Number of loci for Unreduced gametes production (male or females)
# Npop = demographic population size (constant)
# Nb_rep = number of repetitions to perfom for the parameter set
# U = haploid genomic mutation rate
# pleiotropy = number of pleiotrop loci (has to be logical with the number of loci selected before)
# pleio_type = type of pleiotropy, if 0, no pleiotroy, if 1 pleiotropy on male and female UG production, if 2, pleiotropy on all traits
# env_effect = effect of new environment on UG production
# var.add.eff = variance of mutational effects
# om_2 = strength of stabilizing selection (the smaller the stronger the selection)
# dosage = To control for the gigas effect of polyploids

test = simulation_model(Locus_trait = 100,Locus_UG = 100, Npop= 200, Nb_rep= 100, U=0.005, pleiotropy = 0, pleio_type = 0, env_effect = 0.2, var.add.eff = 0.05, om_2 = 9, dosage = 1)



