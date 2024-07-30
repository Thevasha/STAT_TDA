
#Required libraries


library(ggplotify)
library(readxl)
library("TDA")
library(ggplot2)
library(reshape2)
library(matrixStats)
library("ggpubr")
library("fda.usc")
library(plotrix)
library(tidyverse)
library(plotly)
library(plyr) 
library("numbers")
library(distributions3)
library("transport")
library('numbers')
library('twiddler')
theme_set(theme_pubr())
library('tibble')
library(gapminder)
library("TDAkit")
library(gganimate)
library("ggpubr")
library("DescTools")
library(igraph)
library(ggraph)
library(tidygraph)
library(sfsmisc)



#Inputs for the simulation of the model

Cell_density=0.5
count=200


ncell= count
#Square domain with periodic boundary conditions.
width_x <- sqrt(ncell/Cell_density)
width_y <- sqrt(ncell/Cell_density)

#Number of cells in the domain

#Number of simulations
Sample_size<-100

M<-1 #Multiple of total_time or number of times the parameter changes in the simulation

#Total number of time steps for each simulation
total_time=1000
delta_t=0.02 #Time Increment
total_steps<-total_time/delta_t
snap<-16 #Number of snapshots to be stored
Time_frame_unit=total_steps/(snap)

#Propulsion force components 
# Magnitude between 0.009-0.025
Prop_values<-c(0.009)  #First parameter
# Direction of the propulsion force varies with offset values
Initial_Prop_dir<-2*pi*(runif(1, min =0, max = 1) ) #Initial direction of the propulsion force
uniform_repolarization_time_steps=2500
Cell_offset<- sample(1:500, ncell, replace=T)


#Interaction force components 
#Magnitude depends on adhesion, epsilon neighbourhood and other fixed constants
Adhesion_values=c(0.09)#second parameter alpha to change attraction-repulsion strength
epsilon=1.5      
LA=14
LR=0.5
eta=1

d=1 #dimension of persistent homology



#Required function file
source("Prop_force_file.R")
source("Grad.R")
source("Interaction_force.R")
source("Plot_cell_by_frame.R")
source("Compute_homology.R")
source("landscape_plot.R")
source("plotcomparisonmatrix.R")
source("Lp.R")



############### Complete code to generate simulations for varied parameter choice and perform statistical test###############

#for (num in 1:10)
#{

dir.create(file.path("New_folder"))

start<-Sys.time()  


alpha_idx=1
prop_idx=1



#create folders for each parameter choice
dir.create(file.path("New_folder/Simulations"))


#############################################################################################################################
#Code part 1: Generate simulations 

##Implement samples of simulations for specific parameters using the self-propelled particle model

start1<-Sys.time()  
Sample_list<-list() #This list contains all cell positions data at each time for each sample

## LOOP ACROSS DIFFERENT SIMULATIONS
for (sample_id in 1:Sample_size){
  
  Cell_positions_input<-matrix(0,nrow=ncell,ncol=2)
  for (i in 1:ncell){Cell_positions_input[i,1]<-((runif(1,min=0,max=width_x)))}
  for (j in 1:ncell){Cell_positions_input[j,2]<-((runif(1,min=0,max=width_y)))}
  
  # initialize trajectory storage  
  
  Cells_trajectory<-list()
  
  # generate the realization of the prop. force vectors
  Propulsion_dir<-Prop_dir_file(Initial_Prop_dir,(M*total_steps),Cell_offset,uniform_repolarization_time_steps,ncell)
  
  # store the initial condition as the "current" step
  Cell_positions_at_t<-Cell_positions_input
  # store the initial condition as the first position in our trajectory
  Cells_trajectory[[1]]  <- Cell_positions_at_t
  
  ## MAIN SIMULATION LOOP
  for (t in 1:(M*total_steps)) 
  {
    if (t %% total_steps == 1) {prop_idx=(((t-1)/total_steps)+1) }
    Prop_magnitude=Prop_values[prop_idx]
    if (t %% total_steps == 1) {alpha_idx=(((t-1)/total_steps)+1) }
    alpha_magnitude=Adhesion_values[alpha_idx]
    # compute the position update forces
    Interact_force_at_t_matrix<-Interaction_force(Cell_positions_at_t[,1], Cell_positions_at_t[,2], epsilon,LA,LR,alpha_magnitude)
    Propulsion_force_at_t_matrix<-cbind((Prop_magnitude*(cos(Propulsion_dir[,t]))),(Prop_magnitude*(sin(Propulsion_dir[,t]))))
    
    # update current cell positions
    Cell_positions_at_t <- Cell_positions_at_t + 
      (delta_t/eta)*(Propulsion_force_at_t_matrix+Interact_force_at_t_matrix)
    
    # enforce periodic boundary conditions
    Cell_positions_at_t[,1] <- Cell_positions_at_t[,1] %% width_x
    Cell_positions_at_t[,2] <- Cell_positions_at_t[,2] %% width_y
    
    Cells_trajectory[[t+1]]<-Cell_positions_at_t
    
    
  }
  
  # store the trajectory into simulation database
  Sample_list [[sample_id]]<-Cells_trajectory
  
  # Save plots of  the cell images at each time step in simulation (sample) file directory
  dir.create(file.path(sprintf("New_folder/Simulations/Sample_%d",num,sample_id)))
  plot_cell_evolution(Cells_trajectory,sprintf("New_folder/Simulations/Sample_%d/Plot",sample_id),Time_frame_unit) 
  
}


end1<-Sys.time()
Time_elapsed1<-end1- start1

save(Sample_list, file="New_folder/Simulations/simulation database.Rdata")
save(Time_elapsed1, file="New_folder/Simulations/Time_elapsed_simulations.Rdata")







#############################################################################################################################
#Code part-2: Compute homology
### Obtain persistence landscapes

d=0
start2<-Sys.time()
dir.create(file.path("New_folder/Persistence_landscapes_dim0"))
Sample_landscape_list<-list() #This list contains all persistence landscape data at each time for each sample
Sample_dimension_list<-list()

for (sample_id in 1:Sample_size){
  
  dir.create(file.path(sprintf("New_folder/Persistence_landscapes_dim0/Sample_%d",sample_id)))
  Landscape_list<-list()
  Dimension_of_landscape<-list()
  
  
  for (frame_id in 1:((M*total_steps)+1)){
    
    if (frame_id %% Time_frame_unit == 1) {
      Data<-Sample_list[[sample_id]][[frame_id]] 
      homology  = diagRips( Data)
      landscape = newdiag2landscape(homology,dimension = d)
      tsequence<-landscape[["tseq"]]
      k_landscapes<-as.matrix(landscape[["lambda"]])
      
      Landscape_list[[frame_id]]<-as.matrix(k_landscapes)
      Dimension_of_landscape[[frame_id]]<-ncol(as.matrix(k_landscapes))
    }
    
    if (frame_id %% Time_frame_unit == 1){
      landscape_plot(k_landscapes,tsequence,d)
      ggsave(sprintf("New_folder/Persistence_landscapes_dim0/Sample_%d/time_frame_%d.png",sample_id,frame_id))}
    
  }
  
  Sample_landscape_list[[sample_id]]<-Landscape_list
  Sample_dimension_list[[sample_id]]<-Dimension_of_landscape
}


end2<-Sys.time()
Time_elapsed2<-end2- start2
save(Sample_landscape_list,file ="New_folder/Persistence_landscapes_dim0/Sample_landscape_list.Rdata") 
save(Sample_dimension_list,file ="New_folder/Persistence_landscapes_dim0/Sample_dimension_list.Rdata")



#Code part-2: Compute homology
### Obtain persistence landscapes
d=1
start2<-Sys.time()
dir.create(file.path("New_folder/Persistence_landscapes_dim1"))
Sample_landscape_list<-list() #This list contains all persistence landscape data at each time for each sample
Sample_dimension_list<-list()

for (sample_id in 1:Sample_size){
  
  dir.create(file.path(sprintf("New_folder/Persistence_landscapes_dim1/Sample_%d",sample_id)))
  Landscape_list<-list()
  Dimension_of_landscape<-list()
  
  
  for (frame_id in 1:((M*total_steps)+1)){
    
    if (frame_id %% Time_frame_unit == 1) {
      Data<-Sample_list[[sample_id]][[frame_id]] 
      homology  = diagRips( Data)
      landscape = newdiag2landscape(homology,dimension = d)
      tsequence<-landscape[["tseq"]]
      k_landscapes<-as.matrix(landscape[["lambda"]])
      
      Landscape_list[[frame_id]]<-as.matrix(k_landscapes)
      Dimension_of_landscape[[frame_id]]<-ncol(as.matrix(k_landscapes))
    }
    
    if (frame_id %% Time_frame_unit == 1){
      landscape_plot(k_landscapes,tsequence,d)
      ggsave(sprintf("New_folder/Persistence_landscapes_dim1/Sample_%d/time_frame_%d.png",sample_id,frame_id))}
    
  }
  
  Sample_landscape_list[[sample_id]]<-Landscape_list
  Sample_dimension_list[[sample_id]]<-Dimension_of_landscape
}


end2<-Sys.time()
Time_elapsed2<-end2- start2
save(Sample_landscape_list,file ="New_folder/Persistence_landscapes_dim1/Sample_landscape_list.Rdata") 
save(Sample_dimension_list,file ="New_folder/Persistence_landscapes_dim1/Sample_dimension_list.Rdata")


