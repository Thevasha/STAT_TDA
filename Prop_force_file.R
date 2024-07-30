
Prop_dir_file<-function(Prop_init_cell,total_steps,Cell_offset,uniform_repolarization_time_steps,ncell){
  
  
  # Propulsion force direction storage  
  Prop_direction<-list()
  
  ## LOOP ACROSS DIFFERENT CELLS
  for (cell_id in 1:ncell){
    
    # start with fixed direction at time_step 0
    Repolarization <-Prop_init_cell
   
    #initialize propulsion direction (angle) storage for each cell
    Prop_cell<-list()  
   
    ## LOOP ACROSS DIFFERENT TIMESTEPS 
    for (t in 1:total_steps)
    {
      k<-t%% (uniform_repolarization_time_steps-Cell_offset[cell_id])
      
      if (k==0){Repolarization<-2*pi*(runif(1, min =0, max = 1))}
     
      Prop_cell[[t]]<-Repolarization
    }
    
    Prop_direction[[cell_id ]]<- Prop_cell
  }
    

Direction<-matrix(0,nrow=ncell,ncol=total_steps)


for (time_id in 1:total_steps){
  for (cell_val in 1:ncell){
    
    Direction[cell_val,time_id]<- Prop_direction[[cell_val]][[time_id]][1]
   
    
  }
}

Direction

}
  