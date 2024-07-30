
Interaction_force<- function(x,y,epsilon,LA,LR,alpha){
 
  # vector between pairs of cells
  r_x <- outer( x, x, "-" )
  r_y <- outer( y, y, "-" )
  
  # distance between pairs of cells
  r_length = sqrt( r_x^2 + r_y^2 )

  
  # pairs of cells within epsilon are labeled by "1", otherwise zero
  is_eps_neighbor <- r_length < epsilon
  
  #evaluate L-J potential for every pair of cells
  Grad_U = Grad( r_length, LA, LR, alpha )

  # components of the interaction force vector  
  F_x = -Grad_U * ( ( r_x ) / r_length ) * is_eps_neighbor
  F_y = -Grad_U * ( ( r_y ) / r_length ) * is_eps_neighbor
  
  Force=sqrt((F_x^2) +( F_y)^2 )
  # sum all cell interactions
  # na.rm = TRUE skips NaNs created by dividing by |r| = 0 above
  F_x = rowSums(F_x, na.rm=TRUE )
  F_y = rowSums(F_y, na.rm=TRUE )

  Int_mat<-cbind(F_x, F_y)



Int_mat
}




