Lp_distance<-function(Sample_Average_land_curve,p, total_time,sample_idx){
  
  Lp_distance_matrix<-matrix(0,nrow=(total_time),ncol=(total_time))
for (i in 1:(total_time)){
  for(j in 1:(total_time)){
    Land1<-Sample_Average_land_curve[[sample_idx]][[i]]
    Land2<-Sample_Average_land_curve[[sample_idx]][[j]]
    Lp_distance_matrix[i,j]<-(sum(abs(Land1-Land2)^p)/length(Land1))^(1/p)
  }
  
}
data2 <- Lp_distance_matrix
data2
#data2
}