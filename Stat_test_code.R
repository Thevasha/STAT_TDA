library("dplyr")
library("ggplot2")


#Total number of time steps for each simulation
total_time<-1000
delta_t<-0.02 #Time Increment
total_steps<-total_time/delta_t
snap<-10 #Number of snapshots to be stored
Time_frame_unit=total_steps/(snap)




###################################################
#Code part-3: Statistical inference
####################################################
#Non-parametric functional test to compare PLCs of random vs observed trajectory data
#Load the persistent landscape and dimension of persistent landscape data of the simulated random movement data



load("New_folder/Persistence_landscapes_dim1/Sample_landscape_list.Rdata")
load("New_folder/Persistence_landscapes_dim1/Sample_dimension_list.Rdata")


Sample_dimension_list1<-Sample_dimension_list
Sample_landscape_list1<-Sample_landscape_list

remove(Sample_dimension_list)
remove(Sample_landscape_list)

load("New_folder/Persistence_landscapes_dim1/Sample_landscape_list.Rdata")
load("New_folder/Persistence_landscapes_dim1/Sample_dimension_list.Rdata")


Sample_dimension_list2<-Sample_dimension_list
Sample_landscape_list2<-Sample_landscape_list

Dim1<-min(matrix(unlist(Sample_dimension_list1),ncol=1,byrow=TRUE))
Dim2<-min(matrix(unlist(Sample_dimension_list2),ncol=1,byrow=TRUE))

Dim<-min(Dim1,Dim2)


Global_p_val<-c()
Global_unadj_p_val<-c()
Global_effect_size<-c()
for ( time1 in 1:(snap+1)){
  
  t1=1+((time1-1)*Time_frame_unit)

  list1<-list()
  for (s1 in 1:Sample_size){list1[[s1]]<-rowMeans(Sample_landscape_list1[[s1]][[t1]][,1:Dim])}
  cA<-matrix(unlist(list1), ncol = 1000, byrow = TRUE)
  
  list2<-list()
  for (s2 in 1:Sample_size){list2[[s2]]<-rowMeans(Sample_landscape_list2[[s2]][[t1]][,1:Dim])}
  cB<-matrix(unlist(list2), ncol = 1000, byrow = TRUE)
  
  #First compute observed test stat   
  Combined_samples<-rbind(cA,cB)
  curveA<-Combined_samples[1:Sample_size,]
  curveB<-Combined_samples[(Sample_size+1):(2*Sample_size),]
  
  
  ##Compute pointwise Wilcoxin test statistic 
  Filt_div<-1000#Number of division for filtration
  N<-1000/Filt_div
  
  
  #Step-1:perform the pointwise test at every epsilon (filtration scale) and obtain the unadjusted p-values. The unadjusted p-values ordered from min to max
  
  pointwise_Wilcoxin_test<-rep(0,Filt_div)
  pointwise_test_Wilcoxin_test<-rep(0,Filt_div)
  Effect_size<-rep(0,Filt_div)
  
  for (j in 1:Filt_div){
    
    
    Time_g1<-curveA[,(1+((j-1)*N))]
    Time_g2<-curveB[,(1+((j-1)*N))]
    Wilk_test<-wilcox.test(Time_g1,Time_g2)
    #Compute pointwise_p_values
    pointwise_Wilcoxin_test[j]<- Wilk_test$p.value
    pointwise_test_Wilcoxin_test[j]<- Wilk_test$statistic
    #compute effect size
    Z_stat<-qnorm(Wilk_test$p.value/2)
    
    Effect_size[j]<-abs(Z_stat)/sqrt(Sample_size)
    
    
  }
  
  Test_max<-max(pointwise_test_Wilcoxin_test)
  t_maxT<-which(pointwise_test_Wilcoxin_test==Test_max)
  t_maxT<-t_maxT[1]
  
  Obs_pointwise_Wilcoxin_test<-pointwise_Wilcoxin_test
  #Sorted obs_test_stats
  a<-Obs_pointwise_Wilcoxin_test
  
  a[!is.na(a)]<-sort(a)
  #Sorted Index
  index_order<-Obs_pointwise_Wilcoxin_test
  index_order[!is.na(index_order)]<-order(index_order)
  
  
  
  ##########Compute bootstrapped test statistic, observed p-values for 1000 random tests#########################
  
  #Step2:  Initialize counting variables Ci = 0, i = 1, 2, . . . , k k=100
  Ci<-rep(0,Filt_div)
  
  #Step3: Randomly permute data between the two populations and call the resulting data set as a randomized data set. The p-values computed from a randomized data set are denoted by p*
  
  
  
  for(perm in 1:1000){
    
    #First compute # bootstrapped p-values   
    Random_sample<-Combined_samples[sample(nrow(Combined_samples), 2*Sample_size), ]
    curveA<-Random_sample[1:Sample_size,]  
    curveB<-Random_sample[(Sample_size+1):(2*Sample_size),]
    
    
    
    
    
    boot_pointwise_Wilcoxin_test<-rep(NA,Filt_div)
    boot_pointwise_Wilcoxin_test1<-rep(NA,Filt_div)
    
    boot_Wilk_test_statistic<-rep(0,Filt_div)
    
    
    for (j in 1:Filt_div){
      
      
      bTime_g1<-curveA[,(1+((j-1)*N))]
      bTime_g2<-curveB[,(1+((j-1)*N))]
      bWilk_test<-wilcox.test(bTime_g1,bTime_g2)
      #Compute pointwise_p_values
      boot_pointwise_Wilcoxin_test[j]<- bWilk_test$p.value
      
      
      
      
    }
    
    for (ord in 1:Filt_div){if(is.na(index_order[ord])==FALSE){boot_pointwise_Wilcoxin_test1[ord]<-boot_pointwise_Wilcoxin_test[index_order[ord]]}}
    #Step 4
    boot_pointwise_Wilcoxin_test2<- boot_pointwise_Wilcoxin_test1
    for (k in (Filt_div-1):1){
      if(is.na(boot_pointwise_Wilcoxin_test1[k])==FALSE && is.na(boot_pointwise_Wilcoxin_test1[k+1])==FALSE  ){
        boot_pointwise_Wilcoxin_test2[k]<-min(boot_pointwise_Wilcoxin_test1[k],boot_pointwise_Wilcoxin_test1[k+1])}}
    
    b<- boot_pointwise_Wilcoxin_test2
    
    #Step 5
    for (no in 1:Filt_div){ if (is.na(b[no])==FALSE && is.na(a[no])==FALSE){if (b[no]<=a[no]){Ci[no]<-sum(Ci[no],1)}}}
    Ci<-Ci
  }
  
  #Step 6:  The adjusted p-value is computed as
  pi_N<-Ci/1000
  for (index in 1:Filt_div){if(is.na(index_order[index])==TRUE){pi_N[index]<-NA}}
  
  #STEP 7: Enforce monotonicity using successive maximization:
  
  p_val<-pi_N
  
  for (v in 2:(Filt_div)){
    if(is.na(pi_N[v-1])==FALSE && is.na(pi_N[v])==FALSE  )
    { p_val[v]<-max(pi_N[v-1],pi_N[v])}}
  
  
  
  #Plot adj vs unadj p-values on same graph
  
  x_val=seq(1,Filt_div)
  Adjusted<-rep(0,Filt_div)
  for(a1 in 1:Filt_div){Adjusted[a1]<-p_val[index_order[a1]]}
  Unadjusted<-rep(0,Filt_div)
  for(a2 in 1:Filt_div){Unadjusted[a2]<-a[index_order[a2]]}
  
  
  
  Global_p_val<-c(Global_p_val,Adjusted[t_maxT])
  Global_effect_size<-c(Global_effect_size,Effect_size[t_maxT])
  Global_unadj_p_val<-c(Global_unadj_p_val,Unadjusted[t_maxT])
  
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  df <- structure(list(n =  x_val, Adj_p = Adjusted,Unadj_p=Unadjusted), class = "data.frame", row.names = c(NA, -1000L))
  
  #Inputs
  df1<-df %>% mutate(T1 = ifelse(Adj_p<0.05, n, 0))
  df2<-df1 %>% mutate(T2 = ifelse(Unadj_p<0.05, n, 0))


  save(df2, file=sprintf("Dim_1_Local_p_vals%d.Rdata",time1))
  
}


save(Global_p_val, file=sprintf("Dim_1_Global_p_vals%d.Rdata",time1))
save(  Global_unadj_p_val, file=sprintf("Dim_1_Global_unadj%d.Rdata",time1))
save( Global_effect_size, file=sprintf("Dim_1_Global_effect%d.Rdata",time1))



