
# Plot with 'barcode'
# create a list with required components

# Plot with 'barcode'
# create a list with required components
landscape_plot<-function(Lam,tq,d){
  
  land <- list(lambda =Lam , tseq = tq, dimension = d)
  # name the class appropriately
  class(land) <- "landscape"
  
  opar <- par(no.readonly=TRUE)
  PL<-plot(land,top.k=40,colored=TRUE)
  par(opar)
  
  PL
}