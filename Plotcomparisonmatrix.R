
plotcomparisonmatrix <- function(M, times, pcutoff) {
  
  plevels<-c(0,pcutoff,1.0)
  plabels<-c("DIFF",toString(pcutoff),"SAME")
  timelabel <- function(x) paste("t =", toString(x))
  
  datanames=matrix("0",nrow=1,ncol=ncol(M))
  for ( i in 1:length(times) ){datanames[i]<-timelabel(times[i])}
  
  
  library(pheatmap)
  pheatmap(M, display_numbers=TRUE,number_format = "%.3f", cluster_rows = FALSE, cluster_cols = FALSE, drop_levels = TRUE,
           fontsize_number = 8,breaks=plevels, labels_row = datanames, labels_col = datanames,
           legend_breaks = plevels , legend=TRUE, legend_labels=plabels,
           na_col="white", color=c("#990000","#000099","#000099"),number_color="white",angle_col=0,fontsize=8)
  
}
