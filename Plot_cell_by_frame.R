plot_cell_evolution<- function(cell_trajectory, filename,Time_frame_unit){
  
  tot_steps = length(cell_trajectory)
  
  for (frame_id in 1:tot_steps){
    if (frame_id %% Time_frame_unit == 1) {  
    x<-as.vector(cell_trajectory[[frame_id]][,1])
    y<-as.vector(cell_trajectory[[frame_id]][,2])
    # compute distances
    D = sqrt(outer(x,x,"-")^2 + outer(y,y,"-")^2)
    D[D>1]<-0 #Only choose neighbourhood with less than 1 as connected graph
    
    # create a "tidygraph" data structure.
    G <- as_tbl_graph(D) %>% mutate(x=x, y=y) # mutate can be skipped, but useful if one wants to "attach" coordinates to the graph structure
    
    k=list()
    for (m in  1:(length(V(G)))){k[m]=length(which(get.edgelist(G) == V(G)[m]))}
    k=unlist(k)
    
    G1<-as.undirected(G)
    V(G1)$color <- ifelse(k == 0, "red", "blue")
    E(G1)$width <- 3
    V(G1)$size<- 6
    E(G1)$color <- "blue"
    
 
    if (str_length(filename) > 0) {
      png(file=sprintf("%s_%d.png",filename,frame_id))
      plot(G1,vertex.size=V(G1)$size, edge.width=E(G1)$width, edge.color= E(G1)$color, vertex.color=V(G1)$color,vertex.label=NA)
      title(sprintf("Time_frame: %d",frame_id),cex.main=0.95,col.main="black")
      dev.off()
      
    }
    
  }
  
  }
}

