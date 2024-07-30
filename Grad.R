Grad<-function (r,LA,LR,alpha){
  CA<-alpha
  CR<-0.25*alpha
  gradU<-((CA/LA)*exp((-r)/LA))-((CR/LR)*exp((-r)/LR))

}

