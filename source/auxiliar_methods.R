######################
# auxiliar functions #
######################

# multivariate time series visualization
vis2series<- function(X,Y,dates=NULL,lab_X="X",lab_Y="Y",title="",ncol = 2,legend.position = "none"){
  if(nrow(X)>nrow(Y)){
    X<-X[1:nrow(Y),]
    time<-nrow(Y)
  }
  else{
    Y<-Y[1:nrow(X),]
    time<-nrow(X)
  }
  # data
  data_X <- data.table::melt(as.data.table(X),variable.name="serie", value.name="X")
  data_Y<- data.table::melt(as.data.table(Y),variable.name = "serie", value.name ="Y")
  
  data_XY <- data.table::data.table(serie=data_X$serie,X=data_X$X,Y=data_Y$Y)
  data.table::setnames(data_XY,old=c("serie","X","Y"),new=c("serie",lab_X,lab_Y))
  #printcolnames((data_XY))
  
  data_XY_tidy<-data.table::melt(data_XY)
  data_XY_tidy[,tiempo:=rep(1:time,2*ncol(X))]
  if(!is.null(dates)){
    data_XY_tidy[,dates:=rep(dates,2*ncol(X))]
    # plot
    ggplot2::ggplot(data_XY_tidy,aes(x=dates,y=value))+ggplot2::geom_line(aes(linetype=variable))+ggplot2::facet_wrap(serie~.,ncol=ncol, scales = "free")+ggplot2::labs(x="Date",y=latex2exp::TeX(r'($X_n$)'),color=" ")+ggplot2::scale_color_brewer(palette = "Set1")+
    ggplot2::ggtitle(title)+mytheme+ggplot2::theme(legend.position = legend.position )
  }
  else{
    ggplot2::ggplot(data_XY_tidy,aes(x=tiempo,y=value))+ggplot2::geom_line(aes(linetype=variable))+ggplot2::facet_wrap(serie~.,ncol=ncol, scales = "free")+ggplot2::labs(x="Date",y=latex2exp::TeX(r'($X_n$)'),color=" ")+ggplot2::scale_color_brewer(palette = "Set1")+
    ggplot2::ggtitle(title)+mytheme+ggplot2::theme(legend.position = legend.position)
  }
  
  
}

# mse vectorial function
VNMSE_ <- function(X,X_est){
  f <- function(x){
    return(mean(x^2))
  }
  deviations <-apply(X-X_est,2,f)
  variances <- apply(X,2,var)
  return(deviations/variances)
}
