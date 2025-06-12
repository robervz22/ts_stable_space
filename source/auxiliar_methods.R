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

# visualise estimations
plot_estimates_comparison <- function(Y_obs, Y_est_list, col_names,
                                      labels = NULL, dates = NULL,
                                      method_palette = NULL) {  
  Y_obs <- data.table::as.data.table(Y_obs)
  n <- nrow(Y_obs)
  
  if (is.null(dates)) {
    dates <- 1:n
  }
  
  if (is.null(labels)) {
    labels <- names(Y_est_list)
  }
  
  # Combine estimated values
  plot_data <- rbindlist(lapply(seq_along(Y_est_list), function(i) {
    est <- data.table::as.data.table(Y_est_list[[i]])
    dt_list <- lapply(col_names, function(col) {
      data.table::data.table(
        date = dates,
        value = est[[col]],
        variable = col,
        method = labels[i]
      )
    })
    rbindlist(dt_list)
  }))
  
  # Add observed values
  obs_data <- rbindlist(lapply(col_names, function(col) {
    data.table::data.table(
      date = dates,
      value = Y_obs[[col]],
      variable = col,
      method = "Observed"
    )
  }))
  
  # Combine all
  full_data <- rbind(plot_data, obs_data)
  
  # Default colour palette if not supplied
  if (is.null(method_palette)) {
    method_palette <- c(
      "Observed" = "black",
      "PLS" = "#1b9e77",
      "PCA" = "#d95f02",
      "Johansen" = "#7570b3"
    )
  }

  # Plot
  ggplot(full_data, aes(x = date, y = value, colour = method)) +
    geom_line() +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = method_palette) +
    theme_minimal() +
    labs(title = "",
         x = "Date", y = latex2exp::TeX(r'($X_n$)'), colour = "Method")+mytheme
}

