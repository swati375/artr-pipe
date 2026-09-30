library(ggplot2)
library(ggrepel)

plot_deg<- function(df_feature,y,xlabel,ylabel,main,df_genelist)
{
  print(xlabel)
  print(ylabel)
  ## Top 100 peaks by absolute log2FC
  top <- order(
    abs(df_feature),
    decreasing = TRUE
  )[1:100]
  
  plot(
    df_feature,
    y,
    pch = 16,
    xlab = xlabel,
    ylab = ylabel,
    main = main
  )
  
  ## Highlight top 100 peaks
  points(
    df_feature[top],
    top,
    col = "red",
    pch = 16
  )
  
  ## Add gene labels
  text(
    df_feature[top],
    top,
    labels = df_genelist[top],
    col = "red",
    pos = 4,
    cex = 0.6
  )
  
  abline(v = 0, lty = 2)
}



