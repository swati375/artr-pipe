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

plot_piechart<-function(data_features,name)
{
# Make a data frame with region counts
region_counts <- data_features
colnames(region_counts) <- c("Region", "Count")

# Calculate percentage labels
region_counts$Perc <- round(region_counts$Count / sum(region_counts$Count) * 100, 1)
region_counts$Label <- paste0(region_counts$Region, " (", region_counts$Perc, "%)")

# Pie chart with ggplot2
p<-ggplot(region_counts, aes(x = "", y = Count, fill = Region)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = paste0("Peak annotation distribution for ",name)) +
  geom_text(aes(label = Perc), position = position_stack(vjust = 0.5), size = 4)

return (p)
}