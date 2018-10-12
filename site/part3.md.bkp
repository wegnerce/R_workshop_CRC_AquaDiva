# Steps beyond

## Going beyond the boxplot - correlation plots

One of the strong suits of _ggplot2_ is that it comprises numerous vizualization options, with dot plots and box plots being only two of them. And beyond that there are meanwhile a lot of R packages that stick to the _ggplot2_ syntax that expand R's plotting capabilities even more. 

Ok, what we will try now is to get in a very simple way an idea whether our parameters are correlated with each other.

Quick reminder, part of our data is made up by descriptive variables. So first, we will extract only our measurement variables.

```
# What are the dimensions of our dataframe
dim(aqd_mock)
# Quickly check again which columns are containing descriptive variables
head(aqd_mock)
# Subset the dataframe accordingly
aqd_num <- aqd_mock[7:24]
```
Done. Now we calculate correlations between all measurement variables.

```
# Some more necessary R packages
library("ellipses")
library("RColorBrewer")

# Calculate correlations
aqd_cor = cor(aqd_num)

# A sneek peek at our correlations
aqd_cor

# We want to colorize our planned correlation plot, so lets create a palette
my_colors <- brewer.pal(5, "Spectral")
my_colors=colorRampPalette(my_colors)(100)

# Plot the plot ;-)
ord <- order(aqd_cor[1, ])
aqd_ord = aqd_cor[ord, ord]
plotcorr(aqd_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1)  )
```

Oh hallo, that's pretty, what does it mean? Let's break up these lines.

```
# Calculate correlations
aqd_cor = cor(aqd_num)
```

!!! QUESTION
	* What does cor() do?

That was an easy one. What about:

```
# We want to colorize our planned correlation plot, so lets create a palette
my_colors <- brewer.pal(5, "Spectral")
my_colors=colorRampPalette(my_colors)(100)
```

!!! QUESTION
	* What is a palette? 
	* What do brewer.pal() and colorRampPalette() do?

```
# Plot the plot ;-)
ord <- order(aqd_cor[1, ])
aqd_ord = aqd_cor[ord, ord]
plotcorr(aqd_ord , col=my_colors[aqd_ord*50+50], type = "lower", diag = FALSE, numbers = TRUE , mar=c(1,1,1,1))
```

!!! EXERCISE
	* Take a moment and try to figure out what the figure shows you

## Going beyond the boxplot while going back

Last but not least we want to take a look at interactive plots. Plots do not have to static, interactive plots allow us to dive into data into a much more engaging way. 

Good, what we now try is to turn or facetted box plot into an interactive version.

Luckily, this is extremly easy.

```
p <- ggplot(aqd_long, aes(x=Well, y=value)) + geom_boxplot() +
  geom_point(aes(fill = Season, alpha = 0.6, shape = factor(Aquifer)), colour = "Black", position = position_jitterdodge()) +
  labs(title="Dataset parameters") + theme(axis.text.x=element_text(angle = 25, hjust = 1)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~Parameter, scales = "free") + ylab(label = "mg/L") + 
  guides(size = FALSE, alpha = FALSE, fill=guide_legend(title="Season"), shape = guide_legend(title="Aquifer")) +
  scale_shape_manual(values=c(21,22))

p <- ggplotly(p)
p
```

You are now able to explore your data interactively in the viewer window of RStudio. For more ideas about interactive _ggplot2_ plots check out this [link](https://plot.ly/ggplot2/).

!!! EXERCISE
	Create interactive boxplots:

	* For a single variable boxplot (e.g. Fe)
	* And a subset

