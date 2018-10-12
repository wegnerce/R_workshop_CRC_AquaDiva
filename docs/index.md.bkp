# Data visualization and exploration

## Synopsis

In this part of the workshop you will familiarize yourself with:

* R's basic plotting capabilities
* _ggplot2_, its syntax and fundamentals
* as well as some advanced data visualization

If you have questions ASK, feel free to drop me an [e-mail](mailto:carl-eric.wegner@uni-jena.de) also after the workshop.

## Setting the stage

Plotting is a very personal thing (btw color schemes as well), ask three different people and you will get a variety of feedback with respect to plotting and data visualization in general. A lot of people get away with _Excel_, _SigmaPlot_, _GraphPad_, _Matplotlib_ (if you are a Pythonista, which makes you a good person by default) or something else. 

R as its parent S is a programming language primarily dedicated to data manipulation (in a good sense) and statistical analysis. Due to its modular architecture and the availabilitry of sophisticated libraries for data processing and plotting it is an excellent choice for any kind of data visualization.

When I got in touch with R the first time, I was googling options how to do a complex multi-panel figure. The original query was "plot facets vector graphics". If you do the exact same query now, you will end up with _ggplot2_ being the second hit. Yep, this is how I "met" R and years later I would still lie if I would say that I'm any sort of R expert. Most of me using R is still centered around the following dogma:

1. Define the problem
2. Look up necessary R resources
3. Apply available examples to own data
4. Be happy after A LOT of try and error.

Ok, first things first, fire up RStudio and check whether you have the following packages installed and load them. Throughout the course you can either execute commands from within a script (mark the respective lines of code and hit CTRL+ENTER) or directly in the console.

```
library("ggplot2")
library("dplyer")
library("tidyr")
library("ellipse")
library("RColorBrewer")
```

Ideally you should see no error message popping up, if you see any then you did not prepare properly for the workshop - shame on you, or rather me, because I did not tell you in time.

![alt text](fail.jpg)

## Basic plotting

Enough blabla, this session is about data visualization, where are the plots? Here they come. Create a new R script, type the following two lines and execute them:

```
dotchart(rnorm(250), col = "blue", main = "Quick ugly example")
hist(rnorm(250), col = "blue", main = "Quick ugly example II")
```
Alternatively type them one after another in the console.

Keep an eye on the plot window, what did just happen?

!!! QUESTION

     * Try to figure out what the individual functions and parameters do.

First take home message of the day - the R help is your friend (a good one).
The general syntax for calling R's help is:

`?functionXYZ()`

Let's start exploring R's basic plotting capabilities in a bit more detail, this is a bit of a recap of what you did yesterday. 
``` 
# Let's define two arbitrary vectors
bacteria <- c(10, 30, 60, 5, 90)
archaea <- c(25, 27, 22, 37, 10)

# Plot them both
plot(bacteria, type="o", col = "orange")
lines(archaea, type="o", pch=22, lty=2, col = "blue")

axis(1, at=1:5, lab=(c("March","April","May","June","July")))
``` 
Let's say what we just plotted are relative abundances for archaea and bacteria across different months. But the output looks like garbage, apparently we had at first default x-axis labels, which were overwritten. The result is this "beauty" of a plot.

How to fix that, take a look at the following code:
``` 
plot(bacteria, type="o", col = "orange", axes=FALSE, ann = FALSE)
lines(archaea, type="o", pch=22, lty=2, col = "blue")

axis(1, at=1:5, lab=(c("March","April","May","June","July")))
``` 
Does that make it any better? What is now missing?

!!! EXERCISE

     * Call the help for axis(), box(), title(), and legend. 
     * Add a y-axis, a box around the plot, titles for the plot as well as the two axes, and a legend. 

One possible solution:

```
plot(bacteria, type="o", col = "orange", axes=FALSE, xlab = "Month", ylab = "Rel. abundance [%]", main = "Bac and Arc")
lines(archaea, type="o", pch=22, lty=3, col = "blue")
axis(1, at=1:5, lab=(c("March","April","May","June","July")))
box()
axis(2, las = 2)
legend(1, max(bacteria), c("Bacteria", "Archaea"), cex=0.8, col=c("orange", "blue"), 
       pch=c(21,22), lty=c(1,3))
```

Alright, so far we played with a dataset that we quickly created, as you already learned before you can easily import datasets like this from any delimited file. Imagine a file like this (e.g. table.tsv):

| bacteria | archaea |
|----------|---------|
| 10       | 25      |
| 30       | 27      |
| 60       | 22      |
| 5        | 37      |
| 90       | 10      |

Two columns, tab-delimited, and the columns have names ("bacteria2, "archaea").

We could read this file as outlined below. 

``` 
# Read the table, pay attention to the header and sep arguments
rel_prok <- read.table("table.tsv", header=T, sep="\t")

# Instead we merge our two vectors because you did a lot of importing yesterday...
rel_prok <- data.frame(bacteria, archaea)
colnames(rel_prok) <- c(bacteria, archaea)

# We define colors to be used with our data series, because why not
plot_colors <- c("blue", "orange")

# We initiate a PNG devide to save the output
png(filename="output.png", height=250, width=300, bg="white")

# AND NOW?!
# Adapt your code
# and end it with
dev.off()
# to turn off the PNG device
``` 

!!! EXERCISE

     * Plot the data as before, and save the output as a .png.
     * Do you have to adjust the dimensions?
     * What is the dev.off() function doing?
     
## Base R data visualization options

So far we basically only did line charts, base R provides us however with a whole range of different visualization options. 

Let's take one of our vectors and see how we can create bar charts and how we can visualize both data series by dot charts.

Bar charts:

```
bacteria <- c(10, 30, 60, 5, 90)

# A simple bar plot
barplot(bacteria, main="Bacteria relative abundance", xlab="Month",
   ylab="Rel. abundance [%]", names.arg=c("March","April","May","June","July"))

# Simple, and slightly pimped, pattern fill
barplot(bacteria, main="Bacteria relative abundance", xlab="Month",  
   ylab="Rel. abundance [%]", names.arg=c("March","April","May","June","July"), 
   border="gray", density=c(10,20,30,40,50))

# Add a box around the plot because we like boxes
box()

# This time with colors
barplot(bacteria, main="Bacteria relative abundance", xlab="Month",  
        ylab="Rel. abundance [%]", names.arg=c("March","April","May","June","July"),
        col=rainbow(5))
```

And a dot chart:

```
# Plot the dotchart
dotchart(t(rel_prok), color=c("blue", "red"), main="Dotchart Bacteria and Archaea")
```

We finish this first session with a little exercise.

!!! EXERCISE

     Use the simple dataset to plot a grouped bar chart incl. a legend and a dot chart with months as row names. Export both as .png files. _NOT_ using RStudio's export function.

For more examples of using the plotting capabilities of base R have a look for instance [here](https://www.harding.edu/fmccown/r/).

Solutions:

```
# Grouped bar chart
barplot(as.matrix(rel_prok), main="Bac vs Arc", ylab= "Rel. abundance",
        beside=TRUE, col=rainbow(5))
box()

# Place the legend at the top-left corner with no frame  
# using rainbow colors
legend("topleft", c("March","April","May","June","July"), cex=1, 
       bty="n", fill=rainbow(5))
```
```
# Dot chart with months as labels
row.names(rel_prok) <- c("March","April","May","June","July")
rel_prok
dotchart(t(rel_prok), color=c("blue", "red"), main="Dotchart Bacteria and Archaea", cex = 1)
```
