# From Excel to R

## Preliminary stuff
```
seq(x,y,z)
# creates a vector with the numbers from x to y in step size z, e.g. seq(0,30,5) is equivalent to c(0,5,10,15,20,25,30)
```

## Import data from excel via the clipboard
```
df <- read.table(file="clipboard",sep="\t",header=TRUE,row.names=1, dec=",")
# file= source of the data, in our case the clipboard
# sep= separator between data points, in our case "\t" for tabulator-separated
# header= does the dataset contain headers (column names)? TRUE or FALSE
# row.names= in which column are names of the rows stored? Here 1 for first column, leave away for no row names
# dec= sets the decimal separator, either "," or "."
```

## Set the margins of the plot
```
par(mfrow=, mar=, oma=)
# warning! par sets global parameters, if you change anything in par, it will stay changed until you set it back to the default
# all parameters are optional, you can change one or multiple ones at a time
# mfrow= create multi plot figures; requires a vector of 2, c(2,3) will arrange plots in 2 rows and 3 columns; default is mfrow=c(1,1)
# mar= margins around plot in lines, requires a vector of 4, the order is bottom, left, top, right; default is c(5.1,4.1,4.1,2.1)
# oma= outer margins for multi plot figures; default is c(0,0,0,0)
# if someting seems messed up, run the line below
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1), oma=c(0,0,0,0))
```

## Make a basic x-y-plot
```
plot(x=, y=, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE, main=)
# x=, y= vectors of x and y values; supply only x for categorial values; x=NULL creates an empty plot
# all other parameters are optional
# xlim=, ylim= vectors of length 2 specifying the axis ranges, e.g. xlim=c(0,10) for an x axis from 0 to 10
# xaxs=, yaxs= either "i" or "r", with "r" the axis extends 4% beyond the values specified in xlim/ylim, with "i" not; default is "r"
# xaxt=, yaxt= with "n", no axis is drawn, so the axis can be specified with the axis() function
# xlab=, ylab= title for the axis, e.g. "Length in meters" or "" for no title
# frame.plot= draw a frame around the plot, TRUE or FALSE
# main= string for a plot title, e.g. "My first plot"
# see ?par for more parameters that can be used in plot() and axis(), e.g. regarding line colors, text size, style, color and font
```

## Draw an axis
```
axis(side=, at=, labels=, tck=, lwd=, lwd.ticks=, las=)
# side= specifies which axis to draw, 1 bottom x-axis, 2 left y-axis, 3 top x-axis, 4 right y-axis
# at= vector that specifies where ticks and tick labels should be, best use as at=seq() (see above)
# labels= vector with strings to be drawn at the at= positions; should be same length as at=; not required for numbers
# tck= optional; specifies tick length and direction, default is tck=-0.02 (short ticks to the outside)
# lwd= optional;specifies the line width of the axis line, default is lwd=1
# lwd.ticks= optional; specifies the line with of the tick lines, default is lwd.ticks=1
# las= optional; orientation of the tick labels, las=1 is horizontal, las=3 is vertical
```

## Draw grid lines or other kind of lines (not between datapoints though)
```
abline(a=, b=, h=, v=, coef=, lwd=, col=, lty=)
# best only use a & b or h or v or coef
# a=, b= intercept and slope for a line, e.g. a=1, b=2 for intercept 1, slope 2
# h= coordinates for horizontal lines, e.g. for an axis with at=seq(0,1,0.2), this can be h=seq(0,1,0.2)
# v= same as h but for vertical lines
# coef= as a and b, but intercept and slope are in a vector, e.g. coef=c(1,2) will produce the same line as above
# lwd= optional; see axis
# col= optional; specify the color of the line
# lty= optional; line type, 1 solid, 2 dash, 3 dot, 4 dot dash etc. etc.
```

## Plot points
```
points(x=, y=, pch=, cex=, col=, xpd=)
# x=, y= see plot
# pch= optional; kind of symbols to draw, run next line of code to see pch=0 to pch=25
plot(rep(0,25),pch=0:25)
# cex= optional; size of the symbol, default is 1
# col= optional; color of the symbol
# xpd= optional; draw symbols outside of plot area? TRUE or FALSE, default is FALSE
```

## Plot lines between points
```
lines(x=, y=, lwd=, col=)
# x=, y= see plot, lines will be drawn between the neighboring points
# lwd= optional; line width; default is 1
# col= optional; color of the line
```

## Add a legend
```
legend(x=, y=, legend=, bty="n", pch=, pt.cex=, lwd=, col=, y.intersp=, xpd=)
# x=, y= position of the legend based on the axis of the plot; can be replaced by keywords like "topright" for the top right corner of the plot
# legend= vector of text strings for each legend entry
# bty= optional; box around the legend, bty="n" for no box
# pch= optional; draw points in the legend, see points()
# pt.cex= optional; size of the points; default is 1
# lwd= optional; specify width for a line to be drawn in the legend; leave out for no line
# col= optional; vector of colors for the points and lines
# y.intersp= optional; line spacing in the legend; default is 1
# xpd= optional; draw symols outside of plot area? TRUE or FALSE, default is FALSE
```

## Add error bars to points (example for standard dev SD of y-error-bars) this is a quick way using a function that draws arrows
```
arrows(x0=, y0=, x1=, y1=, length=, angle=90, code=3, col=)
# if x and y are the vectors you used to plot points(), then:
# x0=, y0= vectors for coordinates FROM which to draw bars, e.g. x0=x, y0=y-SD
# x1=, y1= vectors for coordinates TO which to draw bars, e.g. x1=x, y1=y+SD
# length= length of the horizontal lines at the ends of the error bar, default is 0.25
# angle= angle of the horizontal lines at the ends of the error bar, leave at 90 
# code= draw horizontal lines at the top, bottom or both sides of the error bar, leave at 3
# col= color of the error bars
```
 
## Add other lines of text to your plot
```
mtext(text=, side=, line=, font=)
# text= string of text to print in the plot
# side= at which side of the plot and in which orientation the text should be printed, side=1 bottom, =2 left, =3 top, =4 right
# line= on which line the text should be printed, >= 0 means outside of the plot, < 0 means inside of the plot
# font= style of the text; =1 plain, =2 bold, =3 italic, =4 bold italic
# for superscript and subscript, text=expression() can be used; run the lines of code below to see an example
plot(NULL, xlim=c(0,10), ylim=c(0,5), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE)
mtext(text=expression("Normal text"["Subscript text"]*"normal again"^"superscript"*"normal"), side=3, line=-3)
mtext(text=expression("For numbers, this also works without quotation marks "[420]*123^-2), side=3, line=-5)
mtext(text=expression("This can also do that "["bottom"]^"top"), side=3, line=-7)
```

## Bar chart
```
barplot(height=, space=, names.arg=, beside=, horiz=, col=, border=, axes=, axisnames=)
# height= vector or matrix of values for the height of the bars; in a matrix, the values in each column will be stacked into one bar
# when using a dataframe called df as input,  use as.matrix(df); you might also need to transpose the matrix with t(as.matrix(df))
# space= optional; distance between the bars; default is 0.2
# names.arg= optional; vector of strings for bar names to be drawn below x axis; default is NULL
# beside= optional; TRUE or FALSE, if TRUE, values from columns are drawn next to each other and not stacked; default is FALSE
# horiz= optional; TRUE or FALSE, if TRUE, bars are drawn horizontally; default is FALSE
# col= optional; vector of colors for bars, row by row
# border= optional; TRUE or FALSE, if TRUE, a box is drawn around each bar; default is TRUE
# axes= optional; TRUE or FALSE, if TRUE, the y axis is drawn; default is TRUE
# axisnames= optional; TRUE or FALSE, if TRUE, the names of the bars (see names.arg=) are drawn; default is TRUE
```

<sub>Written by Martin Taubert</sub>
<sub>Oct 2018</sub>
