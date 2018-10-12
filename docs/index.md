# Data Manipulation and Selection with R
Today we will practice slicing and dicing a dataframe to grab only the data we are interested in. I am not going to exhaustively cover all the features in R for searching and manipulating data and instead this will give you a snap shot in the most common methods I use.

For a more complete discussion / further investigation I highly recommend you check out [Hadley Wickham's](http://hadley.nz/) post that can be found <http://adv-r.had.co.nz/Subsetting.html>.

Today we'll be covering:

1. Slicing a dataframe using row and column numbers
2. Grabbing specific named columns or rows
4. Identifying data based on operators (conditions)
5. Using grep to search for specific keywords / regex
6. Sorting data
7. Doing all of the above with the amazing 'dplyr' package.

## Slicing a dataframe
### Specifying row and column numbers
R has 5 basic data structures and we'll cover some of them today. I typically default to using data frames, and these are most similar to a two-dimensional table in excel. A matrix is a similar structure, but all the columns need to be the same type.

Here I've named our data frame (df), however, you could pick any variable name you want, although try not to use special characters or already named functions.

Data from a data frame can be accessed with the syntax
dataframe[rows, columns]

Hopefully the following examples will make this more clear.

First import the data that we'll be using today.
```python
df <- read.delim("simulated_dataset.txt", header=T)
```

Try these examples:
```python
#Get the first 3 rows from the dataframe
df[c(1,2,3),]

#Get the first 3 columns from the dataframe
df[,c(1,2,3)]

#Get the first 3 rows and first 3 columns
df[c(1,2,3),c(1,2,3)]

#Similar to the bash command, head is a function that by default prints the first 6 rows
head(df)
head(df, 10)

#Similarly, tail prints the last 6 rows
tail(df)

#You can chain the head function with specifying only some columns
head(df[,c(1,2,3)])

#Use the minus sign to **remove** columns or rows
df[-c(1,2,3),]
df[,-c(1,2,3)]
```

Now it is your turn, try to answer the following questions:

&#49;.   Grab the 1st, 3rd, and 5th rows of our data frame.  
&#50;. Get all data except the 3rd column.  
&#51;. Retrieve the last 6 rows.  
&#52;. Get the "SampleID", "Well", "Month", "PH", "DO", and "DOC" from all samples from the "HTU" aquifer and store the data in a new dataframe named "df_HTU".  

As a quick aside, we can use the rbind() function to combine two dataframes by rows.
Try the following:
```python
df1 <- df[c(1,2,3),]
df2 <- df[c(28, 29, 30),]

#Use rbind to merge these 2 sliced dataframes
df3 <- rbind(df1, df2)

#Or you can get the same result with these 2 commands:
df3 <- df[c(1,2,3,28,29,30),]
df3 <- df[c(seq(1,3), seq(28,30)),]
```
### Using ranges and formulas to slice a data frame
Instead of specifying each row or column number, we can also generate these based on a range of numbers or by using functions in R that generate a vector of numbers

For example:
```python
#Grab columns 3 through 5
df[,3:5]
#Same outcome, but we can specify directly that the range is a vector
df[,c(3:5)]
#We can mix specifying numbers and generating a range
df[c(1,2,10:15,19),]
#Try
df[1,2,10:15,19,]
```

Once you get more familiar with some R functions that generate vectors we can also chain these together to specify rows or columns

```python
df[,seq(1,ncol(df), by=3)]
```
&#53;. What is the seq command?  
&#54;. What is the ncol command?  


```python
df[,seq(1,length(df), by=3)]
```
&#55;. What is the difference between length and ncol here?  
&#56;. How would you select a range of rows that represent the last 5 rows in the dataframe?  
Answer:
```python
df[c((length(df)-5):length(df)),]
```
And again, we can chain these together into a vector:
```python
df[c(1,2,seq(3,11, by=1)),]
```
Try to fix these commands:  
&#49;&#48;. df[1,2,seq(3,11,2),]  
&#49;&#49;. Why do you not need to specify "by"?  
&#49;&#50;. df[c(1:5,seq(7,10)]  

###Splicing a data frame by using column and/or row names
Here we're using R conventions/shortcuts to pull out specific vectors based on their name.
To start off with lets look at functions that view the names
```python
#Look at the column names
colnames(df)

#Look at the row names
rownames(df)

#What does this do?
names(df)
```

Lets try some examples:
```python
# Grab columns "SampleID" and "PH"
df[,c("SampleID", "PH")]

# What is the difference in the following 2 commands?
df[c(1,2,3),]
df[c("1", "2", "3"),]
```
Example Questions:  
&#49;&#51;. Get a dataframe with the SampleID and the Zone only.  
&#49;&#52;. Retrieve a dataframe with the sampleID and all continuous data columns.  
&#49;&#53;. Get a dataframe with Samples 1 through 6 with SampleID, Aquifer, TOC and TIC.  

## Using Logical and Conditional Operators
We can also select only rows or columns that match some condition we want. This can be something like all rows where the PH is less than 5, or all wells that match a specific value / string.

But to start off with, lets go over some terms and syntax.  
We can refer to a column (vector) from a dataframe using the established shortcut dataframe-name$column-name. For example:
```python
df$SampleID
df$PH
```

Lets also quickly cover some of tools we can use to test our conditions  
Relational Operators (Guess what they mean)  
* ==    
* \!=    
* <   
* >  
* <=  
* >=  


Logical Operators  
* &  
* \!  
* \|  


Only checks the first element  
* &&  
* ||  


OK, probably pretty confusing, but lets go over some examples that should make these concepts more clear.
```python
df$PH < 10
which(df$PH < 10)
#Whats the difference?

 
df[df$PH < 10, ]
df[which(df$PH < 10),]
#What is R doing in the above commands?


df$SampleID == "Sample13"
which(df$SampleID == "Sample13")
df[which(df$SampleID == "Sample13"),]

df$PH < 10
df$SampleID =="Sample13"

#Combining conditions
df$PH < 10 & df$SampleID == "Sample13"

#example for '&' vs '&&'
df$PH < 10 && df$SampleID == "Sample13" 
#checks to see if both boolean vectors are the same, returns single value

df$PH < 10 | df$SampleID == "Sample13"

#Think about what each of the following statements is doing. It can be helpful to translate these into English / German in your head or on paper.
df[df$DO > 0 & df$NH4 > 0.1,]
#For dataframe named df, return rows where the column DO values are greater than 0 AND where values from column NH4 are greater than 0.1; return all columns

df[which(df$DO > 0 & df$NH4 > 0.1),]
df[df$DO > 0 & df$NH4 < 0.2 & df$Season == "Summer",]
df[which(df$DO > 0 & df$NH4 < 0.2 & df$Season == "Summer"),]

#What does this command return? Can you translate it?
df[which(df$DO > 0 & df$NH4 < 0.2 | df$Season == "Summer"),]

#Some more complicated examples
df[which(df$PH <= 7.41 & df$TOC < 1.4 & df$TOC > 1.2),]
df[which(df$PH < 7.1 | df$PH > 8.1),]
```
Just as a comment, R is great because you can save a lot of the intermediate steps. If you find yourself getting lost in the statements, separate them out. Make sure that you reference the correct dataframe in each of the steps though!

```python
df_working <- df
df_working <- df_working[df_working$PH <= 7.41,]
#MAKE SURE YOU USE THE CORRECT REFERENCE!!!!

df_working <- df_working[df_working$TOC < 1.4,]
df_working <- df_working[df_working$TOC > 1.2,]
#Did this give you the same answer? How can you check?
```

## Using grep to find / search
grep is a command-line program that search for lines that match a regular expression. We can use the R implementation of grep to search vectors for elements that match a regular expression / keyword.

The syntax is grep(PATTERN, vector).

Lets look at some examples:
```python
#Search for values containing "HTU" in the vector df$Zone
df[grep("HTU", df$Zone),]

#Here is an example of a simple regular expression, it reads search for "HTU4" OR "HTU5" within the vector df$Zone
df[grep("HTU4|HTU5", df$Zone),]

#A more complicated regex, search for a string that ends with a two-digit number within the SampleID
df[grep(".*[0-9]{2}", df$SampleID),]
#The ".*" means match anything.

#What is this one doing?
df[c(grep(".*er", df$Season), grep("Autumn", df$Season)),]
```

A word of caution. I often find myself doing stupid stuff like:
```python
df[grep("Sample2[0-9]",df$SampleID),][df[grep("Sample2[0-9]",df$SampleID),]$DO < 0.1,]
```
Can you figure out what is going on? It is messy, ugly, and hard to read. Seeing code like this in your Rscripts is always sad since it will take you much longer to figure out what you were doing, and it will likely obscure bugs in your code.

It would be better to split your commands up:
```python
df_sample20_up <- df[grep("Sample2[0-9]",df$SampleID),]
df_sample20_up[df_sample20_up$DO < 0.1,]
```
&#49;&#54;. Get all data that was collected during June and July.

## Sorting data in R
Confusingly, you will almost always use the **order()** command to sort your data in R. Order gives you a vector with indices of the sorted data, while sort gives you the actual sorted vector.

For example:
```python
sort(df$SampleID)
order(df$SampleID)
```

In general, we use the vector that order returns to sort the full dataframe.

So we can do stuff like this:
```python
#Sort the dataframe based on the values in column "SampleID"
df[order(df$SampleID),]
#What happened here?

#Sort using a numeric value
df[order(row.names(df), decreasing = T),]
df[order(as.numeric(row.names(df)), decreasing = T),]

#If we didn't have a numeric column to sort on we can split apart a column to sort on
df[order(as.numeric(gsub("Sample", "", df$SampleID)), decreasing = T),]
#Here I am using gsub, the find & replace version of grep, to replace "Sample" with "" (nothing), then I'm telling R that these are numbers & not characters, finally I tell it to order these numbers from highest to lowest.

#Why does R sort these numerically by default?
df[order(df$PH),]

#You can check the type of vector you have with:
str(df$PH)
#Or to check all vectors in a dataframe you can use
str(df)
```

Sometimes you want to sort using multiple columns. Order allows you to specify multiple columns and it will sort your data in that order.

```python
df[order(df$PH, df$DO),]
df[order(df$PH, -df$DO),]
```
This question came up when I was working on this dataset and I think it is a fairly typical question.  
How can we order our dataset based on Months?

Lets check what we are starting with:
```python
df$Month
```
It seems like all are the common 3 letter abbreviations, except July **(There are ALWAYS exceptions in your data)**

Lets fix that:
```python
df[df$Month == "July",]$Month <- "Jul" #Gives error, missing values are not allowed
```
OK, that didn't work. What happened?  
Unfortunately it messed up the whole dataframe by adding NAs. Lets reload the data first.
```python
df <- read.delim(file = "simulated_dataset.txt", header=T)
```
So after googling, I relearned that changing the values of "factors" to a new value doesn't work well. Lets change months to a "character" data type.
```python
df$Month <- as.character(df$Month)
```
Now we can try to make the substitution:
```python
df[df$Month == "July",]$Month <- "Jul"
```
Success!!  

Now we need to tell R how we want the data organized. I couldn't think of a good way to order Months alphabetically in a manner that works. But luckily R has a defined vector of the common 3 letter abbreviations, and it is already sorted.

We can use this defined order, to structure our df$Month column. Here I am telling R to treat Months as a factor again, and then defining the order of the values of that factor. I know this is confusing, I typically google how to do this every time. It is extremely useful though, since many of the graphing packages order your data based on the factor levels.
```python
df$Month <- factor(df$Month, levels = month.abb)
```

Now we can order the full dataframe based on our Month factor:
```python
df[order(df$Month, decreasing = F),]
```
Yes, I still think it is MAGIC that this works. Seriously.

## Summary Questions
1. List the names for each of the data columns in df.
2. What is the name of the 3rd column in df.
3. Get the 9-12 rows and 4-10 columns from df.
4. Which samples from well H41 have TOC values less than 0.8?
5. Make a new dataframe that only contains samples collected from the "HTL" aquifer.
6. Using  this new dataframe, which samples had a higher temperature than 4.5?
7. Make another new dataframe that consists of 10 samples, 5 with the lowest NO3 and 5 with the highest NO3 measurements.
8. Re-order the original dataframe by PO4, with the lowest value at the top.
9. Re-order the original dataframe by PO4 and Mg, with the highest values at the top.
10. Re-order the original dataframe by PO4 with the highest values at the top, and then by Mg with the lowest values at the top.
11. Do the same as above, but reverse the order (Mg first, then PO4).
12. Rename the samples "AquaDivaN" when N is 1 to 30.
13. Make 4 new dataframes, 1 for each season.

14. Any other ways of manipulations data that you want / need?

## Manipulating and Summarizing our data with dplyr
dplyr is a very helpful package that is specifically designed to work on dataframes and work **QUICKLY**. It is the updated version of plyr.
I recommend checking out the following link for some ways to user dplyr.  
<https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html>

We'll go through the basic syntax dplyr uses and some examples.

```python
library(dplyr)
#reload the dataframe if necessary
df <- read.delim("simulated_dataset.txt", header=T)
```

### Filter the dataframe
Here we can grab the rows that match criteria we set. This is similar to the conditional operators we used in the above section.

Note: the %>% symbol can be read as a "pipe" where we are chaining together a series of commands.

```python
df %>% filter(Season == "Summer")
df %>% filter(Season =="Summer" & DO < 0.4)
df %>% filter(Well != "H41" & PH > 7.4 & DOC < 2)

filter(df, Season == "Summer")
```

### Sort the dataframe
We can use the arrange command in dplyr to sort our dataframe. This may be easier to remember than order.
```python
df %>% arrange(desc(DO), desc(NH4), Month)
```

### Select columns from the dataframe
```python
df %>% select(SampleID, DO, PH, Mg, Season)
df %>% select(SampleID, TEMP_W_ES:Na)
df %>% select(-(Well:Aquifer))

df %>% select(SampleID, starts_with("EC"))
```

### Add new columns based on functions of previous columns
Many times you will import a relatively simple dataset into R and then want to perform several calculations using that data. For instance you want to convert date / time column to "Hours since start", or you want to convert you absorbance data to counts based on a standard curve, etc...

Here we can use dplyr to run those calculations and add the result to a new column.  
This example will convert our DO measurements to percent oxygen saturation.
```python
df <- df %>% mutate(TEMP_K = TEMP_W_ES + 273.15)
#DO Saturation Formula assuming salinty is 0:
#  exp(-173.4292 + 249.6339*100/T + 143.3483*ln(T/100) + - 21.8492 * (T/100))

df <- df %>% mutate(DO_sat = exp(-173.4292 + (249.6339*100/TEMP_K) + (143.3483*log(TEMP_K/100)) + (-21.8492 * (TEMP_K/100))))
df <- df %>% mutate(O2_perc = DO / DO_sat * 100)

df %>% select(SampleID:Aquifer, O2_perc)
```

There is also a function in case we only care about the results and don't want to add it to the exisiting dataframe.
```python
#Save only new column
df %>% transmute(O2_perc = DO / DO_sat*100)
```

### Summarise / Summarize
This is usually how I end up using dplyr. It is an incredibly powerful way of quickly calculating summary statistics on your data while you try and explore what is going on. This would be analagous to using pivot tables in excel.

I find the dplyr syntax really nice for reading exactly what you are doing.

For example, we can calculate per Well statistics for a specific variable of interest.
```python
df %>% group_by(Well) %>% summarize(mean_O2_perc = mean(O2_perc), sd_o2_perc = sd(O2_perc))
#Take the dataframe df, split into smaller dataframes on the fly based on the values found in the df$Well column. Run the summarise() function on each of the these smaller dataframes to calculate the mean and standard deviation for O2 percent.
```
We can do something similar and calculate season based statistics:
```python
df %>% group_by(Season) %>% summarize(n = n(), mean_PH=mean(PH), mean_TOC=mean(TOC))
#Split dataframe df into smaller dataframes based on the values in the column "Season", the summarize these smaller dataframes by calculating the number of samples (n), and the mean for PH and TOC
```

Sometimes its useful to split the dataframes based on multiple values. Here we will generate summary statistics for each Well based on the sampling season.
```python
df %>% group_by(Well, Season) %>% summarize(n=n(), mean=mean(O2_perc))
df %>% group_by(Well, Season) %>% summarize(n=n(), mean=mean(O2_perc)) %>% arrange(mean, Well)
```

There are also some nice functions were we can generate summary statistics for all data columns at once, such as:
```python
df %>% select(Well, Season, TEMP_W_ES:Na) %>% group_by(Well, Season) %>% summarize_all(funs(mean))

df_sum <- df %>% group_by(Well, Season) %>% summarize_at(.vars = vars(TEMP_W_ES:Na), .funs = c(n="length", mean="mean", sd="sd"))

#I haven't figured out the best way to order these columns yet, but this is a quick and dirty way.
df_sum %>% select(Well, Season, order(names(df_sum)))
```

## dplyr Summary Questions
Questions 15-18:  
re-do questions 4, 5, 8, and 10 using dplyr functions  
&#49;&#57;. Find the mean SO4 levels for the two aquifers  
&#50;&#48;. Find the mean SO4 levels for the two aquifers each season  
&#50;&#49;. Add a new column with the Hydronium ion concentration for each sample (hint: H3O+ = 10^(-pH) )  
&#50;&#50;. Add a new column with the pOH value for each sample (pH + pOH = 14)  

<sub>Written by Will A. Overholt</sub>
<sub>Oct 2018</sub>
