# Basic Statistical Tests with R
This is going to be fairly fast paced and brief discussion of many commonly used statistical tests and how to run them in R. I barely scratch the surface in terms of the kinds of tests available and which should be used. I **HIGHLY HIGHLY** recommend you check out the Handbook of Biological Statistics and the R companion book (both free online) when you are analyzing your own data!!

## Some extremely excellent statistic references
Handbook of Biological Statistics
<http://www.biostathandbook.com/>

R Companion for Biological Statistics
<https://rcompanion.org/rcompanion/a_02.html>

Intro R Statistics
<https://www.bioinformatics.babraham.ac.uk/training/R_Statistics/Introduction%20to%20Statistics%20with%20R%20manual.pdf>

Probability & Statistics
<https://cran.r-project.org/web/packages/IPSUR/vignettes/IPSUR.pdf>

Handbook of Statistical Analyses
<https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_introduction_to_R.pdf>

## Our plan today
Today we will discuss and then explore the following ideas:  
    - Common assumptions & what these look like.  
    - Students T tests and derivations  
    - ANOVA and derivations  
    - Linear Regression / Correlation  

## Assumptions of common parametric tests:
1. Independence  
The measurements collected are independent from each other.  
Common voiolations include:
    - timeseries data
    - repeated measurements from the same sample
    - measurements from small spatial region
    - Other?
  
2. Normality  
Do the measurements come from a normal or gaussian distribution?  
    - Surprisingly, this is actually not that important since most of the common parametric tests are fairly robust to deviations from normality
    - You can try to transform data to meet assumptions of normality
    - A great idea is to plot the data and see what you are getting yourself into

3. Homoscedasticity  
Or the assumption that different groups have the same standard deviation
    - more of a problem if you have an unbalanced dataset 
    - again best to plot you data and get a feel for what it looks likes

An example:
```python
#Loading the libraries we will be using today
library(ggplot2)
library(dplyr)

#An example dataset
ex_norm <- data.frame(sampleid = seq(1,100), 
                      y = rnorm(n=100, mean=10, sd=2), 
                      treatment=rep(c("A","B"), each=50))

#Our class dataset from yesterday
df <- df <- read.delim("simulated_dataset.txt", header=T)
```

Lets start exploring these data!

```python
#Plotting a histogram to see how your univariate data are arranged
ggplot(ex_norm, aes(y))+
  geom_histogram(binwidth = 0.5)

#Using a scatter plot to look at your data
    ggplot(ex_norm, aes(y))+
      geom_point(aes(x = treatment, y = y))

#Using a box plot
ggplot(ex_norm, aes(x=treatment, y=y))+
  geom_boxplot(aes(y=y))

#Fancier box plot (violin plot)
ggplot(ex_norm, aes(x=treatment, y=y))+
  geom_violin(aes(y=y))
```

Working with some of our data: The joys of real life
```python
ggplot(df, aes(NO3))+
  geom_histogram()

ggplot(df, aes(NO3))+
  geom_histogram(binwidth = 5)

ggplot(df, aes(NO3))+
  geom_histogram(aes(fill = Well))

ggplot(df, aes(x=Well, y = NO3))+
  geom_point()

ggplot(df, aes(x=Well, y = NO3))+
  geom_boxplot()
```
So what do you think?  
Are these groups normally distributed?  
We can use the Shapiro Wilk normality test to see. Just FYI, this test is quite sensitive to non-normal data. If it fails, you don't have to immediately panic since the following tests we'll talk about are much more robust to non-normal data.

```python  
df %>% select(Well,NO3) 
   %>% group_by(Well) 
   %>% summarise(statistic = shapiro.test(NO3)$statistic, 
                 p.value = shapiro.test(NO3)$p.value)

df %>% select(Well,NO3) 
   %>% group_by(Well) 
   %>% summarise(statistic = shapiro.test(sqrt(NO3))$statistic, 
                 p.value = shapiro.test(sqrt(NO3))$p.value)

#This test is running the following command on each of the Wells
shapiro.res <- shapiro.test(df[df$Well == "H43",]$NO3)
shapiro.res$statistic
shapiro.res$p.value

```
Is the sampling scheme balanced?
```python
df %>% select(Well,NO3) 
   %>%group_by(Well) 
   %>% summarise(n=n())
```

What can we do?
```python
#squareroot transform
ggplot(df, aes(x = Well, y = sqrt(NO3)))+
  geom_point()+
  stat_summary(fun.y = "mean", geom = "point", color="blue")

#log10 transform
ggplot(df, aes(x = Well, y = log10(NO3+0.5)))+
  geom_point()+
  stat_summary(fun.y = "mean", geom = "point", color="blue")
```

Which looks better? 
```python
#Note, this test is sensitive to non-normal data as well as heteroscedastic data.
#It is nice for picking the best transform though.
#raw data
bartlett.test(data=df, NO3 ~ Well)

#squareroot transform
bartlett.test(data=df, sqrt(NO3) ~ Well)

#log10 transform
bartlett.test(data=df, log10(NO3+1) ~ Well)
```
OK, it is your turn. Individually or in small groups, check other variables in the data frame
Ask yourself if:
    - normal
    - homoscedastic
    - independent
    - transform
Make sure you look at some of the variables associated with the aquifer.

## Student's t tests
### One-sample t tests
These are used to compare the mean of your population to a theoretical mean

```python
#Lets see if the Nitrate from well H43 has a mean of 0.
t.test(df[df$Well == "H43",]$NO3, mu = 0, conf.level = 0.95)

#We know that nitrate violates some assumptions of normality. 
#How does the more normal square root transformed data appear?
t.test(sqrt(df[df$Well == "H43",]$NO3), mu = 0, conf.level = 0.95)

#Always a good idea to see if you stats agree with your intuition of the data
ggplot(df[df$Well == "H43",], aes(x=Well, y=NO3))+
  geom_point(position = position_jitter())

#What happens with our example dataset?
t.test(ex_norm$y, mu = 10, conf.level = 0.95)
```
OK, now it is your turn, test to see if the lower aquifer (HTL) is anaerobic.


###Two-sample t tests
These are much more common (in my opinion). Here we are comparing the means between 2 groups.  
Lets see if the 1st half of our simulated normal data is significantly different from the last half.
```python
t.test(ex_norm[seq(1,50),]$y, ex_norm[seq(51,100),]$y)
#What is this actually doing?

#Lets plot our data again
ggplot(df, aes(x = Aquifer, y = NO3))+
  geom_boxplot()
#Check for equal variances
bartlett.test(data=df, NO3 ~ Aquifer)

#Run the two-sample t.test based on the aquifer
t.test(df[df$Aquifer == "HTL",]$NO3, df[df$Aquifer == "HTU",]$NO3)
```

In general though, if you have a large enough sample size (approximately > 5) you want to use a Welchs t test (default in R and what we've been using). This variation does not assume equal variance. In fact, many places recommend to always use this test over a students t test. If the samples have the same variance you can turn that on in this way
```python
t.test(df[df$Aquifer == "HTL",]$NO3, df[df$Aquifer == "HTU",]$NO3, var.equal = TRUE)
```

If you do not think you samples are normally distributed you can use a Mann Whitney U test. 
```python
wilcox.test(df[df$Aquifer == "HTL",]$NO3, df[df$Aquifer == "HTU",]$NO3)
#Note: We get the "cannot compute exact p-values since we have multiple measurements 
#of the same value (0) in the HTU dataset.
```
In practice, I almost always use the default Welch's t.test for un-equal variance.


## ANOVA
### One-way ANOVA
We can use an analysis of variance (ANOVA) test to compare the means of 2 or more groups. If it is only 2 groups, then a one-way anova and a two-sample t test are identical.

There are A LOT of different ways to do ANOVA tests and we will quickly cover a few. As with everything in this section, I highly recommend you double check references before you get started.

```python
no3_aov <- aov(data=df, sqrt(NO3) ~ Well)
summary(no3_aov)
```
Understanding the output:  
ANOVA essentially tests if the variance between your group means is different from the variance within the groups. Under the null hypothesis, all the groups have a same mean, so the variance between those group means will be the same as the average within group variation. There are 100s of books written on ANOVA and I recommend checking some of them out. Might also be good to manually calculate the values yourself once if you are interested.  
The ratio of these variances under the null fits an F-distribution (hence the F value). The p.value can then be determined from how different our tests F value is from the expected F value under the null hypothesis.  

The residuals in an R anova table refer to the within group statistics. The degrees of freedom are "number observations" - "number of groups". The group variable is specified (between-group means) and the degrees of freedom here are calculated as the "number of groups" - 1. 

Ideally you report the results as Nitrate concentrations in the wells were significantly different (one-way ANOVA, F<sub>2,27</sub> = 29, p < 1e-7).  

Note, ANOVA assumes your group data are normally distributed, with equal standard deviations in the groups, and are the measurements are independent. However ANOVA is fairly robust to non-normal data and if your sampling scheme is balanced it is robust to differences in standard deviations between the groups as well. If your experimental plan is not balanced, you need to be more careful as you can get a lot more false positives (p values < 0.05, but no actual differences between the means). That said, if you have a very significant p value (ours is 1.71e-7) your means are probably different.  
A non-parametric version of ANOVA can also be used if your data are not normal.
This would be the Kruskal Wallis test. It will have a lot less power and should be picked if your data are not normal AND you dont expect them to be normal. Interestingly, some people (John McDonald of the Handbook for Biological Statistics) recommend never using this test.
```python
kruskal.test(data=df, sqrt(NO3) ~ Well)
```
There is also a version of ANOVA if you think your data are normal, but you do not expect equal variances.
```python
oneway.test(data=df, sqrt(NO3) ~ Well)
```
There are a couple of built in ways to check your assumptions after running the ANOVA test.
```python
#Check afterward your ANOVA for whether the groups have similar variances
#Here we are plotting the residuals vs the group means. 
#Remember that the residuals are the "data values" - "group mean". 
#There should be no relationship between the residuals and the means (flat line). 
#The spread should be approximately equal as well.
plot(no3_aov,1)

#We can also check to see if the group data are normal.
plot(no3_aov, 2)
#Here we are plotting quantiles of the residuals against the quantiles of a normal distribution.
#Normal data will match the theoretical normal distribution giving a straight 1:1 line.
```

OK great, now we are pretty sure that NO3 values are different between our wells. However, we want to know **WHICH** wells are different from each other.  
Remember, ANOVA only tells you if there is a diffence in your groups means and NOT which groups are different.  

We can do post-hoc tests to see which groups are significantly different. The most common option that I've seen is to use Tukeys HSD (honest significant differences) test. This assumes equal variances between the groups, and corrects for multiple tests.
```python
TukeyHSD(no3_aov)
```
There are no base R functions for unbalanced / unequal variances for post-hoc testing. There are other packages you can try though (multcomp, lsmeans, DescTools, agricolae).

```python
#You can also do two-way t tests and correct the p.values for multiple tests.
t.test(sqrt(df[df$Well == "H52",]$NO3), sqrt(df[df$Well == "H41",]$NO3))$p.value
t.test(sqrt(df[df$Well == "H52",]$NO3), sqrt(df[df$Well == "H43",]$NO3))$p.value
t.test(sqrt(df[df$Well == "H43",]$NO3), sqrt(df[df$Well == "H41",]$NO3))$p.value
p.adjust(c(0.008654862, 0.00662454, 1.876758e-10), method="bonferroni")
```

Pick a different variable, and test the means. Make sure to check how badly it violates the assumptions.

### Two-way ANOVA
You can use a two-way ANOVA when there are two parameters (factors) associated with your data. In our case, an example would be "Well" and "Season". We have sampled our 3 wells multiple times throughout the year and over multiple years. We can use a two-way ANOVA to test whether our measured data is significantly different between well, between season, and if there is a well specific seasonal effect (intereaction between our factors).  

Two-way ANOVA is often used for repeated measurements from the same individuals and you want to control for differences between those individuals.  

Two-ANOVA with replication tests 3 NULL hypotheses. 1. The means within one of our factors are the same, 2. The means within the other factor are the same. 3. There is no interaction between the factors. You can also do a two-way ANOVA test without replication (one data point per group), but you need to assume that there is no interaction term and the test is much weaker.

Generally, if the interaction between your variables is significant then you need to change your analysis plan since it is not fair to test for significance in only one of the variables. Instead, you can split your data up and do 2 one-way ANOVA tests. 

Lets take a look for what we are testing:
```python
ggplot(df)+
  geom_boxplot(aes(x = Well, y = sqrt(NO3), color = Season))
ggplot(df)+
  geom_boxplot(aes(x = Season, y = sqrt(NO3), color = Well))
```
Thoughts?

```python
#Run the anova in R and save the result
no3_2fac_aov <- aov(data=df, sqrt(NO3) ~ Season + Well + Season:Well)

#view information about the ANOVA
summary(no3_2fac_aov)

#Example of how R specifies ANOVA formulas
summary(aov(data=df, sqrt(NO3) ~ Season + Well))
summary(aov(data=df, sqrt(NO3) ~ Season*Well))
```

Lets look at the assumptions:
```python
plot(no3_2fac_aov, 1)
plot(no3_2fac_aov, 2)
#Actually looks quite a bit better than I expected...
#Might be a good idea to remove the outliers and re-test though.
```

OK, so this is not really fair. But as an example (assuming no interaction was found) we can still do follow up testing between groups.
```python
TukeyHSD(no3_2fac_aov, which = "Well")
```

Take some time and run a two-way ANOVA test with a different dependent variable.

## Linear Regression and Correlation
Linear regression and correlation are used to investigate the relationship between two continuous variables. In general you want to know if two variables are related, how closely they are related, and mathematically describe this relationship.

I am not going into details into the differences between linear regression and correlation. Check out this page for more info http://www.biostathandbook.com/linearregression.html

Also, R is perfectly happy to let you run linear regression with 1 measurement variable and factors (discrete / categorical variables). In this case, it actually runs ANOVA/MANOVA in the background. I find this both convenient (you can use the same syntax and still be correct) and confusing (you may be running a different test than you expect).

This is also why you see a lot of examples where people set up their ANOVA test with a linear regression model.
For example:
```python
linear_model <- lm(data=df, sqrt(NO3) ~ Well)
summary(aov(linear_model))
```

Anways, lets see if Calcium and Sodium concentrations are related.
```
ggplot(df, aes(x = Ca, y = Na))+
  geom_point()
```
This looks pretty clear, but lets test it.
By default, R uses a pearson correlation test which assumes that the data are linearly related and that the residuals are normally distributed. The null hypothesis is that the variables are not related (slope = 0). 
```python
cor.test(~ Ca + Na, data=df)
```
A nonparametric correlation test commonly used is the Spearman rank correlation test. It does not assume a distribution for the variables of if they are related. 
```python
cor.test(~ Ca + Na, data=df, method="spearman")
```
The correlation coefficient and the p value tend to be reported. The correlation coefficient (rho) goes from -1 to 1, where -1 indicates the variables are perfectly negatively correlated, and 1 indicates they are perfectly positively correlated. The p value is calculated using a t distribution with n-2 degrees of freedom.


If we want to generate a linear regression model for these measurements, we can use the built in function lm().
```python
lm(Ca ~ Na, data = df)
```
The formula format R uses is Y ~ X1  

The "Y" variable is traditionally assigned as the dependent variable while the "X" variable(s) are the independent ones.
Our example is not great since we are looking at two dependent variables that are probably not independent. In R the lm() function fits the ordinary-least-squares regression line, or the line that minimizes the distance between your observations and the line itself. This works best for analysing two continuous variables with both independent and dependent variables. It is not valid for 2 dependent variables (assuming there is no cause and effect relationship). There are other techniques in this instance, however, they are not valid if you want to predict unmeasured values.

Anyways, I strongly recommend you read the literature or use methods typically employed in your field before running these analyses on your own data.

But for us, lets assume everything is fine and explore our linear regression results.

```python
linear_model1 <- lm(Ca ~ Na, data = df)
summary(linear_model1)
#Note that the pvalue & the t.value for the slope are the same as we got from the Pearson correlation analysis.
cor.test(~ Ca + Na, data=df)
```

The r^2 value (coefficient of determination) is a measure of how well your regression line fits the data, or how much of the variance of your Y variable (dependent) is explained by the X variable. Values close to 1 indicate your observed Y values fall on your regression line while values close to 0 indicate there is no relationship between your variables.

Like with ANOVA, we can explore data characteristics from the linear regression results.
```python
#Check for homoscedascity & bias (line should be flat around 0)   
plot(linear_model1,1)

#Check for normality
plot(linear_model1,2)
```
OK, lets try another example.
We know that oxygen concentrations are physically dependent on temperature. Lets see if there is a relationship between these values.

Which is the dependent and which is the independent variable?

```python  
ggplot(df, aes(x = TEMP_W_ES, y=DO))+
  geom_point()
cor.test(~ TEMP_W_ES + DO, data=df)
summary(lm(DO ~ TEMP_W_ES, data=df))
```
But, wait. Lets think about this? Does a value of 0 depend on temperature in this case?
 
```python
ggplot(df, aes(x = TEMP_W_ES, y=DO))+
  geom_point(aes(color=Well))

#Lets look at only wells with oxygen.
cor.test(~ TEMP_W_ES + DO, data=df[df$Well == "H41",])
lm_o2 <- lm(data=df[df$Well == "H41",], DO ~ TEMP_W_ES)
summary(lm_o2)
plot(lm_o2, 1)
plot(lm_o2, 2)
```

A very quick aside and primer on multiple linear regression (univariate).  
Here we expect multiple continuous independent variables to effect our measured dependent variable. For example, maybe both temperature and PH affect DOC concentrations. If your goal is predictive and you have a bunch of independent variables that you would like to use to predict unmeasured dependent variables this can be a powerful method. If you are trying to explain cause & effect, you need to be very careful.  
In general, your X variables should not be correlated. If they are, then it is difficult to parse out the effect (which one matters). If you just want to predict a value then it doenst matter, but adding highly correlated variables wont really improve your predictions.
This is a huge and actively growing area of statistics (along with multivariate multiple regression) so I will just talk about how we could do simple multiple regression in R.
Also, linear regression **does not** show causation. It might provide supportive evidence, but you need to be very careful with interpreting your data.

```python
#Lets check how our dependent variable is related to our independent variables
ggplot(df, aes(x = TEMP_W_ES, y = DOC))+
  geom_point()
ggplot(df, aes(x = PH, y = DOC))+
  geom_point()
#Lets see if our independent variables are correlated
cor.test(data=df, ~ PH + TEMP_W_ES)

#Setting up the linear model:
lm_DOC <- lm(data=df, DOC ~ TEMP_W_ES*PH)
#What does the * mean here?
  
#Results from the model:
summary(lm_DOC)
```

## Parting Words
To wrap this up. I hope you're now convinced of how useful this free programing language is. We did not touch upon many of the incredible pacakges that have been developed. Suffice to say, I have never found a problem that I had in R that someone else hadn't already solved and posted aobut (although sometimes it takes awhile to figure out how to ask Google the right question).  
I know I'm sounding like a broken record, but it is always a good idea to familarize yourself with the fundamentals when you are performing some data analysis and I've founded the references at the top to be extremely useful.  

Also, I have had extremely good luck with stackedoverflow posts when I'm stuck or trying to figure out new methods. I encourage you to search the internet for examples of problems that you have. Start simple and refine your analysis as you go.  
But remember, it is best to consider how you will answer your research question **before** you design the experiment. It's impossible to retroactively fix a poorly designed experimental plan and from personal experience spending more time at the beginning stages will make all your future analysis much easier.

<sub>Written by Will A. Overholt</sub>
<sub>Oct 2018</sub>

