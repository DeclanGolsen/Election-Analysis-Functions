# Election Analysis Function #1: Proportions Approximation Function
## A first draft of an R function that approximates solutions to a system of linear equations composed entirely of proportions
### Introduction
As a statistics major, I took pretty much every opportunity I was given to analyze election results and regional demographics,
two fields where proportions play a significant role. The idea of being able to easily solve systems of equations composed entirely
of proportions has seemed interesting to me since I first took linear algebra, though for a long time I didn't know exactly how to
go about it. Since the proportions I was interested in were often rounded percentages, traditional techniques for solving a system of
linear equations proved ineffective. 

As an example of the type of linear equation I'm talking about, take a theoretical country composed of three regions, where three
languages are spoken. The following table provides the proportion of residents in each region that speak each language:

| **Region** | Region 1 | Region 2 | Region 3 | Total |
|------------|----------|----------|----------|-------|
| Language 1 | 0.57     | 0.40     | 0.36     | 0.40  |
| Language 2 | 0.23     | 0.41     | 0.11     | 0.25  |
| Language 3 | 0.20     | 0.19     | 0.53     | 0.35  |
| Population | ?        | ?        | ?        | 1.00  |

We have a linguistic breakdown by region and nationally, but we don't have the weight of each region; in other words, we have the
following system of linear equations:

$0.57x_1 + 0.40x_2 + 0.36x_3 = 0.40$

$0.23x_1 + 0.41x_2 + 0.11x_3 = 0.25$

$0.20x_1 + 0.19x_2 + 0.53x_3 = 0.35$

Where $x_1$ is the percent of the population living in Region 1, $x_2$ is the percent of the population living in Region 2, and $x_3$ is the percent of the population living in Region 3. If no rounding was involved, it would be possible to find a unique solution to this system of linear equations; since these proportions are rounded, however, this proves to be impossible, and this problem stymied me for quite a while. After taking a class involving the use of Monte Carlo methods, however, I was inspired to instead go about *approximating* the solutions to such systems, using iterations wherein random values for $x_1$, $x_2$, ..., $x_n$ are chosen and a distance formula is used to calculate the accuracy of that particular iteration. In order to do this, I wrote the following R function:

``` R
    proportions_approx <- function(x,y,n=1000){
  d <- data.frame(matrix(0,nrow=n,ncol=(ncol(x)+1)))
  for (m in (1:n)) {
    remain <- 1
    xn <- data.frame(matrix(0,nrow=ncol(x),ncol=1))
    for (j in (1:(ncol(x)-1))) {
      if (remain >= xn[j,1]) {
        xn[j,1] <- runif(1,0,remain)
        d[m,(1+j)] <- xn[j,1]
        remain <- remain - xn[j,1]
      }
      else {
        break
      }
    }
    xn[ncol(x),1]<-remain
    d[m,(ncol(x)+1)]<-remain
    xn
    dm <- data.frame(matrix(0,nrow=nrow(x),ncol=1))
    for (i in (1:nrow(x))) {
      di<-0
      dn<-0
      for (j in (1:ncol(x))) {
        di<-x[i,j]*xn[j,1]
        dn<-dn+di
      }
      dm[i,1] <- (y[i,1] - dn)^2
    }
    d[m,1] <- sum(dm)
  }
  print(d[which(d$X1==min(d$X1)),])
}
```
This function takes a matrix (*x*), a column matrix (*y*), and an optional integer (*n*) as inputs. The *x* matrix corresponds to
proportions of the value of interest in each region; in the example provided above, that would be:

$$\begin{array}{ccc}
0.57 & 0.40 & 0.36\\
0.23 & 0.41 & 0.11\\
0.20 & 0.19 & 0.53\\
\end{array}$$

The *y* column matrix corresponds to results across categories, i.e. the overall proportions for the value of interest. In the given
example this would be:

$$\begin{array}{c}
0.40\\
0.25\\
0.35\\
\end{array}$$

Finally, *n* corresponds to the number of iterations to occur when approximating solutions. The default quantity is 1,000, though 5,000
or even 10,000 iterations ought to be used when there are a large number of columns in the *x* matrix (regions in the example). Once
these inputs have been provided, the function assigns random values to $x_1$, $x_2$, ..., $x_n$ that add up to 1 for each iteration, then and stores the distance between the predicted value in the corresponding *y* matrix and the actual corresponding value in the
*y* matrix. In this example, if the values chosen for  $x_1$, $x_2$, and $x_3$ were 0.36, 0.21, and 0.43, respectively, the predicted
values of the *y* matrix would be 0.44, 0.22, and 0.34, and thus the distance (as defined in this function) between the predicted and actual values would be
$(0.44-0.40)^2$ + $(0.22-0.25)^2$ + $(0.34-0.35)^2$ = 0.0026. The function then finds the iteration that yielded the smallest distance,
and outputs that distance, and then the assigned values of $x_1$, $x_2$, ..., $x_n$. To use the function on the example data, we would
type the following into R:

``` R
set.seed(300)
x<-t(matrix(c(0.57, 0.40, 0.36, 0.23, 0.41, 0.11, 0.20, 0.19, 0.53),nrow=3,ncol=3))
y<-matrix(c(0.40, 0.25, 0.35),nrow=3)
proportions_approx(x=x,y=y)
```

When the seed is set at 300, the 367th iteration is deemed the best, where $x_1$, $x_2$, and $x_3$ are set as 0.1161024, 0.4142341, and
0.4696635, respectively, and the distance is 4.851886 x $10^{-6}$. This approximation is quite close to the true values used to generate the
*y* matrix, which were $x_1$ = 0.13, $x_2$ = 0.40, and $x_3$ = 0.47, and if more iterations were used, the approximation would be even
closer.

### Using real-world examples

To test the predictive ability of this function, I looked at some data from the 2020 US Presidential Election, namely at results by county in Delaware and Nevada. I elected to look at Delaware first since it contains only three counties, and could provide basic information on the effectiveness of the function.

For Delaware, I used the data in the following table (taken from the Wikipedia page for the [2020 United States Presidential Election in Delaware](https://en.wikipedia.org/wiki/2020_United_States_presidential_election_in_Delaware)) to create the *x* and *y* matrices:

| **County** | Kent   | New Castle | Sussex | Total  |
|------------|--------|------------|--------|--------|
| Democrat   | 0.5120 | 0.6781     | 0.4382 | 0.5878 |
| Republican | 0.4712 | 0.3072     | 0.5507 | 0.3980 |
| Other      | 0.0168 | 0.0147     | 0.0111 | 0.0142 |
| Population | ?      | ?          | ?      | 1.00   |

And used the following code to create a prediction:

``` R
set.seed(300)
xDel<-matrix(c(0.5120, 0.4712, 0.0168, 0.6781, 0.3072, 0.0147,0.4382, 0.5507, 0.0111),nrow=3,ncol=3)
yDel<-matrix(c(0.5878, 0.3980, 0.0142),nrow=3)
proportions_approx(xDel,yDel)
```

Using one thousand iterations, the values for $x_1$, $x_2$, and $x_3$ that resulted in the smallest distance were 0.0483, 0.6103, and 0.3414, respectively; in other words, the function predicted that 4.83% of Delaware voters lived in Kent County, 61.03% lived in New Castle County, and 34.14% lived in Sussex County. Ordinally, these values are accurate, but this prediction has significant room for improvement, as seen in the table below.

| **County**         | Kent       | New Castle | Sussex     | Total    |
|--------------------|------------|------------|------------|----------|
| Democrat           | 0.5120     | 0.6781     | 0.4382     | 0.5878   |
| Republican         | 0.4712     | 0.3072     | 0.5507     | 0.3980   |
| Other              | 0.0168     | 0.0147     | 0.0111     | 0.0142   |
| **Predicted 1000** | **0.0483** | **0.6103** | **0.3414** | **1.00** |
| Actual Pop.        | 0.1727     | 0.5707     | 0.2566     | 1.00     |

Most notably, the population of Kent County was significantly underpredicted, while the population of Sussex County was overpredicted. Using 5,000 iterations instead, however, the prediction is rather accurate, as seen in the table below the code:

``` R
set.seed(300)
proportions_approx(xDel,yDel,5000)
```
| **County**         | Kent       | New Castle | Sussex     | Total    |
|--------------------|------------|------------|------------|----------|
| Democrat           | 0.5120     | 0.6781     | 0.4382     | 0.5878   |
| Republican         | 0.4712     | 0.3072     | 0.5507     | 0.3980   |
| Other              | 0.0168     | 0.0147     | 0.0111     | 0.0142   |
| Predicted 1000     | 0.0483     | 0.6103     | 0.3414     | 1.00     |
| **Predicted 5000** | **0.1904** | **0.5666** | **0.2430** | **1.00** |
| Actual Pop.        | 0.1727     | 0.5707     | 0.2566     | 1.00     |

Accordingly it seems that, with a large enough sample size and using the first four digits of all the proportions involved, a reasonable prediction for the weight of each category can be achieved with this function.

I chose to also look at Nevada's counties, given that the state has a greater number of them (16 counties and 1 independent city), and because the population of Nevada is extremely consentrated in two counties. Analyzing the state ought to allow the accuracy of the function to be determined when using bigger matrices, as well as accuracy across larger and smaller categories. Given the size of the *x* matrix, I used 5,000 iterations off the bat for Nevada's counties. I also used a CSV instead of inputing data by hand, taken from the Wikipedia page for the [2020 United States Presidential Election in Nevada](https://en.wikipedia.org/wiki/2020_United_States_presidential_election_in_Nevada). My results were determined using the following code:

```R
nevada <- read.csv("~/Documents/nevada.csv", header=FALSE)
set.seed(300)
yNV<-t(nevada[18,])
xNV<-t(nevada1)
proportions_approx(xNV,yNV,5000)
```

Since the table of results is pretty big, I’ll describe some of the findings instead. The proportion of the voting population found within Clark County, home to Las Vegas, was predicted almost perfectly (the prediction was 0.6949, while the actual proportion was 0.6920). For the remaining sixteen counties and county-equivalents, however, the predictions seem to have been rather random; Carson City was predicted to have a proportion of 0.2836, whereas in actuality it was 0.0212, and Washoe County was predicted to have a proportion of 1.772 x $10^{-11}$, when in actuality it was 0.1794. These inaccuracies are almost certainly the results of the randomization technique I employed. The method will be explained further in the **Areas for improvement** section, but essentially the expected value of a randomized proportion decreases the further down its category is in the list. For this dataset counties are listed alphabetically, and Washoe county is the second to last; this means the expected value of the randomized proportion assigned to it is 0.1526 x $10^{-5}$, so it is extremely unlikely that the actual proportion, 0.1794, would be assigned to it, and it would take tens of thousands of iterations to reasonably believe that this prediction would be generated. I also tried out the function with 10,000 iterations, using the following code:

```R
set.seed(300)
proportions_approx(xNV,yNV,10000)
```
And the result was less accurate than the 5,000 iteration version, assigning the proportion 0.8654 to Clark County. Thus, for an *x* matrix with a lot of columns, this version of the function is only useful for identifying larger trends. After the randomization technique is updated, hopefully the function will be more useful for larger datasets.

### Applications

As mentioned in the introduction section, this function can be used for any system of linear equations composed entirely of proportions. Thus this function could be used on any table with missing categorical weights and that shows proportions, whether they pertain to election results, demographics, virology, or any other subject. For hypothetical tables of proportions, where the weights are chosen last, this function would allow for the generation of reasonable approximations, and would also provide some insight as to how feasible it is that there would be any solution to the system of linear equations. Looking back at the language table involving a theoretical country, for example, say the *y* column matrix was rearranged and the table instead looked like this:

| **Region** | Region 1 | Region 2 | Region 3 | Total |
|------------|----------|----------|----------|-------|
| Language 1 | 0.57     | 0.40     | 0.36     | 0.25  |
| Language 2 | 0.23     | 0.41     | 0.11     | 0.40  |
| Language 3 | 0.20     | 0.19     | 0.53     | 0.35  |
| Population | ?        | ?        | ?        | 1.00  |

The total difference between Language 2 and Language 1 is given as 0.15, but Language 2 is only more common than Language 1 in one region, and only by a single percentage point. Thus the function will provide the best solution it can, but the distance, the first number reported by the function, will reflect the accuracy of the solution provided. This should be taken into account in all uses of the function, but especially when the existence of a solution is uncertain.

### Areas for improvement

As mentioned in the **Using real-world examples** section, the most significant area of improvement for this function is the randomization technique. Each proposed value of $x_i$ is uniformly distributed, but the upward bound of the range of each distribution is lowered by draws for earlier values; for example, if 0.43 is selected for $x_1$, and 0.31 is selected for $x_2$, then $x_3$ can only exist between 0 and 0.26. Accordingly, the expected value of each $x_i$ corresponds to $0.5^i$, since the expected value of a uniformly distributed variable between 0 and 1 is 0.5. This method ensures that the values for $x_1$, $x_2$, ..., $x_n$ add up to 1, but means that categories that appear later in the matrix are much less likely to be accurately predicted. I believe the solution to this problem lies in randomizing the assignment of predictions for each iteration; this way the order of rows in a matrix would not impact the prediction for each row. The expected value of draws would still diminish rather quickly, but at least as far as US Counties go, it often is the case that a handful of counties are home to most a state’s population, while most others have a very small population. Prior information about the sort of distribution the weights of categories follows could inform the sort of randomization technique employed, but in any case the one currently in use will be changed at some point in the near future.

Another potential area of improvement would be allowing for incomplete information regarding category weight to included, i.e. if the weight of the population of one region was known, but that of three others was not. In any case, an updated version of the function will likely be created within the next few months.

### Conclusion

The proportions approximation function, proportions_approx() in R, accurately approximates the solutions to a system of linear equations if all values of the *x* matrix and the *y* column matrix are proportions whose columns add up to 1, and if there are relatively few columns in the *x* matrix (say five or less). When there are more columns, the function will be likely to pick up on general trends (i.e. whether a category accounts for a large or small percent of the total), but not more subtle ones. Once the randomization technique is improved upon the function will make better predictions for larger *x* matrices, but accuracy will still be lost as the number of columns increases. Thus, the proportions approximation function has both immediate applications and the prospect of serious improvement.
