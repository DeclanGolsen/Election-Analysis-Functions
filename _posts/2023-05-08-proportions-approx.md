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
these inputs have been provided, the function assigns random values to $x_1$, $x_2$, ..., $x_n$ that add up to 1 for each iteration, then calculates and
calculates and stores the distance between the predicted value in the corresponding *y* matrix and the actual corresponding value in the
*y* matrix. In this example, if the values chosen for  $x_1$, $x_2$, and $x_3$ were 0.36, 0.21, and 0.43, respectiveely, the predicted
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
0.4696635, respectively, and the distance is $4.851886x10^-6$. This approximation is quite close to the true values used to generate the
*y* matrix, which were $x_1$ = 0.13, $x_2$ = 0.43, and $x_3$ = 0.47, and if more iterations were used, the approximation would be even
closer.

### Using real-world examples

To be added soon!
