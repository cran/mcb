\name{mcb.compare}
\alias{mcb.compare}
\title{Comparisons of Model Confidence Bounds for Different Variable selection Methods}
\description{
This function is a supplement of the function mcb. It is used to compare different variable selection methods and would return all the MUCs on same canvas. A good variable selection method’s MUC will tend to arch towards the upper left corner.
}
\usage{
mcb.compare(x, y, B=200, lambdas=NA, methods=NA, level=0.95, seed=122)
}
\arguments{
\item{x}{input matrix presenting independent variables as in mcb.}
\item{y}{response vector as in mcb.}
\item{B}{number of bootstrap replicates to perform; Default value is 200.}
\item{lambdas}{A vector of penalty tuning parameters for each variable selection method. The default values are the optimal choices for each selection method computed automatically.}
\item{methods}{a vector including all variable selection methods the user wants to test and compare. The default value is c ('aLasso', 'Lasso', 'SCAD', 'MCP', 'stepwise', 'LAD', 'SQRT')}
\item{level}{user-defined confidence level as in mcb; Default value is 0.95.}
\item{seed}{Default value is 122.}
}
\value{
The mcb.compare method returns an object of class “mcb.compare”
An object of class "mcb.compare " is a list containing at least the following components:
\item{mcb}{a list containing the bootstrap coverage rate and the corresponding model confidence bound for all user-given variable selection methods.}
\item{mucplot}{plot of the model uncertainty curves for all variable selection methods and could be used to choose the best method.}
\item{mcbframe}{a list containing all the information about MCBs for all variable selection methods under all available bootstrap coverage rates.}
}
\references{
Li,Y., Luo,Y., Ferrari,D., Hu,X. and Qin,Y. (2019) Model Confidence Bounds for Variable Selection. Biometrics, 75:392-403.
}
\examples{
data(Diabetes) # load data
x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')]
y <- Diabetes[,c('Y')]
x <- data.matrix(x)
y <- data.matrix(y)
\dontrun{result <- mcb.compare(x=x, y=y)}
# plot of the model uncertainty curves for all variable selection methods
\dontrun{result$mucplot}
# a list containing the bootstrap coverage rate and mcb which based on Lasso
\dontrun{result$mcb$Lasso}
# a dataframe containing all the information about MCBs which based on Lasso
\dontrun{result$mcbframe$Lasso}
}
