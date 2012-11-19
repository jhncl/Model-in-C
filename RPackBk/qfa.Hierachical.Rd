\name{qfa.Hierachical}
\alias{qfa.Hierachical}
\title{Fit a hierachical logistic curve model}
\description{
Fits a hierachical logistic curve model. Custom models may be used
}
\usage{
qfa.Hierachical(experiment,Scaling,iter,upd,thin,inits,PlotOutput=TRUE,work,CustomModel=FALSE)
}
\arguments{
  \item{experiment}{A data set obtained from Rod.Read() for one experiment.}
  \item{Scaling}{Do you require the data set to be scaled from (0-1)? Are your priors set for a scaled data set? TRUE or FALSE.}
  \item{iter}{The number of iterations to output for the posterior sample.}
  \item{upd}{The number of iterations used for updating the model.}
  \item{thin}{Every number of iterations is stored from the posterior sample. Thinning the number of iterations of the posterior sample down from the choosen number iter.}
  \item{PlotOutput}{TRUE or FALSE argument, whether to create plots that are then stored to .pdf files in the working directory.}
  \item{work}{A unique name to identify the functions output.}
  \item{CustomModel}{To fit a custom model, set as custom model file name, otherwise CustomModel should remain defined as FALSE.}
  \item{inits}{!!!!!!!}
}
\value{
An R list that includes the posterior sample means, the data set... This can be then fed into the functions qfaplots.H() or qfa.Interaction().
}
\keyword{qfaBayes}
