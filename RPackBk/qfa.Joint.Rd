\name{qfa.Joint}
\alias{qfa.Joint}
\title{Fit a joint hierachical logistic curve model}
\description{
Fits a joint hierachical logistic curve model. Custom models may be used
}
\usage{
qfa.Joint(Control,Query,Scaling,iter,upd,thin,PlotOutput=TRUE,work,CustomModel=FALSE)
}

\arguments{
  \item{Control}{A data set obtained from Rod.Read() for the Control experiment.}
  \item{Query}{A data set obtained from Rod.Read() for the Control experiment.}
  \item{Scaling}{Do you require the data set to be scaled from (0-1)? Are your priors set for a scaled data set? TRUE or FALSE.}
  \item{iter}{The number of iterations to output for the posterior sample.}
  \item{upd}{The number of iterations used for updating the model.}
  \item{thin}{Every number of iterations is stored from the posterior sample. Thinning the number of iterations of the posterior sample down from the choosen number iter.}
  \item{PlotOutput}{TRUE or FALSE argument, whether to create plots that are then stored to .pdf files in the working directory.}
  \item{work}{A unique name to identify the functions output.}
  \item{CustomModel}{To fit a custom model, set as custom model file name, otherwise CustomModel should remain defined as FALSE.}
}
\value{
An R list that includes the posterior sample means, the data set... This can be then fed into the functions qfaplots.I().
}
\keyword{qfaBayes}
