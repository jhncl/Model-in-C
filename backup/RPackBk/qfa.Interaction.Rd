\name{qfa.Interaction}
\alias{qfa.Interaction}
\title{Fit an interaction model}
\description{
Fits an interaction model. Custom models may be used
}
\usage{
qfa.Interaction(Control,Query,iter,upd,thin,PlotOutput=TRUE,work,CustomModel=FALSE,Priors=FALSE,CustomInteractionDef=FALSE)
}
\arguments{
  \item{Control}{A data set obtained from Rod.Read() for the Control experiment.}
  \item{Query}{A data set obtained from Rod.Read() for the Control experiment.}
  \item{iter}{The number of iterations to output for the posterior sample.}
  \item{upd}{The number of iterations used for updating the model.}
  \item{thin}{Every number of iterations is stored from the posterior sample. Thinning the number of iterations of the posterior sample down from the choosen number iter.}
  \item{PlotOutput}{TRUE or FALSE argument, whether to create plots that are then stored to .pdf files in the working directory.}
  \item{work}{A unique name to identify the functions output.}
  \item{CustomModel}{To fit a custom model, set as custom model file name, otherwise CustomModel should remain defined as FALSE.}
  \item{CustomInteractionDef}{To fit a custom interaction defintion, set as custom interaction defintion file name, otherwise CustomInteractionDef should remain defined as FALSE.}
  \item{Priors}{!!!!!!!!}
}
\value{
An R list that includes the posterior sample means, the data set... This can be then fed into the functions qfaplots.I().
}
\keyword{qfaBayes}
