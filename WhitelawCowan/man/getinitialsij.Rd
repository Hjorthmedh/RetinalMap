\name{getinitialsij}
\alias{getinitialsij}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate initial conditions
%%  ~~function to do ... ~~
}
\description{
This function returns the matrix of inital conditions for the system
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getinitialsij(Nr, Nt, type, smallnumber = 1, min = 5e-05, max = 0.00011)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nr}{
%%     ~~Describe \code{Nr} here~~
}
  \item{Nt}{
%%     ~~Describe \code{Nt} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{smallnumber}{
%%     ~~Describe \code{smallnumber} here~~
}
  \item{min}{
%%     ~~Describe \code{min} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Nr,Nt,type,smallnumber=1,min=0.00005,max=0.00011){
	switch(type,
		zero=mat.or.vec(Nr,Nt),
		small=smallnumber*(mat.or.vec(Nr,Nt)+1),
		rand=matrix(runif(Nr*Nt,min,max),nrow=Nr)
	)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
