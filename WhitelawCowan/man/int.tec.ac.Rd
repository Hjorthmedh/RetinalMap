\name{int.tec.ac}
\alias{int.tec.ac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inter-tecal-connections
%%  ~~function to do ... ~~
}
\description{
This function returns the matrix of inter-tectal connections
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
int.tec.ac(Nt, dist.sc, bij, dist.bij)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nt}{
%%     ~~Describe \code{Nt} here~~
}
  \item{dist.sc}{
%%     ~~Describe \code{dist.sc} here~~
}
  \item{bij}{
%%     ~~Describe \code{bij} here~~
}
  \item{dist.bij}{
%%     ~~Describe \code{dist.bij} here~~
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
function(Nt,dist.sc,bij,dist.bij){
	mat.bij<-mat.or.vec(Nt,Nt)
	q<-which(dist.sc<=dist.bij,arr.ind=TRUE)
	mat.bij[q]=1
	return(mat.bij)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
