\name{getcij}
\alias{getcij}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate chemospecificity matrix
%%  ~~function to do ... ~~
}
\description{
This function converts the raw gradients values of the input into the chemospecificity matrix.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getcij(chem, Ai, aj, Bi, bj, Nr, Nt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chem}{
%%     ~~Describe \code{chem} here~~
}
  \item{Ai}{
%%     ~~Describe \code{Ai} here~~
}
  \item{aj}{
%%     ~~Describe \code{aj} here~~
}
  \item{Bi}{
%%     ~~Describe \code{Bi} here~~
}
  \item{bj}{
%%     ~~Describe \code{bj} here~~
}
  \item{Nr}{
%%     ~~Describe \code{Nr} here~~
}
  \item{Nt}{
%%     ~~Describe \code{Nt} here~~
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
function(chem,Ai,aj,Bi,bj,Nr,Nt){
	EPSA=1
	EPSB=1
	EPS=1
	
	cij<-mat.or.vec(Nr,Nt)
	if(chem=="add"){
		for(i in 1:Nr){
			for(j in 1:Nt){
				cij[i,j]=EPSA*Ai[i]*(max(aj)-aj[j])+EPSB*Bi[i]*bj[j]
			}
		}
	}else if(chem=="mult"){
		for(i in 1:Nr){
			for(j in 1:Nt){
				cij[i,j]=EPS*Ai[i]*(max(aj)-aj[j])*Bi[i]*bj[j]
			}			
		}
	}
	cij=1+cij
	return(cij)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
