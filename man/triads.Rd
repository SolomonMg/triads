\name{triadcensus}
\alias{triadcensus}
\title{ Node-level triad census for networks }
\description{ 
  \code{triadcensus} calculates the triad census positions of 
 *each node* in a graph object. Currently, it works via an igraph
 object, utilizing the igraph function get.adjlist() to return the data
 in the format necessary to input into Moody's functions.  
}

\usage{
triadcensus(g)
}

\arguments{
  \item{g}{ igraph graph object }
}
\value{
  \item{triads}{a matrix of triad membership}}

\author{Solomon Messing <messing~at~stanford...>}

\note{ These functions replicate SAS module itpcen.mod, programmed
 by James Moody (last updated 2004). The program was written for 
 Moody's dissertation: The Structure of Adolescent Social Relations: 
 Modeling Friendship in Dynamic Social Settings. University of 
 North Carolina at Chapel Hill, Department of Sociology. 
 Peter S. Bearman, Chair. This is based on routines originally 
 developed by Ron Burt and implemented in his "STRUCTURE" program.
 A link to the manual is available at \url{http://faculty.chicagobooth.edu/ronald.burt/teaching/STRUCmanual.pdf}.}
\examples{
##

  n <- 20
  g <- barabasi.game(n, power = 1, m = NULL) 
  triadcensus <- triadcensus(g)
}
\keyword{ ~sna }

