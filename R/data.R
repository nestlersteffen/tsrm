#' Zero acquaintance round robin data from Richard Rau
#'
#' Data from a round robin study conducted by Richard Rau at the University of Leipzig in 2017.
#'
#' @format
#' A data frame with 2,162 measurements from one large round-robin group on the following 4 
#' round robin variables (taken on unnumbered 6-point rating scales with higher numbers indicating a
#' higher values): \cr
#' \code{liking}: rating of dimension liking \cr
#' \code{dominant}: rating of dimension dominant \cr
#' \code{affectionate}: rating of dimension affectionate \cr
#' \code{knowing}: rating of whether person i knows j \cr
#'
"Leipzig"

#' Zero acquaintance round robin data from David Kenny
#'
#' Data from a Albright et al. (1988) Study 2
#'
#' @format
#' A data frame with 124 measurements from 7 round-robin groups on the following 5 round robin
#' variables (taken on unnumbered 7-point rating scales with higher numbers indicating a
#' higher value of the trait): \cr
#' \code{sociable}: rating of dimension sociable \cr
#' \code{irritable}: rating of dimension good-natured \cr
#' \code{responsible}: rating of dimension responsible \cr
#' \code{anxious}: rating of dimension calm \cr
#' \code{intellectual}: rating of dimension intellectual \cr
#' The data frame also contains the gender (\code{actor.sex}; \code{1} = F,
#' \code{2} = M) of the participants and their self-ratings on the five assessed traits
#' (\code{actor.sociable} and so on).
#'
#' @references
#' Albright, L., Kenny, D. A., & Malloy, T. E. (1988). Consensus in 
#' personality judgments at zero acquaintance. \emph{Journal of Personality 
#' and Social Psychology}, \emph{55}(3), 387–395. 
#' \doi{10.1037/0022-3514.55.3.387}
#'
#'@source
#' \url{http://davidakenny.net/srm/srmdata.htm}
#'
"Kenzer"

#' Round robin data reported in Warner et al. (1979)
#'
#' Data from a Warner et al. (1979) Table 7
#'
#' @format
#' A data frame with 56 measurements from one round-robin group. Each variable represents the
#' proportion of time speaking by person i to j on one of three days.
#'
#' @references
#' Warner, R. M., Kenny, D. A., & Stoto, M. (1979). A new round robin 
#' analysis of variance for social interaction data. \emph{Journal of Personality 
#' and Social Psychology}, \emph{37}(10), 1742–1757.
#' \doi{10.1037/0022-3514.37.10.1742}
#'
"Warner"

#' Zero acquaintance round robin data from Thomas Malloy
#'
#' Data from Albright et al. (1988) Study 1
#'
#' @format
#' A data frame with 216 measurements from 12 round-robin groups on the following 5 round-robin
#' variables (taken on unnumbered 7-point rating scales with higher numbers indicating a
#' higher value of the trait with the exception for good and calm): \cr
#' \code{sociable}: rating of dimension sociable \cr
#' \code{irritable}: rating of dimension good-natured \cr
#' \code{responsible}: rating of dimension responsible \cr
#' \code{anxious}: rating of dimension calm \cr
#' \code{intellectual}: rating of dimension intellectual \cr
#'
#' The data frame also contains the gender (\code{actor.sex}; \code{1} = F,
#' \code{2} = M) of the participants and their self-ratings on the five assessed traits
#' (\code{actor.sociable} and so on).
#'
#' @references
#' Albright, L., Kenny, D. A., & Malloy, T. E. (1988). Consensus in 
#' personality judgments at zero acquaintance. \emph{Journal of Personality 
#' and Social Psychology}, \emph{55}(3), 387–395. 
#' \doi{10.1037/0022-3514.55.3.387}
#'
#'@source
#' \url{http://davidakenny.net/srm/srmdata.htm}
#'
"Malzer"

#' Round robin data reported in Lashley and Bond (1997) 
#'
#' Data from Lashley and Bond (1997) Table 1
#'
#' @format
#' A data frame with 30 measurements from one round-robin group. A single number is the
#' number of times an actor initiated aggression toward a partner.
#'
#' @references
#' Lashley, B. R., & Bond, C. F., Jr. (1997). Significance testing for round robin data. 
#' \emph{Psychological Methods}, \emph{2}(3), 278–291. 
#' \doi{10.1037/1082-989X.2.3.278}
#'
"LashleyBond"

#' Round robin data from Hallmark and Kenny
#'
#' Data from Kenny et al. (1994)
#'
#' @format
#' A data frame with 802 measurements from 30 round-robin groups on the following 7 round-robin
#' variables (taken on unnumbered 7-point rating scales with higher numbers indicating a
#' higher value of the trait): \cr
#'
#' \code{calm}: rating of dimension calm-anxious \cr
#' \code{sociable}: rating of dimension sociable-withdrawn \cr
#' \code{liking}: rating of dimension like-do not like \cr
#' \code{careful}: rating of dimension careful-careless \cr
#' \code{relaxed}: rating of dimension relaxed-tense \cr
#' \code{talkative}: rating of dimension talkative-quiet \cr
#' \code{responsible}: rating of dimension responsible-undependable \cr
#'
#' The data frame also contains the gender (\code{actor.sex}; \code{1} = F,
#' \code{2} = M) of the participants and their age in years (\code{actor.age}).
#' Note that the data was assessed in two conditions: odd round robin group numbers indicate
#' groups in which participants rated all traits for a person at a time whereas even numbers
#' refer to groups in which participants rated all the people for each trait.
#'
#' @references
#' Kenny, D. A., Albright, L., Malloy, T. E., & Kashy, D. A. (1994). Consensus in
#' interpersonal perception:  Acquaintance and the big five.
#' \emph{Psychological Bulletin, 116}(2), 245-258.
#' \doi{10.1037/0033-2909.116.2.245}
#'
#'@source
#' \url{http://davidakenny.net/srm/srmdata.htm}
#'
"HallmarkKenny"

#' Triadic round robin data reported in Bond, Horn, and Kenny (1997) 
#'
#' Data from Bond, Horn, and Kenny (1997) Table 3
#'
#' @format
#' A data frame with 120 measurements from one round-robin group. A single number is a
#' judge's perception of an actor's liking of a partner. j.id is the id of the judge, p.id and 
#' t.id are the ids of the perceiver (actor) and target (partner), respectively.
#'
#' @references
#' Bond, C. F., Jr., Horn, E. M., & Kenny, D. A. (1997). A model for triadic relations. 
#' \emph{Psychological Methods}, \emph{2}(1), 79–94. 
#' \doi{10.1037/1082-989X.2.1.79}
#'
"BondHornKenny"

#' Triadic round robin data from Curry and Emerson (1970) 
#'
#' Data from a study by Curry and Emerson (1970) used in Kenny et al. (1996)
#'
#' @format
#' A data frame with 2,688 measurements from eight round-robin groups. A single number is a
#' judge's perception of an actor's liking of a partner. 
#'
#' @references
#' Kenny, D. A., Bond, C. F., Mohr, C. D., & Horn, E. M. (1996). Do we now how much people like
#' another? \emph{Journal of Personality and Social Psychology}, \emph{71}, 928–936. 
#' \doi{10.1037/0022-3514.71.5.928}
#'
"CurryEmerson"