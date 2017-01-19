###################################################
##
## Test the methylInheritanceSimInternalMethods functions
##
###################################################


###################################################
## estBetaAlpha() functions
###################################################

test.estBetaAlpha_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlpha(c(0.5, 0.2), 0.00001)
    
    exp <- 0.125
    
    message <- paste0("test.estBetaAlpha_good_01() ",
                    "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## estBetaBeta() functions
###################################################

test.estBetaBeta_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBeta(c(0.3, 0.2), 0.00001)
    
    exp <- 0.035
    
    message <- paste0("test.estBetaBeta_good_01() ",
                      "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}