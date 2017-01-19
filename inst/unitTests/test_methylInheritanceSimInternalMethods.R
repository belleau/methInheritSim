###################################################
##
## Test the methylInheritanceSimInternalMethods functions
##
###################################################


###################################################
## estBetaAlpha() functions
###################################################

test.estBetaAlpha_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlpha(c(0.5,0.1, 0.22, 0.21))
    
    exp <- 0.75
    
    message <- paste0("test.estBetaAlpha_good_01() ",
                    "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

