# Issue 7
# Gabriel Odom and Fernanda Veitzman
# 2022-01-21

# We received an issue on GitHub:
#   https://github.com/TransBioInfoLab/coMethDMR/issues/7
# We resolved it, but I need to confirm that our solution is correct and write
#   a testing script for it. Also, this is probably related to Issue #11 because
#   it's caused by problems in the input data.


###  Current example  ###
data(betaMatrix_ex1)
CreateRdrop(data = betaMatrix_ex1, method = "pearson")


###  Add All Missing Sample  ###
beta2_mat <- rbind(betaMatrix_ex1, rep(NA_real_, ncol(betaMatrix_ex1)))
CreateRdrop(data = beta2_mat, method = "pearson")
# No error; no warning
CreateRdrop(data = beta2_mat, method = "pearson", use = "everything")
# All r_drop values = 0; Warning message: "Missing correlation values detected.
#   These are set to 0."


###  Add All Missing Probe  ###
beta3_mat <- cbind(betaMatrix_ex1, cgTest = rep(NA_real_, nrow(betaMatrix_ex1)))
CreateRdrop(data = beta3_mat, method = "pearson")
# Error in cor(data_i, data_no_i_mean, method = method, use = use) :
#   no complete element pairs
CreateRdrop(data = beta3_mat, method = "pearson", use = "everything")
# All r_drop values = 0; Warning message: "Missing correlation values detected.
#   These are set to 0."


###  Add Partially Missing Sample  ###
beta4_mat <- rbind(
  betaMatrix_ex1,
  c(
    NA_real_,
    runif(n = ncol(betaMatrix_ex1) - 1)
  )
)
CreateRdrop(data = beta4_mat, method = "pearson")
# No error; no warning
CreateRdrop(data = beta4_mat, method = "pearson", use = "everything")
# All r_drop values = 0; Warning message: "Missing correlation values detected.
#   These are set to 0."


###  Add Partially Missing Probe  ###
beta5_mat <- cbind(
  betaMatrix_ex1,
  cgTest = c(
    NA_real_,
    runif(n = nrow(betaMatrix_ex1) - 1)
  )
)
CreateRdrop(data = beta5_mat, method = "pearson")
# No error; no warning
CreateRdrop(data = beta5_mat, method = "pearson", use = "everything")
# All r_drop values = 0; Warning message: "Missing correlation values detected.
#   These are set to 0."
