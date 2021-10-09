# Lizhong Liu
context("lmmTestAllRegions")

data(betasChr22_df)
data(pheno_df)

# CpGisland_ls <- readRDS(
#   system.file(
#     "extdata", "CpGislandsChr22_ex.rds",
#     package = 'coMethDMR',
#     mustWork = TRUE
#   )
# )
# 
# coMeth_ls <- CoMethAllRegions(
#   dnam = betasChr22_df,
#   rDropThresh_num = 0.5,
#   CpGs_ls = CpGisland_ls,
#   arrayType = "450k",
#   returnAllCpGs = FALSE
# )

coMeth_ls <- list(
  `chr22:18268062-18268249` = c("cg12460175", "cg14086922", "cg21463605"),
  `chr22:18531243-18531447` = c("cg25257671", "cg06961233", "cg08819022")
)

test_that("lmmTestAllRegions returns df with correct classes", {

  expect_s3_class(
    lmmTestAllRegions(
      betas = betasChr22_df,
      region_ls = coMeth_ls,
      pheno_df,
      contPheno_char = "stage",
      covariates_char = c("age.brain", "sex"),
      modelType = "randCoef",
      arrayType = "450k"
    ),
    "data.frame"
  )

})

