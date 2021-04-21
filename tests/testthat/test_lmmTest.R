# Lizhong Liu
context("lmmTest")

data(betasChr22_df)
data(pheno_df)

CpGsChr22_char <- c(
  "cg02953382", "cg12419862", "cg24565820", "cg04234412", "cg04824771",
  "cg09033563", "cg10150615", "cg18538332", "cg20007245", "cg23131131",
  "cg25703541"
)
coMethCpGs <- CoMethSingleRegion(CpGsChr22_char, betasChr22_df)
coMethBetaMatrix <- betasChr22_df[coMethCpGs$CpGsSubregions[[1]], ]

test_that("lmmTest returns df with correct classes", {

  expect_s3_class(
    lmmTest(
      betaOne_df = coMethBetaMatrix,
      pheno_df,
      contPheno_char = "stage",
      covariates_char = c("age.brain", "sex"),
      modelType = "randCoef",
      arrayType = "450k"
    ),
    "data.frame"
  )

})

