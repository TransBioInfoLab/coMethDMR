# Gabriel Odom
# 2022-03-29

context("CpGsInfoAllRegions")


###  Setup  ###
data(betasChr22_df)
data(pheno_df)

AllRegionNames_char <- c("chr22:18267969-18268249", "chr22:18531243-18531447")


###  Actual  ###
out_df <- CpGsInfoAllRegions(
  AllRegionNames_char,
  betas_df = betasChr22_df,
  pheno_df = pheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex")
)


###  Expected  ###
expected_df <- data.frame(
  Region = c(
    rep("chr22:18267969-18268249", 4), rep("chr22:18531243-18531447", 3)
  ),
  cpg = c(
    "cg18370151", "cg12460175", "cg14086922", "cg21463605", "cg25257671",
    "cg06961233", "cg08819022"
  ),
  chr = factor(
    x = rep("chr22", 7),
    levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "*")
  ),
  pos = c(
    18267969L, 18268062L, 18268239L, 18268249L, 18531243L, 18531385L, 18531447L
  ),
  slopeEstimate = c(
    -0.0078, -0.0393, -0.0779, -0.1005, -0.0617, -0.0512, -0.1009
  ),
  slopePval = c(0.6911, 0.323, 0.0586, 0.0264, 0.1918, 0.2291, 0.1108),
  UCSC_RefGene_Name = rep("", 7),
  UCSC_RefGene_Accession = rep("", 7),
  UCSC_RefGene_Group = rep("", 7)
)


###  Tests  ###
test_that("CpGsInfoAllRegions() returns a data frame", {
  expect_equal(is.data.frame(out_df), TRUE)
})

test_that("CpGsInfoAllRegions() output has 9 columns", {
  expect_equal(ncol(out_df), 9)
})

test_that("CpGsInfoAllRegions() calculated estimates correctly", {
  expect_equal(out_df$slopeEstimate, expected_df$slopeEstimate)
})

test_that("CpGsInfoAllRegions() calculated p-values correctly", {
  expect_equal(out_df$slopePval, expected_df$slopePval)
})

