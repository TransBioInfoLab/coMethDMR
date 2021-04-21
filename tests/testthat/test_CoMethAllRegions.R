# Lizhong Liu
context("CoMethAllRegions")

data(betasChr22_df)
CpGsChr22_ls <- readRDS(
  system.file ("extdata",
               "CpGislandsChr22_ex.RDS",
               package = 'coMethDMR',
               mustWork = TRUE
  )
)

out_ls <- CoMethAllRegions(
  dnam = betasChr22_df,
  CpGs_ls = CpGsChr22_ls,
  arrayType = "450k",
  returnAllCpGs = FALSE
)

test_that("CoMethAllRegions() returns a list", {
  expect_equal(is.list(out_ls), TRUE)
})

test_that("CoMethAllRegions() returns an S3 list", {
  expect_type(out_ls, "list")
})

test_that("CoMethAllRegions() has 3 elements", {
  expect_equal(length(out_ls), 3)
})

