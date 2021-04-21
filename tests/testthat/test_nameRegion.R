# Lizhong Liu
context("nameRegion")

CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
CpGsOrdered_df <- OrderCpGsByLocation(CpGs_char,
                                      genome = "hg19",
                                      arrayType = c("EPIC"),
                                      output = "dataframe")

test_that("nameRegion returns vector with correct value", {

  expect_equal(
    object = NameRegion(CpGsOrdered_df),
    expected = "chr10:100028236-100028499"
  )

})

