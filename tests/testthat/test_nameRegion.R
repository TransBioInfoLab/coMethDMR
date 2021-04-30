# Lizhong Liu
context("nameRegion")

CpGsOrdered_df <- data.frame(
  chr = c("chr10", "chr10", "chr10", "chr10"),
  pos = c(100028236L, 100028320L, 100028468L, 100028499L),
  cpg = c("cg20214853", "cg04677227", "cg11632906", "cg07146435")
)

test_that("nameRegion returns vector with correct value", {

  expect_equal(
    object = NameRegion(CpGsOrdered_df),
    expected = "chr10:100028236-100028499"
  )

})

