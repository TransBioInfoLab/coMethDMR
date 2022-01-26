# Test Preliminary Data Flags
# Gabriel Odom
# 2022-01-26

context("MarkMissing")

data(betasChr22_df)
test_df <- betasChr22_df


###  Add two bad rows  ###
P <- ncol(test_df)
badProbes_mat <- rbind(
  matrix(
    rep(NA_real_, P), nrow = 1
  ),
  matrix(
    c( runif(n = P - 10) , rep(NA_real_, 10) ),
    nrow = 1
  )
)
rownames(badProbes_mat) <- c("cgTest1", "cgTest2")
colnames(badProbes_mat) <- colnames(test_df)

test_df <- rbind(betasChr22_df, badProbes_mat)


###  Add two bad columns  ###
N <- nrow(test_df)
badSamps_df <- data.frame(
  GSMtest1 = rep(NA_real_, N),
  GSMtest2 = c( runif(N - 6500), rep(NA_real_, 6500) )
)
test_df <- cbind(test_df, badSamps_df)


###  Execute Tests  ###
mark1_ls <- MarkMissing(test_df)
test_that(
  "Probes are dropped correctly",
  expect_equal(
    object = mark1_ls$dropProbes,
    expected = c("cgTest1", "cgTest2")
  )
)
test_that(
  "Samples are dropped correctly",
  expect_equal(
    object = mark1_ls$dropSamples,
    expected = c("GSMtest1", "GSMtest2")
  )
)

