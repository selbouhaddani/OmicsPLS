library(O2PLS)
context("O2PLS main function checking")

test_that("Normal input goes without error", {
  expect_error(o2m(1:10, 1:10, 1, 0, 0), NA)
  expect_error(o2m(diag(1,6,5), diag(1,6,6), 1, 1, 1), NA)
  expect_error(o2m(diag(1,6,5), diag(1,6,6), 1, 1, 1, stripped = TRUE), NA)
  expect_error(o2m(diag(1,6,5), diag(1,6,6), 1, 1, 1, p_thresh = 1, q_thresh = 1), NA)
  expect_error(o2m(diag(1,6,5), diag(1,6,6), 1, 1, 1, p_thresh = 1, q_thresh = 1, stripped = TRUE), NA)
})

test_that("Errors in o2m are thrown", {
  expect_error(o2m(diag(3),diag(4),1,0,0),                   "rows")
  expect_error(o2m(matrix(1:9,3),matrix(1:9,3),0,1,1),       "n")
  expect_error(o2m(matrix(c(1:8,NA),3),matrix(1:9,3),1,1,1), "NA")
  expect_error(o2m(matrix(1:9,3),matrix(c(1:8,NaN),3),1,1,1),"NaN")
  expect_error(o2m(matrix(1:9,3),matrix(c(1:8,Inf),3),1,1,1),"finite")
  expect_error(o2m(diag(4),diag(4),1.5,0,0),                 "n")
  expect_error(o2m(diag(4),diag(4),1,1.5,0),                 "nx")
  expect_error(o2m(diag(4),diag(4),1,0,1.5),                 "ny")
})

test_that("size, ratios and names are correct", {
  expect_equal(o2m(diag(4),diag(4),1,1,0)$R2X,      0.5)
  expect_equal(o2m(diag(4),diag(4),1,0,1)$R2Y,      0.5)
  expect_equal(nrow(o2m(diag(4),diag(4),1,0,0)$W.), 4)
  expect_equal(nrow(o2m(diag(5),diag(5),1,0,0)$C.), 5)
  expect_equal(rownames(o2m(data.frame(a=1:10,b=2:11,c=3:12),orth(1:10),1,0,0)$W.), letters[1:3])
  expect_equal(rownames(o2m(orth(1:10),data.frame(a=1:10,b=2:11,c=3:12),1,0,0)$C.), letters[1:3])
})