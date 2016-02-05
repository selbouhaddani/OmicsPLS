library(O2PLS)
context("O2PLS main function checking")

test_that("Normal input goes without error", {
  expect_error(o2m(diag(5), diag(5), 1, 1, 1), NA)
  expect_error(o2m(diag(5), diag(5), 1, 1, 1, stripped = TRUE), NA)
  expect_error(o2m(diag(5), diag(5), 1, 1, 1, p_thresh = 1, q_thresh = 1), NA)
  expect_error(o2m(diag(5), diag(5), 1, 1, 1, p_thresh = 1, q_thresh = 1, stripped = TRUE), NA)
})

test_that("Errors in o2m are thrown", {
  expect_error(o2m(diag(3),diag(4),1,0,0),                   "nrow")
  expect_error(o2m(matrix(1:9,3),matrix(1:9,3),0,1,1),       "joint")
  expect_error(o2m(matrix(c(1:8,NA),3),matrix(1:9,3),0,1,1), "is.na")
  expect_error(o2m(matrix(1:9,3),matrix(c(1:8,NaN),3),0,1,1),"is.na")
  expect_error(o2m(matrix(1:9,3),matrix(c(1:8,Inf),3),0,1,1),"is.finite")
  expect_error(o2m(diag(4),diag(4),1.5,0,0),                 "round")
  expect_error(o2m(diag(4),diag(4),1,1.5,0),                 "round")
  expect_error(o2m(diag(4),diag(4),1,0,1.5),                 "round")
})

test_that("size and ratios are correct", {
  expect_equal(o2m(diag(4),diag(4),1,1,0)$R2X,      0.5)
  expect_equal(o2m(diag(4),diag(4),1,0,1)$R2Y,      0.5)
  expect_equal(nrow(o2m(diag(4),diag(4),1,0,0)$W.), 4)
  expect_equal(nrow(o2m(diag(5),diag(5),1,0,0)$C.), 5)
})
