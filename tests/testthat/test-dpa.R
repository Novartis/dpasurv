test_that("dpa works correctly with bootstrap", {
  data(simdata)

  # results of dynamic path analysis
  set.seed(1234)
  res_dpa <- dpa(survival::Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=100)

  # test to verify length and names
  expect_equal(length(res_dpa), 4)
  expect_equal(names(res_dpa), c("coefs","boot.coefs","meta","aalen"))

  expect_equal(names(res_dpa$coefs), c("outcome","M"))
  expect_equal(names(res_dpa$coefs$outcome), c("times","Intercept","M","x"))
  expect_equal(names(res_dpa$coefs$M), c("times","(Intercept)","x"))

  expect_equal(dim(res_dpa$coefs$outcome)[1],length(unique(sort(simdata[["stop"]][simdata[["event"]]==1]))))
  expect_equal(dim(res_dpa$coefs$M)[1],length(unique(sort(simdata[["stop"]][simdata[["event"]]==1]))))

  expect_equal(names(res_dpa$boot.coefs), c("outcome","M"))
  expect_equal(names(res_dpa$boot.coefs$outcome), c("times","Intercept","M","x","boot.id"))
  expect_equal(length(unique(res_dpa$boot.coefs$outcome$boot.id)), 100)

  expect_equal(names(res_dpa$boot.coefs$M), c("times","(Intercept)","x","boot.id"))
  expect_equal(length(unique(res_dpa$boot.coefs$M$boot.id)), 100)

  expect_equal(names(res_dpa$meta), c("outcome","mediator","variables"))

})
