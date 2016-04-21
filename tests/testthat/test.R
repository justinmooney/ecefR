library(ecefR)
context("Number Checking")

test_that("ECEF2LLA is correct", {
    ecef <- c(-1472.238, 5005.461, 3692.592)
    lla <- c(35.46943, 106.39, 20.9998)
    expect_equal(xyz2lla(ecef), lla, tolerance = 1e-6)
    expect_equal(lla2xyz(lla), ecef, tolerance = 1e-6)
})