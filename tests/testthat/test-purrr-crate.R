testthat::test_that("crated purrr function does not throw error", {
  fn = carrier::crate(function(x) return(paste0("test", x)))
  # try(pak::pkg_remove("mirai"), silent = TRUE)
  testthat::expect_equal(
    purrr::map_chr(1:5, fn),
    c("test1", "test2", "test3", "test4", "test5")
  )
})
