# Introductory notes ------------------------------------------------------

# Analytic-GRIMMER (A-GRIMMER) was developed by Aurélien Allard
# (https://aurelienallard.netlify.app/post/anaytic-grimmer-possibility-standard-deviations/).
# His original algorithm received some modifications here, for these reasons:
# -- Tapping scrutiny's infrastructure for implementing error detection
# techniques; for example, functions like `reround()` and
# `decimal_places_scalar()`.
# -- Changing the return value to logical, which is the expected output from the
# basic implementation of any consistency test within scrutiny.
# -- Adjusting variable names to the tidyverse style guide and scrutiny's
# domain-specific conventions.
# -- Adding support for multi-item scales (see `items` argument) via
# `rsprite2::GRIMMER_test()`.

# Translation of variable names -------------------------------------------

# original           --> scrutiny
# ********               ********
#
# aGrimmer           --> grimmer_scalar
# mean               --> x
# SD                 --> sd
# decimals_mean      --> digits_x  (argument removed; counting internally)
# decimals_SD        --> digits_sd (argument removed; counting internally)
# realmean           --> x_real
# realsum            --> sum_real
# effective_n        --> n_items
# Lsigma             --> sd_lower
# Usigma             --> sd_upper
# Lowerbound         --> sum_squares_lower
# Upperbound         --> sum_squares_upper
# Possible_Integers  --> integers_possible
# Predicted_Variance --> var_predicted
# Predicted_SD       --> sd_predicted
# Matches_Oddness    --> matches_parity
# FirstTest          --> pass_test1
# Matches_SD         --> matches_sd (to which a `pass_test2` object was added)
# Third_Test         --> pass_test3

# # Example inputs 1:
# x <- "1.03"
# sd <- "0.41"
# n <- 40
# items <- 1
# show_reason <- TRUE
# rounding <- "up_or_down"
# threshold <- 5
# symmetric <- FALSE
# tolerance <- .Machine$double.eps^0.5
# # decimals_mean <- 2
# # decimals_SD <- 2

# # Example inputs 2:
# # (actually derived from this distribution: c(1, 1, 2, 3, 3, 4, 4, 4, 4, 5))
# x <- "3.10"
# sd <- "1.37"
# n <- 10
# items <- 1
# show_reason <- TRUE
# rounding <- "up_or_down"
# threshold <- 5
# symmetric <- FALSE
# tolerance <- .Machine$double.eps^0.5
# # decimals_mean <- 2
# # decimals_SD <- 2

# # Example inputs 3:
# # (edge case from `pigs5`)
# x <- "2.57"
# sd <- "2.57"
# n <- 30
# items <- 1
# show_reason <- TRUE
# rounding <- "up_or_down"
# threshold <- 5
# symmetric <- FALSE
# tolerance <- .Machine$double.eps^0.5
# # decimals_mean <- 2
# # decimals_SD <- 2

# Implementation ----------------------------------------------------------

grimmer_scalar <- function(
  x,
  sd,
  n,
  items = 1,
  show_reason = FALSE,
  rounding = "up_or_down",
  threshold = 5,
  symmetric = FALSE,
  tolerance = .Machine$double.eps^0.5
) {
  check_type(x, "character")
  check_type(sd, "character")

  # # Provisional warning about the test 3 false-positive bug
  # cli::cli_warn(c(
  #   "False-positive GRIMMER results possible.",
  #   "!" = "I became aware of a bug in the `grimmer*()` functions.",
  #   "x" = "GRIMMER's test 3 can flag consistent values as inconsistent.",
  #   "x" = "For now, please use `show_reason = TRUE` and interpret results \
  #     of test 3 with care. (The first two tests and GRIM are not affected.)",
  #   "i" = "The next version of scrutiny will provide a fix. Many apologies.",
  #   "i" = "For more information, see \
  #   {.url https://github.com/lhdjung/scrutiny/issues/80}"
  # ))

  x_orig <- x

  # put GRIM TEST forward, since it doesn't require any of the variables computed for GRIMMER
  # GRIM TEST: It says `x_orig` because the `x` object has been coerced from
  # character to numeric, but `grim_scalar()` needs the original number-string.
  # Similarly, since this function also gets `items` passed down, it needs the
  # original `n`, not `n_items`.
  pass_grim <- grim_scalar(
    x = x_orig,
    n = n,
    items = items,
    rounding = rounding,
    threshold = threshold,
    symmetric = symmetric,
    tolerance = tolerance
  )

  if (!pass_grim) {
    if (show_reason) {
      return(list(FALSE, "GRIM inconsistent"))
    }
    return(FALSE)
  }

  digits_x <- decimal_places_scalar(x)
  digits_sd <- decimal_places_scalar(sd)
  x <- as.numeric(x)
  sd <- as.numeric(sd)

  n_items <- n * items

  p10x <- 10^(digits_x + 1)
  p10x_frac <- 5 / p10x
  # x (mean) bounds, lower and upper:
  if (x < p10x_frac) {
    x_lower <- 0
  } else {
    x_lower <- x - p10x_frac
  }
  x_upper <- x + p10x_frac

  p10 <- 10^(digits_sd + 1)
  p10_frac <- 5 / p10
  # SD bounds, lower and upper:
  if (sd < p10_frac) {
    sd_lower <- 0
  } else {
    sd_lower <- sd - p10_frac
  }
  sd_upper <- sd + p10_frac

  # create a vector of all possible integers between the lower and upper bounds of the sum:
  sum_lower <- x_lower * n_items
  sum_upper <- x_upper * n_items
  sum_possible <- round(sum_lower):round(sum_upper)
  x_real_possible <- sum_possible / n_items

  test1_results <- logical(length(x_real_possible))
  test2_results <- logical(length(x_real_possible))
  test3_results <- logical(length(x_real_possible))
  for (ii in seq_along(x_real_possible)) {
    x_real <- x_real_possible[ii]
    sum_real <- x_real * n_items
    

    # Sum of squares bounds, lower and upper:
    sum_squares_lower <- ((n - 1) * sd_lower^2 + n * x_real^2) * items^2
    sum_squares_upper <- ((n - 1) * sd_upper^2 + n * x_real^2) * items^2

    # TEST 1: Check that there is at least one integer between the lower and upper
    # bounds (of the reconstructed sum of squares of the -- most likely unknown --
    # values for which `x` was reported as a mean). Ceiling the lower bound and
    # flooring the upper bound determines whether there are any integers between
    # the two. For example:
    # -- If `sum_squares_lower` is 112.869 and `sum_squares_upper` is 113.1156,
    # `ceiling(sum_squares_lower)` and `floor(sum_squares_upper)` both return
    # `113`, so there is an integer between them, and `<=` returns `TRUE`.
    # -- TODO: add an example where there is no integer in between; and thus,
    # `pass_test1` is `FALSE`.
    test1_results[ii] <- ceiling(sum_squares_lower) <= floor(sum_squares_upper)
    # if test1_results[ii] is FALSE, no need to run the rest of the tests, fill test2_results[ii] and test3_results[ii] with FALSE and continue to the next iteration of the loop
    if (!test1_results[ii]) {
      test2_results[ii] <- FALSE
      test3_results[ii] <- FALSE
      next
    }
    

    # Create a vector of all possible integers between the lower and upper bounds
    # of the sum of squares:
    integers_possible <- ceiling(sum_squares_lower):floor(sum_squares_upper)

    # Create the predicted variance and SD:
    var_predicted <- (integers_possible / items^2 - n * x_real^2) / (n - 1)
    sd_predicted <- sqrt(var_predicted)

    # Reconstruct the SD:
    sd_rec_rounded <- reround(
      x = sd_predicted,
      digits = digits_sd,
      rounding = rounding,
      threshold = threshold,
      symmetric = symmetric
    )

    # Introduce a small numeric tolerance to the reported and reconstructed SD
    # values before comparing them. This helps avoid false-negative results of the
    # comparison (i.e., treating equal values as unequal) that might occur due to
    # spurious precision in floating-point numbers.
    sd <- dustify(sd)
    sd_rec_rounded <- dustify(sd_rec_rounded)

    # Check the reported SD for near-equality with the reconstructed SD values;
    # i.e., equality within a very small tolerance. This test is applied via
    # `purrr::map_lgl()` because `dustify()` doubled the length of both vectors.
    matches_sd <- purrr::map_lgl(
      .x = sd,
      .f = function(x) any(dplyr::near(x, sd_rec_rounded, tol = tolerance))
    )

    # TEST 2: If none of the reconstructed SDs matches the reported one, the
    # inputs are GRIMMER-inconsistent.
    test2_results[ii] <- any(matches_sd[!is.na(matches_sd)])
    # if test2_results[ii] is FALSE, no need to run test 3, fill test3_results[ii] with FALSE and continue to the next iteration of the loop
    if (!test2_results[ii]) {
      test3_results[ii] <- FALSE
      next
    }
    

    # Determine if any integer between the lower and upper bounds has the same
    # parity (i.e., the property of being even or odd) as the reconstructed sum:
    matches_parity <- sum_real %% 2 == integers_possible %% 2

    matches_sd_and_parity <- purrr::map_lgl(
      .x = matches_parity,
      .f = function(x) any(x & matches_sd)
    )

    # TEST 3: At least one none of the reconstructed SDs has to match the reported
    # one, and the corresponding reconstructed sums have to match in parity.
    test3_results[ii] <- any(matches_sd_and_parity)

  }
  # convert test1_results, test2_results, test3_results to a dataframe of 3 rows (test1, test2, test3) and length(x_real_possible) columns
  test_results <- rbind(test1_results, test2_results, test3_results)

  # if any of the columns has all TRUE values, then the inputs are GRIMMER-consistent
  pass3 <- apply(test_results, 2, function(x) all(x))
  if (any(pass3)) {
    if (show_reason) {
      return(list(TRUE, "Passed all"))
      } 
      return(TRUE)
  }
  # if we get here, it means that there is no column with all TRUE values, so the inputs are GRIMMER-inconsistent. We can check which test(s) failed.
  # if the best column has test 1 and test 2 as TRUE, but not test 3, then the inputs are GRIMMER-inconsistent because of test 3.
  pass2 <- apply(test_results, 2, function(x) x[1] && x[2] && !x[3])
  if (any(pass2)) {
    if (show_reason) {
      return(list(FALSE, "GRIMMER inconsistent (test 3)"))
    }
    return(FALSE)
  }
  # If the best column has only test 1 as TRUE, but not test 2 (we don't care about test 3 in this case), then the inputs are GRIMMER-inconsistent because of test 2.
  pass1 <- apply(test_results, 2, function(x) x[1] && !x[2])
  if (any(pass1)) {
    if (show_reason) {
      return(list(FALSE, "GRIMMER inconsistent (test 2)"))
    }
    return(FALSE)
  }
  # if we get here, it means that there isn't any column with 3 TRUEs
  # there isn't any column with test 1 and test 2 as TRUE
  # there isn't any column with only test 1 as TRUE
  # all the columns have test 1 as FALSE, we can conclude that the inputs are GRIMMER-inconsistent because of test 1.
  if (show_reason) {
    return(list(FALSE, "GRIMMER inconsistent (test 1)"))
  }
  return(FALSE)
}


#' The GRIMMER test (granularity-related inconsistency of means mapped to error
#' repeats)
#'
#' @description `grimmer()` checks if reported mean and SD values of integer
#'   data are mathematically consistent with the reported sample size and the
#'   number of items that compose the mean value. It works much like [`grim()`].
#'
#'   The function is vectorized, but it is recommended to use [`grimmer_map()`]
#'   for testing multiple cases.
#'
#' @param x String. The reported mean value.
#' @param sd String. The reported standard deviation.
#' @param n Integer. The reported sample size.
#' @param items Integer. The number of items composing the `x` and `sd` values.
#'   Default is `1`, the most common case.
#' @param show_reason Logical. For internal use only. If set to `TRUE`, the
#'   output is a list of length-2 lists which also contain the reasons for
#'   inconsistencies. Don't specify this manually; instead, use `show_reason` in
#'   [`grimmer_map()`]. See there for explanation. Default is `FALSE`.
#'
#' @inheritParams grim
#'
#' @return Logical. `TRUE` if `x`, `sd`, `n`, and `items` are mutually
#'   consistent, `FALSE` if not.

#' @details GRIMMER was originally devised by Anaya (2016). The present
#'   implementation follows Allard's (2018) refined Analytic-GRIMMER algorithm.
#'   It uses a variant of Analytic-GRIMMER first implemented in
#'   \href{https://lukaswallrich.github.io/rsprite2/reference/GRIMMER_test.html}{`rsprite2::GRIMMER_test()`}
#'   that can be applied to multi-item scales.
#'
#'   The scrutiny version embeds GRIMMER in the broader system of consistency
#'   testing, as laid out in
#'   \href{https://lhdjung.github.io/scrutiny/articles/consistency-tests-in-depth.html}{*Consistency
#'   tests in depth*}. The `grimmer()` function
#'   is a vectorized (multiple-case) version of this basic implementation. For
#'   more context and variable name translations, see the top of the R/grimmer.R
#'   source file.

#' @references Allard, A. (2018). Analytic-GRIMMER: a new way of testing the
#'   possibility of standard deviations.
#'   https://aurelienallard.netlify.app/post/anaytic-grimmer-possibility-standard-deviations/
#'
#'   Anaya, J. (2016). The GRIMMER test: A method for testing the validity of
#'   reported measures of variability. *PeerJ Preprints.*
#'   https://peerj.com/preprints/2400v1/

#' @export
#'
#' @examples
#' # A mean of 5.23 is not consistent with an SD of 2.55
#' # and a sample size of 35:
#' grimmer(x = "5.23", sd = "2.55", n = 35)
#'
#' # However, mean and SD are consistent with a
#' # sample size of 31:
#' grimmer(x = "5.23", sd = "2.55", n = 31)
#'
#' # For a scale composed of two items:
#' grimmer(x = "2.74", sd = "0.96", n = 63, items = 2)

grimmer <- Vectorize(grimmer_scalar)
