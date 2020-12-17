testthat::test_that("Test that plot_amplicon_len produces a plot",
                    testthat::expect_is(plot_amplicon_len(blasto_example),
                                        "ggplot")
)

testthat::test_that("Test that display_tree produces a plot",
                    testthat::expect_is(display_tree(blasto_example),
                                        "gtable")
                    )

testthat::test_that("Test that display_wordcloud produces a plot",
                    testthat::expect_is(display_wordcloud(blasto_example),
                                        "ggplot")
)

testthat::test_that("Test that plot_distance produces a plot",
                    testthat::expect_is(plot_distance(blasto_example),
                                        "gtable")
)

testthat::test_that("Test that plot_primer_mismatch produces a plot",
                    testthat::expect_is(
                      plot_primer_mismatch_locs(blasto_example),
                      "gtable")
)
