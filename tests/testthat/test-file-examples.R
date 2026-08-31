example_design <- function()
{
  sample.names <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  matrix(c(1, 0, 0,
           1, 0, 0,
           0, 1, 0,
           0, 1, 0,
           0, 0, 1,
           0, 0, 1), nrow = 6, byrow = TRUE,
    dimnames = list(sample.names, c("A", "B", "C")))
}

test_that("packaged file example runs the standard RNentropy workflow", {
  fixture <- system.file("extdata", "RNentropy_example.tsv",
    package = "RNentropy")
  sample.names <- rownames(example_design())

  results <- RNentropy(fixture, tr.col = 1, design = example_design(),
    header = FALSE, skip.col = c(2, 3), col.names = sample.names)
  results <- RN_select(results)

  expect_true(nzchar(fixture))
  expect_identical(dim(results$expr), c(13L, 6L))
  expect_setequal(rownames(results$selected),
    c("tr_a", "tr_b", "tr_c", "tr_ab", "tr_ac", "tr_bc",
      "tr_switch_1", "tr_switch_2"))
  expect_identical(sort(unique(unlist(results$selected[c("A", "B", "C")]))),
    c(-1, 0, 1))
})

test_that("packaged file example runs the isoform-switch workflow", {
  fixture <- system.file("extdata", "RNentropy_example.tsv",
    package = "RNentropy")
  sample.names <- rownames(example_design())

  results <- RNentropy_iso_switch(fixture, tr.col = 1, gene.col = 2,
    design = example_design(), header = FALSE, skip.col = 3,
    col.names = c("GENE_ID", sample.names))
  results <- RN_iso_select(results)

  expect_identical(unname(results$gene_status["gene_switch"]), "TESTED")
  expect_identical(unname(results$gene_status["gene_low"]), "NOTEST(EXP)")
  expect_identical(sum(results$gene_status == "NOTEST(ISO)"), 9L)
  expect_identical(rownames(results$selected), "gene_switch")
})
