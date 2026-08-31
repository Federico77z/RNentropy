iso_select_fixture <- function()
{
  lpv <- matrix(c(
    3, 1,
    2, 2,
    10, 10,
    NA, 4
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("gene_", 1:4),
    c("ISO_LPV_sample_1", "ISO_LPV_sample_2")))
  gene.status <- c(gene_1 = "TESTED", gene_2 = "TESTED",
    gene_3 = "NOTEST(EXP)", gene_4 = "TESTED")
  sample.status <- matrix("TESTED", nrow = 4, ncol = 2,
    dimnames = dimnames(lpv))
  sample.status["gene_3", ] <- NA
  sample.status["gene_4", 1] <- "NO_CURRENT_EXPRESSION"

  list(
    expr = data.frame(gene = rep(rownames(lpv), each = 2)),
    design = diag(2),
    lpv = lpv,
    gene_status = gene.status,
    sample_status = sample.status,
    res = data.frame(gene_status = unname(gene.status), lpv,
      row.names = rownames(lpv), check.names = FALSE)
  )
}

test_that("iso-switch correction uses all eligible gene-sample tests", {
  input <- iso_select_fixture()
  results <- RN_iso_select(input)

  eligible <- matrix(input$gene_status == "TESTED", nrow = nrow(input$lpv),
    ncol = ncol(input$lpv)) & is.finite(input$lpv)
  expected <- matrix(NA_real_, nrow = nrow(input$lpv), ncol = ncol(input$lpv),
    dimnames = dimnames(input$lpv))
  expected[eligible] <- -log10(p.adjust(10 ^ -input$lpv[eligible], method = "BH"))
  colnames(expected) <- c("CORR_ISO_LPV_sample_1",
    "CORR_ISO_LPV_sample_2")

  expect_equal(results$lpv_adj, expected)
  expect_identical(rownames(results$selected), c("gene_4", "gene_1"))
  expect_identical(names(results$selected), c("gene_status",
    colnames(input$lpv), colnames(expected)))
  expect_identical(results$expr, input$expr)
  expect_identical(results$design, input$design)
  expect_identical(results$lpv, input$lpv)
  expect_identical(results$gene_status, input$gene_status)
  expect_identical(results$sample_status, input$sample_status)
  expect_identical(results$res, input$res)
})

test_that("threshold is inclusive and correction method is configurable", {
  input <- iso_select_fixture()

  for(method in p.adjust.methods)
  {
    results <- RN_iso_select(input, method = method)
    eligible <- is.finite(input$lpv) &
      matrix(input$gene_status == "TESTED", nrow = nrow(input$lpv),
        ncol = ncol(input$lpv))
    expect_equal(results$lpv_adj[eligible],
      -log10(p.adjust(10 ^ -input$lpv[eligible], method = method)),
      info = method)
  }

  uncorrected <- RN_iso_select(input, lpv_t = 0.01, method = "none")
  expect_identical(rownames(uncorrected$selected),
    c("gene_4", "gene_1", "gene_2"))

  bonferroni <- RN_iso_select(input, lpv_t = 0.0005, method = "bonferroni")
  expect_identical(rownames(bonferroni$selected), "gene_4")
  expect_equal(bonferroni$lpv_adj["gene_4", 2], -log10(0.0005))
})

test_that("non-testable values are excluded and empty selections are valid", {
  input <- iso_select_fixture()
  input$lpv["gene_1", 1] <- Inf

  results <- RN_iso_select(input, lpv_t = 1e-10)

  expect_true(is.na(results$lpv_adj["gene_1", 1]))
  expect_true(all(is.na(results$lpv_adj["gene_3", ])))
  expect_s3_class(results$selected, "data.frame")
  expect_identical(nrow(results$selected), 0L)
  expect_identical(names(results$selected), c("gene_status",
    colnames(input$lpv), colnames(results$lpv_adj)))
})

test_that("ties in corrected signal retain input order", {
  input <- iso_select_fixture()
  input$lpv["gene_1", ] <- c(4, 1)
  input$lpv["gene_4", ] <- c(NA, 4)

  results <- RN_iso_select(input, lpv_t = 1, method = "none")

  expect_identical(rownames(results$selected),
    c("gene_1", "gene_4", "gene_2"))
})

test_that("invalid selector inputs are rejected", {
  input <- iso_select_fixture()

  expect_error(RN_iso_select(NULL), "output of RN_iso_calc")
  expect_error(RN_iso_select(input, lpv_t = 0), "greater than 0")
  expect_error(RN_iso_select(input, lpv_t = 1.1), "less than or equal")
  expect_error(RN_iso_select(input, lpv_t = NA_real_), "single p-value")
  expect_error(RN_iso_select(input, method = "unknown"),
    "p.adjust.methods")

  invalid <- input
  invalid$gene_status <- invalid$gene_status[-1]
  expect_error(RN_iso_select(invalid), "one value for each row")

  invalid <- input
  names(invalid$gene_status)[1] <- "other_gene"
  expect_error(RN_iso_select(invalid), "must match")

  invalid <- input
  invalid$lpv["gene_1", 1] <- -1
  expect_error(RN_iso_select(invalid), "non-negative")
})

test_that("complete S7 selection agrees with a direct adjusted-p calculation", {
  data("RN_IsoSwitch_Example_S7", package = "RNentropy")
  input <- RN_iso_calc(RN_IsoSwitch_Example_S7, gene.col = "GENE_ID")
  results <- RN_iso_select(input)

  eligible <- matrix(input$gene_status == "TESTED", nrow = nrow(input$lpv),
    ncol = ncol(input$lpv)) & is.finite(input$lpv)
  adjusted <- matrix(NA_real_, nrow = nrow(input$lpv), ncol = ncol(input$lpv),
    dimnames = dimnames(input$lpv))
  adjusted[eligible] <- p.adjust(10 ^ -input$lpv[eligible], method = "BH")
  expected <- rownames(input$lpv)[apply(adjusted, 1,
    function(x) any(x <= 0.01, na.rm = TRUE))]

  expect_gt(length(expected), 0L)
  expect_setequal(rownames(results$selected), expected)
  expect_true(all(results$selected$gene_status == "TESTED"))
  expect_true(all(apply(results$lpv_adj[expected, , drop = FALSE], 1,
    function(x) any(x >= 2, na.rm = TRUE))))
})
