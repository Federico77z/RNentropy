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

simes_reference <- function(lpv, gene_status)
{
  pvalues <- rep(NA_real_, nrow(lpv))
  names(pvalues) <- rownames(lpv)
  for(gene in which(gene_status == "TESTED"))
  {
    p <- sort(10 ^ -lpv[gene, is.finite(lpv[gene, ])])
    if(length(p) > 0)
    {
      pvalues[gene] <- min(length(p) * p / seq_along(p))
    }
  }
  pvalues
}

test_that("Simes combination matches hand-computed values", {
  expect_equal(RNentropy:::.RN_simes(0.03), 0.03)
  expect_equal(RNentropy:::.RN_simes(c(0.1, 0.001)), 0.002)
  expect_equal(RNentropy:::.RN_simes(c(0.01, 0.01)), 0.01)
  expect_equal(RNentropy:::.RN_simes(c(0.04, 0.02, 0.03)), 0.04)
  expect_equal(RNentropy:::.RN_simes(c(1, 1, 1)), 1)
})

test_that("gene-level p-values are corrected across TESTED genes", {
  input <- iso_select_fixture()
  results <- RN_iso_select(input)

  expect_equal(results$gpv, -log10(c(gene_1 = 0.002, gene_2 = 0.01,
    gene_3 = NA, gene_4 = 1e-4)))
  expected <- rep(NA_real_, 4)
  names(expected) <- rownames(input$lpv)
  expected[-3] <- -log10(p.adjust(c(0.002, 0.01, 1e-4), method = "BH"))
  expect_equal(results$gpv_adj, expected)
  expect_null(results$lpv_adj)

  expect_identical(rownames(results$selected), c("gene_4", "gene_1", "gene_2"))
  expect_identical(names(results$selected), c("gene_status", "ISO_GPV",
    "CORR_ISO_GPV", colnames(input$lpv)))
  expect_equal(unname(results$selected$CORR_ISO_GPV),
    unname(expected[c("gene_4", "gene_1", "gene_2")]))
  expect_equal(as.matrix(results$selected[colnames(input$lpv)]),
    input$lpv[c("gene_4", "gene_1", "gene_2"), ])
  expect_identical(results$expr, input$expr)
  expect_identical(results$design, input$design)
  expect_identical(results$lpv, input$lpv)
  expect_identical(results$gene_status, input$gene_status)
  expect_identical(results$sample_status, input$sample_status)
  expect_identical(results$res, input$res)
})

test_that("threshold is inclusive and correction method is configurable", {
  input <- iso_select_fixture()
  simes <- simes_reference(input$lpv, input$gene_status)
  family <- !is.na(simes)

  for(method in p.adjust.methods)
  {
    results <- RN_iso_select(input, method = method)
    expect_equal(unname(results$gpv_adj[family]),
      unname(-log10(p.adjust(simes[family], method = method))), info = method)
    expect_true(is.na(results$gpv_adj["gene_3"]), info = method)
  }

  uncorrected <- RN_iso_select(input, gpv_t = 1e-4, method = "none")
  expect_identical(rownames(uncorrected$selected), "gene_4")

  bonferroni <- RN_iso_select(input, gpv_t = 0.006, method = "bonferroni")
  expect_identical(rownames(bonferroni$selected), c("gene_4", "gene_1"))
  expect_equal(unname(bonferroni$gpv_adj["gene_1"]), -log10(0.006))
})

test_that("non-testable values are excluded and empty selections are valid", {
  input <- iso_select_fixture()
  input$lpv["gene_1", 1] <- Inf

  results <- RN_iso_select(input, gpv_t = 1e-10)

  expect_equal(unname(results$gpv["gene_1"]), 1)
  expect_true(is.na(results$gpv["gene_3"]))
  expect_s3_class(results$selected, "data.frame")
  expect_identical(nrow(results$selected), 0L)
  expect_identical(names(results$selected), c("gene_status", "ISO_GPV",
    "CORR_ISO_GPV", colnames(input$lpv)))
})

test_that("TESTED genes without finite sample p-values are excluded", {
  input <- iso_select_fixture()
  input$lpv["gene_1", ] <- NA

  results <- RN_iso_select(input, method = "bonferroni")

  expect_true(is.na(results$gpv["gene_1"]))
  expect_true(is.na(results$gpv_adj["gene_1"]))
  expect_equal(unname(results$gpv_adj[c("gene_2", "gene_4")]),
    -log10(c(0.02, 2e-4)))
  expect_false("gene_1" %in% rownames(results$selected))
})

test_that("empty RN_iso_calc results can be selected", {
  input <- RN_iso_calc(data.frame(gene = character(), sample = numeric()),
    gene.col = "gene")

  results <- RN_iso_select(input)

  expect_length(results$gpv, 0L)
  expect_length(results$gpv_adj, 0L)
  expect_identical(nrow(results$selected), 0L)
  expect_identical(names(results$selected),
    c("gene_status", "ISO_GPV", "CORR_ISO_GPV", "ISO_LPV_sample"))
})

test_that("ties in corrected signal are broken by raw signal, then input order", {
  input <- iso_select_fixture()
  input$gene_status["gene_3"] <- "TESTED"
  input$lpv["gene_1", ] <- c(-log10(0.006), NA)
  input$lpv["gene_2", ] <- c(NA, -log10(0.004))
  input$lpv["gene_3", ] <- c(NA, -log10(0.5))
  input$lpv["gene_4", ] <- c(-log10(0.004), NA)

  results <- RN_iso_select(input, gpv_t = 0.1)

  expect_equal(unname(results$gpv_adj[c("gene_1", "gene_2", "gene_4")]),
    rep(-log10(0.008), 3))
  expect_identical(rownames(results$selected),
    c("gene_2", "gene_4", "gene_1"))
})

test_that("invalid selector inputs are rejected", {
  input <- iso_select_fixture()

  expect_error(RN_iso_select(NULL), "output of RN_iso_calc")
  expect_error(RN_iso_select(input, gpv_t = 0), "greater than 0")
  expect_error(RN_iso_select(input, gpv_t = 1.1), "less than or equal")
  expect_error(RN_iso_select(input, gpv_t = NA_real_), "single p-value")
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

test_that("complete S7 selection agrees with a direct Simes calculation", {
  data("RN_IsoSwitch_Example_S7", package = "RNentropy")
  input <- RN_iso_calc(RN_IsoSwitch_Example_S7, gene.col = "GENE_ID")
  results <- RN_iso_select(input)

  simes <- simes_reference(input$lpv, input$gene_status)
  family <- !is.na(simes)
  adjusted <- p.adjust(simes[family], method = "BH")
  expected <- names(adjusted)[adjusted <= 0.01]

  expect_identical(sum(family), 4058L)
  expect_equal(results$gpv_adj[family], -log10(adjusted))
  expect_identical(nrow(results$selected), 34L)
  expect_setequal(rownames(results$selected), expected)
  expect_true(all(results$selected$gene_status == "TESTED"))
  expect_false(is.unsorted(rev(results$selected$CORR_ISO_GPV)))
})
