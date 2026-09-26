iso_select_fixture <- function()
{
  pv <- matrix(c(
    3, 1,
    2, 2,
    10, 10,
    NA, 4
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("gene_", 1:4),
    c("ISO_PV_sample_1", "ISO_PV_sample_2")))
  gene.status <- c(gene_1 = "TESTED", gene_2 = "TESTED",
    gene_3 = "NOTEST(EXP)", gene_4 = "TESTED")
  sample.status <- matrix("TESTED", nrow = 4, ncol = 2,
    dimnames = dimnames(pv))
  sample.status["gene_3", ] <- NA
  sample.status["gene_4", 1] <- "NO_CURRENT_EXPRESSION"

  list(
    expr = data.frame(gene = rep(rownames(pv), each = 2)),
    design = diag(2),
    pv = pv,
    gene_status = gene.status,
    sample_status = sample.status,
    res = data.frame(gene_status = unname(gene.status), pv,
      row.names = rownames(pv), check.names = FALSE)
  )
}

simes_reference <- function(pv, gene_status)
{
  pvalues <- rep(NA_real_, nrow(pv))
  names(pvalues) <- rownames(pv)
  for(gene in which(gene_status == "TESTED"))
  {
    p <- sort(10 ^ -pv[gene, is.finite(pv[gene, ])])
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

test_that("gene-level Simes p-values are compared with the threshold", {
  input <- iso_select_fixture()
  results <- RN_iso_select(input)

  expect_equal(results$gene_pv, -log10(c(gene_1 = 0.002, gene_2 = 0.01,
    gene_3 = NA, gene_4 = 1e-4)))
  expect_null(results$gpv)
  expect_null(results$gpv_adj)
  expect_null(results$lpv_adj)

  expect_identical(rownames(results$selected), c("gene_4", "gene_1", "gene_2"))
  expect_identical(names(results$selected), c("gene_status", "ISO_GENE_PV",
    colnames(input$pv)))
  expect_equal(unname(results$selected$ISO_GENE_PV),
    unname(results$gene_pv[c("gene_4", "gene_1", "gene_2")]))
  expect_equal(as.matrix(results$selected[colnames(input$pv)]),
    input$pv[c("gene_4", "gene_1", "gene_2"), ])
  expect_identical(results$expr, input$expr)
  expect_identical(results$design, input$design)
  expect_identical(results$pv, input$pv)
  expect_identical(results$gene_status, input$gene_status)
  expect_identical(results$sample_status, input$sample_status)
  expect_identical(results$res, input$res)
})

test_that("threshold is inclusive", {
  input <- iso_select_fixture()

  expect_identical(rownames(RN_iso_select(input, gene_pv = 1e-4)$selected),
    "gene_4")
  expect_identical(rownames(RN_iso_select(input, gene_pv = 0.002)$selected),
    c("gene_4", "gene_1"))
  expect_identical(nrow(RN_iso_select(input, gene_pv = 9e-5)$selected), 0L)
})

test_that("non-testable values are excluded and empty selections are valid", {
  input <- iso_select_fixture()
  input$pv["gene_1", 1] <- Inf

  results <- RN_iso_select(input, gene_pv = 1e-10)

  expect_equal(unname(results$gene_pv["gene_1"]), 1)
  expect_true(is.na(results$gene_pv["gene_3"]))
  expect_s3_class(results$selected, "data.frame")
  expect_identical(nrow(results$selected), 0L)
  expect_identical(names(results$selected), c("gene_status", "ISO_GENE_PV",
    colnames(input$pv)))
})

test_that("TESTED genes without finite sample p-values are excluded", {
  input <- iso_select_fixture()
  input$pv["gene_1", ] <- NA

  results <- RN_iso_select(input, gene_pv = 1)

  expect_true(is.na(results$gene_pv["gene_1"]))
  expect_equal(unname(results$gene_pv[c("gene_2", "gene_4")]), c(2, 4))
  expect_identical(rownames(results$selected), c("gene_4", "gene_2"))
})

test_that("empty RN_iso_calc results can be selected", {
  input <- RN_iso_calc(data.frame(gene = character(), sample = numeric()),
    gene.col = "gene")

  results <- RN_iso_select(input)

  expect_length(results$gene_pv, 0L)
  expect_identical(nrow(results$selected), 0L)
  expect_identical(names(results$selected),
    c("gene_status", "ISO_GENE_PV", "ISO_PV_sample"))
})

test_that("ties in gene-level signal retain input order", {
  input <- iso_select_fixture()
  input$pv["gene_1", ] <- c(3, NA)
  input$pv["gene_2", ] <- c(NA, 5)
  input$pv["gene_4", ] <- c(NA, 3)

  results <- RN_iso_select(input)

  expect_identical(rownames(results$selected),
    c("gene_2", "gene_1", "gene_4"))
})

test_that("invalid selector inputs are rejected", {
  input <- iso_select_fixture()

  expect_error(RN_iso_select(NULL), "output of RN_iso_calc")
  expect_error(RN_iso_select(input, gene_pv = 0), "greater than 0")
  expect_error(RN_iso_select(input, gene_pv = 1.1), "less than or equal")
  expect_error(RN_iso_select(input, gene_pv = NA_real_), "single p-value")

  invalid <- input
  invalid$gene_status <- invalid$gene_status[-1]
  expect_error(RN_iso_select(invalid), "one value for each row")

  invalid <- input
  names(invalid$gene_status)[1] <- "other_gene"
  expect_error(RN_iso_select(invalid), "must match")

  invalid <- input
  invalid$pv["gene_1", 1] <- -1
  expect_error(RN_iso_select(invalid), "non-negative")
})

test_that("complete S7 selection agrees with a direct Simes calculation", {
  data("RN_IsoSwitch_Example_S7", package = "RNentropy")
  input <- RN_iso_calc(RN_IsoSwitch_Example_S7, gene.col = "GENE_ID")
  results <- RN_iso_select(input)

  simes <- simes_reference(input$pv, input$gene_status)
  testable <- !is.na(simes)
  expected <- names(simes)[testable & simes <= 0.01]

  expect_identical(sum(testable), 4058L)
  expect_equal(results$gene_pv[testable], -log10(simes[testable]))
  expect_identical(nrow(results$selected), 111L)
  expect_setequal(rownames(results$selected), expected)
  expect_true(all(results$selected$gene_status == "TESTED"))
  expect_false(is.unsorted(rev(results$selected$ISO_GENE_PV)))
})
