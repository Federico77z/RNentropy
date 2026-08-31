test_that("iso-switch results reproduce the compact C++ S7 benchmark", {
  input <- read.delim(system.file("extdata", "isoswitch_s7.tsv",
    package = "RNentropy"), check.names = FALSE)
  expected <- read.delim(system.file("extdata", "isoswitch_s7_cpp_expected.tsv",
    package = "RNentropy"), check.names = FALSE, colClasses = "character",
    na.strings = c("", "NA"))

  results <- RN_iso_calc(input, gene.col = "GENE_ID")

  expect_identical(unname(results$gene_status), expected$STATUS)
  expect_identical(names(results$gene_status), expected$GENE_ID)
  expect_identical(dim(results$lpv), c(4L, 6L))

  for(gene in expected$GENE_ID[expected$STATUS == "TESTED"])
  {
    cpp <- unlist(expected[expected$GENE_ID == gene, -(1:2)], use.names = FALSE)
    current.zero <- cpp == "-"
    numeric.values <- !is.na(cpp) & !current.zero

    expect_identical(unname(results$sample_status[gene, current.zero]),
      rep("NO_CURRENT_EXPRESSION", sum(current.zero)))
    expect_equal(unname(results$lpv[gene, numeric.values]), as.numeric(cpp[numeric.values]),
      tolerance = 5e-5)
  }
})

test_that("complete S7 package data retain the historical input structure", {
  data("RN_IsoSwitch_Example_S7", package = "RNentropy")

  expect_s3_class(RN_IsoSwitch_Example_S7, "data.frame")
  expect_identical(dim(RN_IsoSwitch_Example_S7), c(81314L, 7L))
  expect_identical(names(RN_IsoSwitch_Example_S7),
    c("GENE_ID", "brain_s7_1", "heart_s7_1", "kidney_s7_1",
      "liver_s7_1", "lung_s7_1", "muscle_s7_1"))
  expect_identical(length(unique(RN_IsoSwitch_Example_S7$GENE_ID)), 28426L)
  expect_identical(anyDuplicated(row.names(RN_IsoSwitch_Example_S7)), 0L)
  expect_false(any(is.na(RN_IsoSwitch_Example_S7)))
  expect_true(all(vapply(RN_IsoSwitch_Example_S7[-1], is.numeric, logical(1))))
  expect_true(all(as.matrix(RN_IsoSwitch_Example_S7[-1]) >= 0))
})

test_that("eligibility uses individual values strictly above min_expr", {
  input <- data.frame(
    gene = c("boundary", "boundary", "eligible", "eligible"),
    sample_1 = c(1, 0, 2, 0),
    sample_2 = c(0, 1, 0, 2)
  )

  results <- RN_iso_calc(input, gene.col = "gene", min_expr = 1)

  expect_identical(unname(results$gene_status), c("NOTEST(EXP)", "TESTED"))
})

test_that("replicates are excluded together and zero states are distinguished", {
  input <- data.frame(
    gene = c("gene", "gene"),
    sample_a1 = c(2, 0), sample_a2 = c(0, 2),
    sample_b = c(0, 0), sample_c = c(0, 0)
  )
  design <- matrix(c(1, 0, 0,
                     1, 0, 0,
                     0, 1, 0,
                     0, 0, 1), nrow = 4, byrow = TRUE)

  results <- RN_iso_calc(input, "gene", design)

  expect_identical(unname(results$gene_status), "TESTED")
  expect_identical(unname(results$sample_status["gene", ]),
    c("NO_REFERENCE_EXPRESSION", "NO_REFERENCE_EXPRESSION",
      "NO_CURRENT_EXPRESSION", "NO_CURRENT_EXPRESSION"))
  expect_true(all(is.na(results$lpv["gene", ])))
})

test_that("large statistics are capped at 300", {
  input <- data.frame(gene = c("gene", "gene"),
    sample_1 = c(1e8, 0), sample_2 = c(0, 1e8))

  results <- RN_iso_calc(input, "gene")

  expect_identical(unname(results$lpv["gene", ]), c(300, 300))
})

test_that("file wrapper accepts names and original numeric column positions", {
  fixture <- system.file("extdata", "isoswitch_s7.tsv", package = "RNentropy")

  named <- RNentropy_iso_switch(fixture, tr.col = "TR_ID", gene.col = "GENE_ID")
  numbered <- RNentropy_iso_switch(fixture, tr.col = 1, gene.col = 2)

  expect_equal(named$lpv, numbered$lpv)
  expect_identical(named$gene_status, numbered$gene_status)
})

test_that("file wrapper resolves skipped and renamed columns", {
  input <- data.frame(
    transcript = c("tr_1", "tr_2"),
    annotation = c("a", "b"),
    gene = c("gene", "gene"),
    sample_1 = c(2, 0), sample_2 = c(0, 2)
  )
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  write.table(input, path, sep = "\t", quote = FALSE, row.names = FALSE)

  results <- RNentropy_iso_switch(path, tr.col = "transcript", gene.col = 3,
    skip.col = 2, col.names = c("gene_id", "first", "second"))

  expect_identical(names(results$expr), c("gene_id", "first", "second"))
  expect_identical(unname(results$gene_status), "TESTED")
})

test_that("isoform rows do not need to be consecutive", {
  interleaved <- data.frame(
    transcript = c("a_1", "b_1", "a_2", "b_2"),
    gene = c("a", "b", "a", "b"),
    sample_1 = c(10, 8, 2, 3),
    sample_2 = c(2, 3, 10, 8)
  )
  grouped <- interleaved[c(1, 3, 2, 4), ]

  direct.interleaved <- RN_iso_calc(interleaved, gene.col = "gene")
  direct.grouped <- RN_iso_calc(grouped, gene.col = "gene")

  expect_identical(direct.interleaved$gene_status,
    direct.grouped$gene_status)
  expect_identical(direct.interleaved$sample_status,
    direct.grouped$sample_status)
  expect_equal(direct.interleaved$lpv, direct.grouped$lpv)

  paths <- tempfile(fileext = c(".interleaved.tsv", ".grouped.tsv"))
  on.exit(unlink(paths))
  write.table(interleaved, paths[1], sep = "\t", quote = FALSE,
    row.names = FALSE)
  write.table(grouped, paths[2], sep = "\t", quote = FALSE,
    row.names = FALSE)

  file.interleaved <- RNentropy_iso_switch(paths[1], tr.col = "transcript",
    gene.col = "gene")
  file.grouped <- RNentropy_iso_switch(paths[2], tr.col = "transcript",
    gene.col = "gene")

  expect_identical(file.interleaved$gene_status, file.grouped$gene_status)
  expect_identical(file.interleaved$sample_status, file.grouped$sample_status)
  expect_equal(file.interleaved$lpv, file.grouped$lpv)
})

test_that("invalid iso-switch inputs are rejected", {
  valid <- data.frame(gene = c("gene", "gene"),
    sample_1 = c(2, 0), sample_2 = c(0, 2))

  missing <- valid
  missing$sample_1[1] <- NA
  negative <- valid
  negative$sample_1[1] <- -1

  expect_error(RN_iso_calc(missing, "gene"), "finite")
  expect_error(RN_iso_calc(negative, "gene"), "negative")
  expect_error(RN_iso_calc(valid, "gene", min_expr = -1), "non-negative")
  expect_error(RN_iso_calc(valid, "gene", pseudocount = 0), "positive")
  expect_error(RN_iso_calc(valid, "unknown"), "does not identify")
})
