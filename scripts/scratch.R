generateseq <- function(length, whicherr, whereerr) {
  s <- rep("0", length)
  s[whereerr] <- c(as.character(1:9), LETTERS)[whicherr]
  paste(s, collapse = "")
}

generateseqs <- function(n, length, p, sub = TRUE, indel = FALSE) {
  nerr <- rbinom(n = n, size = length, prob = p)
  whereerr <- whicherr <- rep(list(integer()), n)
  haserr <- which(nerr > 0)
  whicherr[haserr] <- lapply(nerr[haserr], sample.int, n = 3 * sub + 5 * indel, replace = TRUE)
  whereerr[haserr] <- lapply(nerr[haserr], sample.int, n = 300, replace = FALSE)
  list(
    nerr = nerr,
    whereerr = whereerr,
    whicherr = whicherr,
    seqs = mapply(generateseq, length = length, whicherr = whicherr, whereerr = whereerr)
  )
}

illumina <- generateseqs(10000, 300, 0.005, sub = TRUE, indel = FALSE)
pacbio <- generateseqs(1000, 1500, 0.005, sub = TRUE, indel = TRUE)

sum(illumina$nerr == 0)
sum(pacbio$nerr == 0)
sum(illumina$nerr == 1)
sum(pacbio$nerr == 1)

nerr <- rbinom(n = 1000, size = 1500, prob = 0.005)
where
