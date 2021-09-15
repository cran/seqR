## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE---------------------------------------------------------
library(seqR)

## -----------------------------------------------------------------------------
count_kmers(c("AAA", "ACDA"),
            k=2)

## -----------------------------------------------------------------------------
count_kmers(list(c("A", "A", "A"), c("A", "C", "D", "A")),
            k=2)

## -----------------------------------------------------------------------------
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            k=5)

## -----------------------------------------------------------------------------
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            kmer_gaps=c(2, 0)) # gapped 3-mer template: X__XX

## -----------------------------------------------------------------------------
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            k=5,
            positional=TRUE)

## -----------------------------------------------------------------------------
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            kmer_gaps=c(2, 0), # gapped 3-mer template: X__XX
            positional=TRUE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
get_random_string_vector <- function(seq_num, seq_len, alphabet, seed=123) {
  set.seed(seed)
  sapply(1:seq_num,
         function(n) paste0(sample(alphabet, seq_len, replace=TRUE), collapse=""))
}

r <- count_kmers(get_random_string_vector(300, 1000, LETTERS),
                 k=4)
print(pryr::object_size(r))
print(pryr::object_size(as.matrix(r)))

## -----------------------------------------------------------------------------
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA"),
            k=4)

## -----------------------------------------------------------------------------
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA"),
            k=4,
            kmer_alphabet = c("A", "K"))

## -----------------------------------------------------------------------------
sequences <- get_random_string_vector(seq_num = 200, seq_len = 1000, alphabet = LETTERS)
microbenchmark::microbenchmark(
  multithreaded = count_kmers(sequences, k=5, batch_size = 200),
  singlethreaded = count_kmers(sequences, k=5, batch_size = 1),
  times=11
)

## -----------------------------------------------------------------------------
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA", "ASDDA"),
            k=3,
            batch_size=1,
            verbose = TRUE)

## -----------------------------------------------------------------------------
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_counts=FALSE)

## -----------------------------------------------------------------------------
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_counts=TRUE)

## -----------------------------------------------------------------------------
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_names=FALSE)

## -----------------------------------------------------------------------------
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_names=TRUE)

## -----------------------------------------------------------------------------
count_multimers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
                k_vector=c(1,2,3,4),
                verbose=TRUE)

## -----------------------------------------------------------------------------
mA <- count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
                  k=1)
mB <- count_kmers(c("VVVVVAAVA", "ADSDDD", "AAAAAV"),
                  k=1)

rbind_columnwise(mA, mB)

