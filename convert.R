## Sequence file converter from FASTA and NEXUS to sequential PHYLIP
## to accompany codonmodeltest.
## John S. S. Denton

require(stringr, quiet = TRUE)
require(ape, quiet = TRUE)

seq.check <- function(seqfile) 

{
	seq.name <- str_sub(seqfile, 1, unique(as.numeric(str_locate(seqfile, "\\."))) - 1)
	stored <- paste(seq.name, ".phy", sep="")

	if (isTRUE(str_detect(seqfile, pattern = "\\.phy"))) {
		seqfile <- stored
		}

	if (isTRUE(str_detect(seqfile, pattern = "\\.fasta"))) {
		seq.tmp <- read.FASTA(seqfile)
		write.dna(seq.tmp, file = stored, format = "sequential", nbcol= -1, colsep="")
		seqfile <- stored
		}

	if (isTRUE(str_detect(seqfile, pattern = "\\.nex"))) {
		seq.tmp <- read.nexus.data(seqfile)
		write.dna(seq.tmp, file = stored, format = "sequential", nbcol= -1, colsep="")
		seqfile <- stored
		}
}
