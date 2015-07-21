#' GRanges to BED
#'
#' Saves an GRanges-object as a BED-file and returns the file location.
#'
#' @param GR GRanges-object.
#' @return Path to file.
#' @author Malte Thodberg
#' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths. The BED file is save to temp_dir().
#' @import GenomicRanges
#' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' @export
GR_to_BED <- function(GR){
	# Get bed
	bed <- data.frame(seqnames(GR),
									 start(GR)-1,
									 end(GR),
									 1:length(GR),
									 width(GR),
									 strand(GR))

	# Random file name
	bed_fname <- tempfile()

	# Write to file
	write.table(x=bed,
							file=bed_fname,
							quote=FALSE,
							sep="\t",
							row.names=FALSE,
							col.names=FALSE)

	# Return name of file
	bed_fname
}

