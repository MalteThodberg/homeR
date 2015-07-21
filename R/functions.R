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

#' Motif Enrichment with Homer
#'
#' Call the findMotifsGenome.pl from Homer from inside R.
#'
#' @param pos_file See findMotifsGenome.pl
#' @param genome See findMotifsGenome.pl
#' @param mask See findMotifsGenome.pl
#' @param bg See findMotifsGenome.pl
#' @param len See findMotifsGenome.pl
#' @param size See findMotifsGenome.pl
#' @param S See findMotifsGenome.pl
#' @param mis See findMotifsGenome.pl
#' @param norevopp See findMotifsGenome.pl
#' @param nomotif See findMotifsGenome.pl
#' @param rna	See findMotifsGenome.pl
#' @return Homer-output as data.frame.
#' @author Malte Thodberg
#' @details Simple R-wrapper for Homer's findMotifsGenome.pl. Saves all temporary files to tempdir().
#' @seealso \code{\link{GR_to_BED}} \code{\link{tempdir}}
#' @export
call_homer <- function(pos_file, genome, # Mandatory
											 mask=NULL, bg=NULL, len=NULL, size=NULL, S=NULL, mis=NULL, norevopp=NULL, nomotif=NULL, rna=NULL){ # Basic options

	# Build basic commond line
	cline <- sprintf("findMotifsGenome.pl %s %s %s -nofacts",
									 pos_file, genome, tempdir())

	# Add extra settings with values
	if(!is.null(mask)){cline <- paste(cline, "-mask", bg)}
	if(!is.null(bg)){cline <- paste(cline, "-bg", bg)}
	if(!is.null(len)){cline <- paste(cline, "-len", len)}
	if(!is.null(size)){cline <- paste(cline, "-size", size)}
	if(!is.null(S)){cline <- paste(cline, "-S", S)}
	if(!is.null(mis)){cline <- paste(cline, "-mis", mis)}
	if(!is.null(norevopp)){cline <- paste(cline, "-norevopp", norevopp)}
	if(!is.null(rna)){cline <- paste(cline, "-rna", rna)}

	# Add TRUE/FALSE settings.
	if(!is.null(nomotif)){cline <- paste(cline, "-nomotif")}

	# Execute!
	system(cline)

	# Read results back into R
	# To bed

	# Return cline
	cline
}
