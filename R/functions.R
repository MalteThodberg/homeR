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

#' Parse known motifs from Homer
#'
#' Parse known motifs form a Homer analysis from tempdir(). Should not normally be called by the user.
#'
#' @return A dataframe of known motif results.
#' @author Malte Thodberg
#' @details Parses from the "knownResults.txt" file of a Homer analysis.
#' @seealso \code{\link{call_homer}} \code{\link{tempdir}}
#' @import magrittr dplyr
#' @export
parse_known <- function(){
	# Read simple table
	known_motifs <- read.table(file=file.path(tempdir(), "knownResults.txt"), sep="\t", header=T, comment.char="")

	# Reformat
	colnames(known_motifs) <- c("Name", "Consensus", "Pval", "Log-pval", "Qval", "Tcount", "T", "Bcount", "B")

	known_motifs %>%
		mutate(T=gsub(pattern="%", replacement="", x=T) %>% as.numeric %>% divide_by(100),
					 B=gsub(pattern="%", replacement="", x=B) %>% as.numeric %>% divide_by(100)) %>%
		select(-Tcount, -Bcount)

	# Return
	known_motifs
}

#' Parse de-novo motifs from Homer
#'
#' Parse de-novo motifs from a Homer analysis from tempdir(). Should not normally be called by the user.
#'
#' @return A list containing a dataframe of de-novo motif results and a list of corresponding PWMs.
#' @author Malte Thodberg
#' @details Parses from the "homerResults/" folder and "homerMotifs.all.motifs" file of a Homer analysis and merges the results.
#' @seealso \code{\link{call_homer}} \code{\link{tempdir}}
#' @import magrittr stringr tidyr dplyr
#' @export
parse_homer <- function(){
	### Folder motifs

	# Files
	motif_fnames <- list.files(file.path(tempdir(), "homerResults"), full.names=TRUE)
	motif_fnames <- motif_fnames[grepl(pattern=".motif", x=motif_fnames, fixed=TRUE)]

	# Seperate PWMS
	homer_PWMs <- lapply(motif_fnames, read.table, skip=1, col.names=c("A", "C", "G", "T"))

	# Reformat info
	folder_motifs <- lapply(motif_fnames, read.table, nrows=1, sep="\t") %>%
		Reduce(rbind, .) %>%
		as.data.frame()
	colnames(folder_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence")

	# Split names
	folder_motifs <- separate(folder_motifs, Name, into=c("Name", "Guess"), sep=",", extra="drop")

	# Rename PWMs
	names(homer_PWMs) <- folder_motifs$Name

	### File motifs

	# Read from file starting with >
	heap <- readLines(file.path(tempdir(), "homerMotifs.all.motifs"))
	heap <- heap[grep(pattern=">", x=heap)]

	# Reformat info
	file_motifs <- heap %>%
		str_split(pattern="\t") %>%
		as.data.frame %>%
		t %>%
		as.data.frame
	colnames(file_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence", "Statistics")
	rownames(file_motifs) <- NULL

	file_motifs <- file_motifs %>%
		mutate(`Log-odds`=as.numeric(`Log-odds`),
					 `Log-pval`=as.numeric(`Log-pval`),
					 Placeholder=as.integer(Placeholder))

	### Merge info

	# Merge frames
	homer_motifs <- left_join(folder_motifs, file_motifs, by=c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence"))

	# Split columns
	homer_motifs <- separate(homer_motifs, Occurence, into=c("T", "B", "P"), sep=",", extra="drop")
	homer_motifs <- separate(homer_motifs, Statistics, into=c("Tpos", "Tstd", "Bpos", "Bstd", "StrandBias", "Multiplicity"), sep=",", extra="drop")

	# Reformat
	homer_motifs <- homer_motifs %>%
		mutate(Consensus=gsub(pattern=">", replacement="", x=Consensus),
					 guess=gsub(pattern="BestGuess:", replacement="", x=Guess),
					 T=str_split(string=T, pattern=":|\\(|\\)|%") %>% sapply(function(x) as.numeric(x[3]) / 100),
					 B=str_split(string=B, pattern=":|\\(|\\)|%") %>% sapply(function(x) as.numeric(x[3]) / 100),
					 P=gsub(pattern="P:", replacement="", x=P) %>% as.numeric,
					 Tpos=gsub(pattern="Tpos:", replacement="", x=Tpos) %>% as.numeric,
					 Tstd=gsub(pattern="Tstd:", replacement="", x=Tstd) %>% as.numeric,
					 Bpos=gsub(pattern="Bpos:", replacement="", x=Bpos) %>% as.numeric,
					 Bstd=gsub(pattern="Bstd:", replacement="", x=Bstd) %>% as.numeric,
					 StrandBias=gsub(pattern="StrandBias:", replacement="", x=StrandBias) %>% as.numeric,
					 Multiplicity=gsub(pattern="Multiplicity:", replacement="", x=Multiplicity) %>% as.numeric) %>%
		mutate(orientation=ifelse(test=grepl(pattern="RV", x=motif_fnames), yes="reverse", no="forward"),
					 length=sapply(homer_PWMs, nrow))

	### Trim

	final_motifs <- !is.na(homer_motifs$Tpos)

	homer_motifs <- subset(homer_motifs, final_motifs, select=-Placeholder)
	homer_PWMs <- homer_PWMs[final_motifs]

	### Return
	list(homer_motifs=homer_motifs, homer_PWMs=homer_PWMs)
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
#' @param mset See findMotifsGenome.pl
#' @param basic See findMotifsGenome.pl
#' @param bits See findMotifsGenome.pl
#' @param nocheck See findMotifsGenome.pl
#' @param mcheck See findMotifsGenome.pl
#' @param noknown See findMotifsGenome.pl
#' @param mknown See findMotifsGenome.pl
#' @param p See findMotifsGenome.pl
#' @return List with output: command line used, and (if selected), known results, de-novo results and de-novo PWMs.
#' @author Malte Thodberg
#' @details Simple R-wrapper for Homer's findMotifsGenome.pl. Saves all temporary files to tempdir().
#' @seealso \code{\link{GR_to_BED}} \code{\link{tempdir}}
#' @export
call_homer <- function(pos_file, genome, # Mandatory
											 mask=NULL, bg=NULL, len=NULL, size=NULL, S=NULL, mis=NULL, norevopp=NULL, nomotif=NULL, rna=NULL, # Basic options
											 mset=NULL, basic=NULL, bits=NULL, nocheck=NULL, mcheck=NULL, noknown=NULL, mknown=NULL, #Known motif options
											 p=NULL){ # Homer options

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
	if(!is.null(mset)){cline <- paste(cline, "-mset", mset)}
	if(!is.null(rna)){cline <- paste(cline, "-rna", rna)}
	if(!is.null(mset)){cline <- paste(cline, "-mset", mset)}
	if(!is.null(mcheck)){cline <- paste(cline, "-mcheck", mcheck)}
	if(!is.null(mknown)){cline <- paste(cline, "-mknown", mknown)}
	if(!is.null(p)){cline <- paste(cline, "-p", p)}

	# Add TRUE/FALSE settings.
	if(!is.null(nomotif)){cline <- paste(cline, "-nomotif")}
	if(!is.null(basic)){cline <- paste(cline, "-basic")}
	if(!is.null(bits)){cline <- paste(cline, "-bits")}
	if(!is.null(nocheck)){cline <- paste(cline, "-nocheck")}
	if(!is.null(noknown)){cline <- paste(cline, "-noknown")}

	# Execute!
	system(cline)

	# Read results back into R
	res <- list(command=cline)

	if(is.null(nomotif)){
		res["known_motifs"] <- parse_known()
	}

	if(is.null(noknown)){
		tmp <- parse_homer()
		res["homer_motifs"] <- tmp[[1]]
		res["homer_PWMs"] <- tmp[[2]]
	}

	# Return cline
	res
}
