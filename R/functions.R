#' GRanges to BED
#'
#' Saves an GRanges-object as a BED-file and returns the file location.
#'
#' @param GR GRanges-object.
#' @param bed_fname Path to where file will be saved. Must end in ".bed". Defaults to tempfile().
#' @param trackline logical: Whether a dummy trackline should be included in the .bed file
#' @return Path to file.
#' @author Malte Thodberg
#' @details Saves a GRanges object as a BED-file. Names are set to increasing integers and scores are set to feature widths (these will be overwritten).
#' @seealso \code{\link{tempdir}} \code{\link{tempfile}}
#' @export
GR_to_BED <- function(GR, bed_fname=NULL, trackline=FALSE){
# 	# Get bed
# 	bed <- data.frame(seqnames(GR),
# 									 start(GR)-1,
# 									 end(GR),
# 									 1:length(GR),
# 									 width(GR),
# 									 strand(GR))

	# Add some info
	GR$name <- 1:length(GR)
	GR$score <- width(GR)

	# Random file name
	if(is.null(bed_fname)){
		bed_fname <- tempfile(pattern="temp_bed_file", fileext=".bed")
	}

	if(trackline == FALSE){
		rtracklayer::export(object=GR, con=bed_fname)
	}else{
		rtracklayer::export(object=GR, con=bed_fname, trackLine=methods::new("BasicTrackLine"))
	}

# 	# Write to file
# 	write.table(x=bed,
# 							file=bed_fname,
# 							quote=FALSE,
# 							sep="\t",
# 							row.names=FALSE,
# 							col.names=FALSE)

	# Return name of file
	bed_fname
}

#' Parse known motifs from Homer
#'
#' Parse known motifs form a Homer analysis from output_dir. Should not normally be called by the user.
#'
#' @param output_dir Path to output dir for Homer analysis. Defaults to tempdir()
#' @return A dataframe of known motif results.
#' @author Malte Thodberg
#' @details Parses from the "knownResults.txt" file of a Homer analysis.
#' @seealso \code{\link{call_homer}} \code{\link{tempdir}}
#' @import magrittr
#' @export
parse_known <- function(output_dir=tempdir()){
	# Read simple table
	known_motifs <- utils::read.table(file=file.path(output_dir, "knownResults.txt"), sep="\t", header=T, comment.char="")

	# Reformat
	colnames(known_motifs) <- c("Name", "Consensus", "Pval", "Log-pval", "Qval", "Tcount", "T", "Bcount", "B")

	known_motifs %>%
		dplyr::mutate(T=gsub(pattern="%", replacement="", x=T) %>% as.numeric %>% divide_by(100),
					 B=gsub(pattern="%", replacement="", x=B) %>% as.numeric %>% divide_by(100)) %>%
		dplyr::select(-Tcount, -Bcount)

	# Return
	known_motifs
}

#' Parse de-novo motifs from Homer
#'
#' Parse de-novo motifs from a Homer analysis from output_dir. Should not normally be called by the user.
#' @param output_dir Path to output dir for Homer analysis. Defaults to tempdir()
#' @return A list containing a dataframe of de-novo motif results and a list of corresponding PWMs.
#' @author Malte Thodberg
#' @details Parses from the "homerResults/" folder and "homerMotifs.all.motifs" file of a Homer analysis and merges the results.
#' @seealso \code{\link{call_homer}} \code{\link{tempdir}}
#' @import magrittr
#' @export
parse_homer <- function(output_dir=tempdir()){
	### Folder motifs

	# Files
	motif_fnames <- list.files(file.path(output_dir, "homerResults"), full.names=TRUE)
	motif_fnames <- motif_fnames[grepl(pattern=".motif", x=motif_fnames, fixed=TRUE)]

	# Seperate PWMS
	homer_PWMs <- lapply(motif_fnames, utils::read.table, skip=1, col.names=c("A", "C", "G", "T"))

	# Reformat info
	folder_motifs <- lapply(motif_fnames, utils::read.table, nrows=1, sep="\t") %>%
		Reduce(rbind, .) %>%
		as.data.frame()
	colnames(folder_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence")

	# Split names
	folder_motifs <- tidyr::separate(folder_motifs, Name, into=c("Name", "Guess"), sep=",", extra="drop")

	# Rename PWMs
	names(homer_PWMs) <- folder_motifs$Name

	### File motifs

	# Read from file starting with >
	heap <- readLines(file.path(output_dir, "homerMotifs.all.motifs"))
	heap <- heap[grep(pattern=">", x=heap)]

	# Reformat info
	file_motifs <- heap %>%
		stringr::str_split(pattern="\t") %>%
		as.data.frame %>%
		t %>%
		as.data.frame
	colnames(file_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence", "Statistics")
	rownames(file_motifs) <- NULL

	file_motifs <- file_motifs %>%
		dplyr::mutate(`Log-odds`=as.numeric(`Log-odds`),
					 `Log-pval`=as.numeric(`Log-pval`),
					 Placeholder=as.integer(Placeholder))

	### Merge info

	# Merge frames
	homer_motifs <- dplyr::left_join(folder_motifs, file_motifs, by=c("Consensus", "Name", "Log-odds", "Log-pval", "Placeholder", "Occurence"))

	# Split columns
	homer_motifs <- tidyr::separate(homer_motifs, Occurence, into=c("T", "B", "P"), sep=",", extra="drop")
	homer_motifs <- tidyr::separate(homer_motifs, Statistics, into=c("Tpos", "Tstd", "Bpos", "Bstd", "StrandBias", "Multiplicity"), sep=",", extra="drop")

	# Reformat
	homer_motifs <- homer_motifs %>%
		dplyr::mutate(Consensus=gsub(pattern=">", replacement="", x=Consensus),
					 Guess=gsub(pattern="BestGuess:", replacement="", x=Guess),
					 T=str_split(string=T, pattern=":|\\(|\\)|%") %>% sapply(function(x) as.numeric(x[3]) / 100),
					 B=str_split(string=B, pattern=":|\\(|\\)|%") %>% sapply(function(x) as.numeric(x[3]) / 100),
					 P=gsub(pattern="P:", replacement="", x=P) %>% as.numeric,
					 Tpos=gsub(pattern="Tpos:", replacement="", x=Tpos) %>% as.numeric,
					 Tstd=gsub(pattern="Tstd:", replacement="", x=Tstd) %>% as.numeric,
					 Bpos=gsub(pattern="Bpos:", replacement="", x=Bpos) %>% as.numeric,
					 Bstd=gsub(pattern="Bstd:", replacement="", x=Bstd) %>% as.numeric,
					 StrandBias=gsub(pattern="StrandBias:", replacement="", x=StrandBias) %>% as.numeric,
					 Multiplicity=gsub(pattern="Multiplicity:", replacement="", x=Multiplicity) %>% as.numeric) %>%
		dplyr::mutate(Orientation=ifelse(test=grepl(pattern="RV", x=motif_fnames), yes="reverse", no="forward"),
					 Length=sapply(homer_PWMs, nrow))

	### Trim

	final_motifs <- !is.na(homer_motifs$Tpos)

	homer_motifs <- subset(homer_motifs, final_motifs, select=-Placeholder)
	homer_PWMs <- homer_PWMs[final_motifs]

	### Return
	list(homer_motifs=homer_motifs, homer_PWMs=homer_PWMs)
}

#' Parse instances of de-novo motifs from Homer
#'
#' Parse de-novo motifs instances from a Homer analysis from output_dir. Should not normally be called by the user.
#' @param output_dir Path to output dir for Homer analysis. Defaults to tempdir()
#' @return A dataframe of the position of motif occurences.
#' @author Malte Thodberg
#' @details Parses from the "homerResults/motifInstances.tab" file of a Homer analysis and merges the results.
#' @seealso \code{\link{call_homer}} \code{\link{find_instances}} \code{\link{tempdir}}
#' @import magrittr
#' @export
parse_instances <- function(output_dir=tempdir()){
	## Read and clean

	# Read data
	i <- utils::read.table(file.path(output_dir, "motifInstances.tab"),
									sep="\t",
									header=TRUE,
									comment.char="",
									quote="")

	# Clean columns
	i <- subset(i, select=c(1:5, 20:ncol(i)))

	info_cols <- colnames(i)[1:7]
	motif_cols <- colnames(i)[8:ncol(i)]

	info_cols[1] <- "PeakID"
	motif_cols <- stringr::str_split(motif_cols, pattern=fixed(".")) %>% sapply(., function(x) x[2])

	colnames(i) <- c(info_cols, motif_cols)

	# Clean rows
	i <- dplyr::arrange(i, PeakID)

	## Make into matrix

	# Matrix of cumbersome strings
	string_mat <- i[,8:ncol(i)]

	# Only get position from strings
	pos_mat <- apply(string_mat, 2, function(x) data.frame(raw=x) %>%
									 	tidyr::separate(col=raw, into=c("pos", "consensus", "strand", "score"), extra="drop", sep=",|\\(|\\)") %>%
									 	.$pos %>%
									 	as.integer)

	# Return
	pos_mat
}

#' Motif Enrichment with Homer
#'
#' Call the findMotifsGenome.pl script from Homer directly from R.
#'
#' @param pos_file <#> (Genomic Ranges object)
#' @param genome <#> (Installed Homer genome) or (path to FASTA)
#' @param output_dir Path to output dir for Homer analysis. Defaults to tempdir()
#' @param mask (mask repeats/lower case sequence, can also add 'r' to genome, i.e. mm9r)
#' @param bg <background position file> (genomic positions to be used as background, default=automatic) removes background positions overlapping with target positions
#' @param chopify (chop up large background regions to the avg size of target regions)
#' @param len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the programto run out of memory - in these cases decrease the number of sequences analyzed (-N), or try analyzing shorter sequence regions (i.e. -size 100)]
#' @param size <#> (fragment size to use for motif finding, default=200) or (i.e. -size -100,50 will get sequences from -100 to +50 relative from center) or given (uses the exact regions you give it)
#' @param S <#> (Number of motifs to optimize, default: 25)
#' @param mis <#> (global optimization: searches for strings with # mismatches, default: 2)
#' @param norevopp (don't search reverse strand for motifs)
#' @param rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)
#' @param mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)
#' @param bits (scale sequence logos by information content, default: doesn't scale)
#' @param mcheck <motif file> (known motifs to check against de novo motifs)
#' @param mknown <motif file> (known motifs to check for enrichment)
#' @param gc (use GC-percentage for sequence content normalization, now the default)
#' @param cpg (use CpG-percentage instead of GC-percentage for sequence content normalization)
#' @param noweight (no CG correction)
#' @param h (use hypergeometric for p-values, binomial is default)
#' @param N <#> (Number of sequences to use for motif finding, default=max(50k, 2x input)
#' @param local <#> (use local background, # of equal size regions around peaks to use i.e. 2)
#' @param redundant <#> (Remove redundant sequences matching greater than # percent, i.e. -redundant 0.5)
#' @param maxN <#> (maximum percentage of N's in sequence to consider for motif finding, default: 0.7)
#' @param maskMotif <motif file1> [motif file 2]... (motifs to mask before motif finding)
#' @param rand (randomize target and background sequences labels)
#' @param ref <peak file> (use file for target and background - first argument is list of peak ids for targets)
#' @param oligo (perform analysis of individual oligo enrichment)
#' @param dumpFasta (Dump fasta files for target and background sequences for use with other programs)
#' @param preparse (force new background files to be created)
#' @param preparsedDir <directory> (location to search for preparsed file and/or place new files)
#' @param keepFiles (keep temporary files)
#' @param fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)
#' @param nlen <#> (length of lower-order oligos to normalize in background, default: -nlen 3)
#' @param nmax <#> (Max normalization iterations, default: 160)
#' @param neutral (weight sequences to neutral frequencies, i.e. 25-percentage, 6.25-percentage, etc.)
#' @param olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)
#' @param p <#> (Number of processors to use, default: 1)
#' @param e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)
#' @param cache <#> (size in MB for statistics cache, default: 500)
#' @param quickMask (skip full masking after finding motifs, similar to original homer)
#' @param minlp <#> (stop looking for motifs when seed logp score gets above #, default: -10)
#' @return List with output: command line used, knowm motifs, Homer motifs (de-novo) and Homer PWMs.
#' @author Malte Thodberg
#' @details
#' Simple R-wrapper for Homer's findMotifsGenome.pl.
#' Instead of flags, it uses R-arguments which are pasted to a Homer command. Flags that modify output format are not implemented:
#' -nomotif, -find, -enhancers, -enhancersOnly, -basic, -nocheck, -noknown, -nofacts, -opt, -peaks, -homer2.
#'
#' Saves all temporary files to output_dir. Note these files are only deleted upon closing the R-session, which can in some cases lead to files from previous runs being reloaded.
#' @seealso \code{\link{GR_to_BED}} \code{\link{tempdir}}
#' @export
call_homer <- function(pos_file, genome, output_dir=tempdir(),# Mandatory
											 mask=NULL, bg=NULL, chopify=NULL, len=NULL, size=NULL, S=NULL, mis=NULL, norevopp=NULL, rna=NULL, # Basic options
											 mset=NULL, bits=NULL, mcheck=NULL, mknown=NULL, #Known motif options
											 gc=NULL, cpg=NULL, noweight=NULL, # Sequence normalization options
											 h=NULL, N=NULL, local=NULL, redundant=NULL, maxN=NULL, maskMotif=NULL, rand=NULL, ref=NULL, oligo=NULL, dumpFasta=NULL, preparse=NULL, preparsedDir=NULL, keepFiles=NULL, fdr=NULL, #Advanced options
											 nlen=NULL, nmax=NULL, neutral=NULL, olen=NULL, p=NULL, e=NULL, cache=NULL, quickMask=NULL, minlp=NULL){ # Homer options

	# GR to temporary file
	pos_file <- GR_to_BED(GR=pos_file)

	if(!is.null(bg)){
		bg <- GR_to_BED(GR=bg)
		}

	# Build basic commond line
	cline <- sprintf("findMotifsGenome.pl %s %s %s -nofacts",
									 pos_file, genome, output_dir)

	if(!is.null(mask)){cline <- paste(cline, "-mask")}
	if(!is.null(bg)){cline <- paste(cline, "-bg", bg)}
	if(!is.null(chopify)){cline <- paste(cline, "-chopify")}
	if(!is.null(len)){cline <- paste(cline, "-len", len)}
	if(!is.null(size)){cline <- paste(cline, "-size", size)}
	if(!is.null(S)){cline <- paste(cline, "-S", S)}
	if(!is.null(mis)){cline <- paste(cline, "-mis")}
	if(!is.null(norevopp)){cline <- paste(cline, "-norevopp")}
	if(!is.null(rna)){cline <- paste(cline, "-rna")}
	if(!is.null(mset)){cline <- paste(cline, "-mset", mset)}
	if(!is.null(bits)){cline <- paste(cline, "-bits")}
	if(!is.null(mcheck)){cline <- paste(cline, "-mcheck", mcheck)}
	if(!is.null(mknown)){cline <- paste(cline, "-mknown", mknown)}
	if(!is.null(gc)){cline <- paste(cline, "-gc")}
	if(!is.null(cpg)){cline <- paste(cline, "-cpg")}
	if(!is.null(noweight)){cline <- paste(cline, "-noweight")}
	if(!is.null(h)){cline <- paste(cline, "-h")}
	if(!is.null(N)){cline <- paste(cline, "-N", N)}
	if(!is.null(local)){cline <- paste(cline, "-local", local)}
	if(!is.null(redundant)){cline <- paste(cline, "-redundant", redundant)}
	if(!is.null(maxN)){cline <- paste(cline, "-maxN", maxN)}
	if(!is.null(maskMotif)){cline <- paste(cline, "-maskMotif", maskMotif)}
	if(!is.null(rand)){cline <- paste(cline, "-rand")}
	if(!is.null(ref)){cline <- paste(cline, "-ref", ref)}
	if(!is.null(oligo)){cline <- paste(cline, "-oligo")}
	if(!is.null(dumpFasta)){cline <- paste(cline, "-dumpFasta")}
	if(!is.null(preparse)){cline <- paste(cline, "-preparse")}
	if(!is.null(preparsedDir)){cline <- paste(cline, "-preparsedDir", preparsedDir)}
	if(!is.null(keepFiles)){cline <- paste(cline, "-keepFiles")}
	if(!is.null(fdr)){cline <- paste(cline, "-fdr", fdr)}
	if(!is.null(nlen)){cline <- paste(cline, "-nlen", nlen)}
	if(!is.null(nmax)){cline <- paste(cline, "-nmax", nmax)}
	if(!is.null(neutral)){cline <- paste(cline, "-neutral")}
	if(!is.null(olen)){cline <- paste(cline, "-olen", olen)}
	if(!is.null(p)){cline <- paste(cline, "-p", p)}
	if(!is.null(e)){cline <- paste(cline, "-e", e)}
	if(!is.null(cache)){cline <- paste(cline, "-cache", cache)}
	if(!is.null(quickMask)){cline <- paste(cline, "-quickMask")}
	if(!is.null(minlp)){cline <- paste(cline, "-minlp", minlp)}

# 	# Add extra settings with values
# 	if(!is.null(bg)){cline <- paste(cline, "-bg", bg)}
# 	if(!is.null(len)){cline <- paste(cline, "-len", len)}
# 	if(!is.null(size)){cline <- paste(cline, "-size", size)}
# 	if(!is.null(S)){cline <- paste(cline, "-S", S)}
# 	if(!is.null(mis)){cline <- paste(cline, "-mis", mis)}
#
# 	if(!is.null(mset)){cline <- paste(cline, "-mset", mset)}
# 	if(!is.null(mcheck)){cline <- paste(cline, "-mcheck", mcheck)}
# 	if(!is.null(mknown)){cline <- paste(cline, "-mknown", mknown)}
#
# 	if(!is.null(N)){cline <- paste(cline, "-N", N)}
# 	if(!is.null(local)){cline <- paste(cline, "-local", local)}
# 	if(!is.null(redundant)){cline <- paste(cline, "-redundant", redundant)}
# 	if(!is.null(maxN)){cline <- paste(cline, "-maxN", maxN)}
# 	if(!is.null(maskMotif)){cline <- paste(cline, "-maskMotif", maskMotif)}
# 	if(!is.null(ref)){cline <- paste(cline, "-ref", ref)}
# 	if(!is.null(preparsedDir)){cline <- paste(cline, "-preparsedDir", preparsedDir)}
# 	if(!is.null(fdr)){cline <- paste(cline, "-fdr", fdr)}
#
# 	if(!is.null(nlen)){cline <- paste(cline, "-nlen", nlen)}
# 	if(!is.null(nlen)){cline <- paste(cline, "-nlen", nlen)}
# 	if(!is.null(p)){cline <- paste(cline, "-p", p)}
#
# 	# Add TRUE/FALSE settings.
# 	if(!is.null(mask)){cline <- paste(cline, "-mask")}
#
# 	if(!is.null(chopify)){cline <- paste(cline, "-chopify")}
# 	if(!is.null(norevopp)){cline <- paste(cline, "-norevopp")}
# 	if(!is.null(rna)){cline <- paste(cline, "-rna")}
#
# 	if(!is.null(bits)){cline <- paste(cline, "-bits")}
#
# 	if(!is.null(gc)){cline <- paste(cline, "-gc")}
# 	if(!is.null(cpg)){cline <- paste(cline, "-cpg")}
#
# 	if(!is.null(h)){cline <- paste(cline, "-h")}
# 	if(!is.null(rand)){cline <- paste(cline, "-rand")}
# 	if(!is.null(oligo)){cline <- paste(cline, "-oligo")}
# 	if(!is.null(dumpFasta)){cline <- paste(cline, "-dumpFasta")}
# 	if(!is.null(preparse)){cline <- paste(cline, "-preparse")}
# 	if(!is.null(keepFiles)){cline <- paste(cline, "-keepFiles")}

	# Execute!
	system(cline)

	# Read results back into R
	res <- list(command=cline)

	res$known_motifs <- parse_known(output_dir=output_dir)

	tmp <- parse_homer(output_dir=output_dir)
	res$homer_motifs <- tmp$homer_motifs
	res$homer_PWMs <- tmp$homer_PWMs

	# Return cline
	res
}

#' Find instances of de-novo Homer motifs
#'
#' Find instances of de-novo Homer motifs in a GRanges using Homer's annotatePeaks.pl script.
#'
#' @param pos_file <#> (Genomic Ranges object)
#' @param genome <#> (Installed Homer genome) or (path to FASTA)
#' @param output_dir Path to output dir for Homer analysis. Defaults to tempdir()
#' @return A tidy dataframe of instances of motifs.
#' @author Malte Thodberg
#' @details Must point to a Homer output directory to find a 'homerMotifs.all.motifs' file. Saves a the raw output file from annotatePeaks.pl as 'motifInstances.tab" in the same folder.
#' @seealso \code{\link{call_homer}} \code{\link{tempdir}}
#' @export
find_instances <- function(pos_file, genome, output_dir=tempdir()){
	# GR to temporary file
	pos_file <- GR_to_BED(GR=pos_file, trackline=TRUE)

	# Look for motifs
	pos_cline <- sprintf("annotatePeaks.pl %s %s -m %s -size given > %s",
											 pos_file,
											 genome,
											 file.path(output_dir, "homerMotifs.all.motifs"),
											 file.path(output_dir, "motifInstances.tab"))

	# Execute!
	system(pos_cline)

	# Read back into R
	parse_instances(output_dir=output_dir)
}
