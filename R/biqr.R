#' @title
#' Plot NOMe-seq data
#' 
#' @description
#' create bubble plot of NOMe-seq data
#' 
#' @details
#' This package implements visualization of NOMeSeq data
#' 
#' @author
#' Wanding Zhou \email{Wanding.Zhou@vai.org}
#'
#' @examples
#' \dontrun{
#' library(biqr)
#' setwd('~/tools/biqr')
#' seq.fns <- paste0('1134071/', list.files('1134071', pattern=".seq"))
#' seq.fns.3F1R <- seq.fns[grepl("3F1R", seq.fns)]
#' seq.fns.2F2R <- seq.fns[grepl("2F2R", seq.fns)]
#' seq.fns.2F2R
#' seq.fns.3F1R
#' 
#' aln <- localaln(seq.fns.2F2R, "chr8", 11664775, 11667135)
#' aln <- localaln(seq.fns.3F1R, chr8, 11664775, 11667135)
#' plotcytosine(aln$basemat[1:20,], 'hcg')
#' plotcytosine(aln$basemat[1:20,], 'gch')
#' }
#' @details
#'
"_PACKAGE"

#' Get reference sequence
#'
#' Get reference sequence from DAS server
#' 
#' @param chrm chromosome
#' @param beg region start
#' @param end region end
#' @return DNAString
#' @importFrom XML xmlToList
#' @importFrom XML xmlParse
#' @import Biostrings
#' @export
get.reference <- function(chrm, beg, end) {
  dna <- xmlToList(xmlParse(paste0('http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=',chrm,':',beg,',',end)))$SEQUENCE$DNA
  DNAString(toupper(gsub("[\n]", "", dna$text)))
}

#' Local align one query
#'
#' Local align one query
#'
#' @param seq.fn sequence file name
#' @param refseq reference
#' @return a list of aligned columns, score and strand
#' @export
localaln1 <- function(seq.fn, refseq) {
  
  ## setup substitution matrix
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  GAmat <- mat
  CTmat <- mat
  GAmat['G','A'] <- 1
  CTmat['T','C'] <- 1
  for (a in c('A','T','C','G')) {
    GAmat[a, 'N'] <- 0
    CTmat[a, 'N'] <- 0
  }
  
  ## read target sequence
  targetseq <- readDNAStringSet(seq.fn)
  targetseq.r <- reverseComplement(targetseq)
  
  ## pairwise alignment
  aln.pw = pairwiseAlignment(targetseq, refseq, type="local", substitutionMatrix=CTmat, gapOpening=3, gapExtension=0.5)
  aln.pc = pairwiseAlignment(targetseq.r, refseq, type="local", substitutionMatrix=GAmat, gapOpening=3, gapExtension=0.5)
  aln.dw = pairwiseAlignment(targetseq.r, refseq, type="local", substitutionMatrix = CTmat, gapOpening=3, gapExtension=0.5)
  aln.dc = pairwiseAlignment(targetseq, refseq, type="local", substitutionMatrix = GAmat, gapOpening=3, gapExtension=0.5)
  
  ## find strand
  strand.index <- which.max(c(aln.pw@score, aln.pc@score, aln.dw@score, aln.dc@score))
  maxaln <- c(aln.pw, aln.pc, aln.dw, aln.dc)[strand.index][[1]]
  if (strand.index == 1 || strand.index == 3) strand = '+';
  if (strand.index == 2 || strand.index == 4) strand = '-';
  ## cat(strand.index, "\t");
  ## cat(maxaln@score, "\n")
  
  list(alncol=strsplit(as.character(aligned(maxaln)),"")[[1]],
       score=maxaln@score,
       strand=strand)
}

#' Local align multiple queries
#'
#' Local align multiple queries
#'
#' @param seq.fns a vector of sequence file names
#' @param chrm target chromosome
#' @param beg target begin
#' @param end target end
#' @param min.score minimum alignment score for a query to be considered
#' @export
localaln <- function(seq.fns, chrm, beg, end, min.score=100) {
  
  ## load reference
  refseq <- get.reference(chrm, beg, end)
  
  ## align each sequence
  seq.alned <- lapply(seq.fns, function(x) localaln1(x, refseq))
  seq.alned <- seq.alned[sapply(seq.alned, function(x) x$score >min.score)]
  
  ## merge to base matrix
  basemat <- rbind(strsplit(as.character(refseq),"")[[1]],
                   do.call(rbind, lapply(seq.alned, function(x) x$alncol)))
  
  ## trim basemat off the plasmid from two sides
  plrange <- which(apply(basemat, 2, function(x) sum(x=="-") / dim(basemat)[1]) < 0.7)
  basemat <- basemat[,min(plrange):max(plrange)]
  
  list(basemat = basemat, 
       strand = sapply(seq.alned, function(x) x$strand))
}

#' Plot cytosine
#'
#' Plot cytosine
#'
#' @param basemat base matrix of the alignment
#' @param target either "hcg", "gch" or "cg"
#' @param margin a vector of margins c(bottom, left, top, right)
#' @param bar.top height of top bar
#' @param alpha transparency
#' @import grid
#' @export
plotcytosine <- function(basemat, target, margin=c(0.05,0.05,0.1,0.05), bar.top=0.05, alpha=0.9) {
  
  mar.bottom <- margin[1]
  mar.left <- margin[2]
  mar.top <- margin[3]
  mar.right <- margin[4]
  
  ## identify cytosine locations
  if (target == "hcg") {
    locs <- gregexpr("[ACT]CG", paste0(basemat[1,],collapse=""))[[1]]+1
    methcol <- 'black'
  } else if (target == "gch") {
    locs <- gregexpr("GC[ACT]", paste0(basemat[1,],collapse=""))[[1]]+1
    methcol <- '#40E0D0'
  } else if (target == "cg") {
    locs <- gregexpr("CG", paste0(basemat[1,],collapse=""))[[1]]+1
    methcol <- 'black'
  }
  ##cgd.locs <- gregexpr("CG[AGT]", paste0(basemat[1,],collapse=""))[[1]]+1
  ##dgc.locs <- gregexpr("[AGT]GC", paste0(basemat[1,],collapse=""))[[1]]+1

  nr <- nrow(basemat)
  nc <- ncol(basemat)
  basemat.c <- basemat[2:nr,locs]

  loc.beg <- min(locs)
  loc.end <- max(locs)
  xlocs <- mar.left + (locs-loc.beg) * (1-mar.left-mar.right) / (loc.end-loc.beg)
  xs <- matrix(rep(xlocs, nr-1), nrow=nr-1, byrow = T)
  ys <- matrix(rep(mar.bottom + (2:nr-1) * (1-mar.top-mar.bottom-bar.top) / nr, 
                   length(locs)), ncol=length(locs))
  
  grid.newpage()
  grid.lines(x=c(0,1),y=1-margin[3], gp=gpar(lwd=5))
  grid.text(sprintf("%d bp", nc), 0.5, 1-mar.top*0.9, just = c('center','bottom'), gp=gpar(fontsize=30))
  grid.segments(xlocs, 1-mar.top, xlocs, 1-mar.top-bar.top, gp=gpar(lwd=5))
  grid.circle(xs, ys, r=0.3/nr, 
              gp = gpar(
                alpha = alpha,
                lty = ifelse(basemat.c=='-', 0, 1),
                fill = ifelse(basemat.c=='C', methcol, ifelse(basemat.c=='T', 0, 'white'))))
}

