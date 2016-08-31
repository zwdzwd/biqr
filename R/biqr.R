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
#' aln <- localaln(seq.fns.3F1R, "chr8", 11664775, 11667135)
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

#' Local align one clone
#'
#' Local align one clone
#'
#' @param seq.fn sequence file name
#' @param refseq reference
#' @return CloneAlignment1
#' @export
clone.aln1 <- function(seq.fn, refseq) {
  
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
  
  structure(list(aln = maxaln,
                 meta = list(strand = strand)),
            class="CloneAlignment1")
}

#' Local align multiple clones
#'
#' Local align multiple clones
#'
#' @param seq.fns a vector of sequence file names
#' @param chrm target chromosome
#' @param beg target begin
#' @param end target end
#' @param min.score minimum alignment score for a query to be considered
#' @param snames a vector of sample names
#' @return CloneAlignment
#' @export
clone.aln <- function(seq.fns, chrm, beg, end, min.score=100, snames=NULL) {
  
  ## load reference
  refseq <- get.reference(chrm, beg, end)
  
  if (missing(snames))
    snames <- basename(seq.fns)
  stopifnot(length(unique(snames)) == length(snames))
  
  ## align each sequence
  alns <- lapply(seq.fns, function(x) clone.aln1(x, refseq))
  alns <- setNames(alns, snames)
  
  ## merge to base matrix
  alns.hi.score <- alns[sapply(alns, function(x) x$aln@score >min.score)]
  m <- rbind(strsplit(as.character(refseq),"")[[1]],
             do.call(rbind, lapply(alns.hi.score, function(x) {
               strsplit(as.character(aligned(x$aln)),"")[[1]]})))
  rownames(m) <- c('reference', names(alns.hi.score))
  
  ## trim m off the plasmid from two sides
  plrange <- which(apply(m, 2, function(x) sum(x=="-") / dim(m)[1]) < 0.7)
  m <- m[,min(plrange):max(plrange)]
  
  structure(list(m = m, alns = alns), class="CloneAlignment")
}

#' Format a base matrix alignment
#' 
#' Format a base matrix alignment
#' 
#' @param m base matrix for alignment
#' @param col.width column width
#' @return a string
#' @export
formatmat <- function(m, col.width=80) {
  nc <- ncol(m)
  nr <- nrow(m)
  snames <- rownames(m)
  s <- ""
  for (i in 1:((nc+col.width)/col.width)) {
    for (j in 1:nr) {
      s <- paste0(s, sprintf("%20s\t%s\n", substr(snames[j],1,20), paste0(m[j,((i-1)*col.width+1):min(i*col.width,nc)],collapse="")))
    }
    s <- paste0(s, "\n")
  }
  cat(s)
  invisible(s)
}

#' Subset cytosine from CloneAlignment
#' 
#' Subset cytosine from CloneAlignment
#' 
#' @param m base matrix
#' @param context "hcg" / "gch" / "cg"
#' @return a matrix of cytosine alignment
#' @export
subsetCytosine <- function(m, context="hcg") {
  if (context == "hcg") {
    locs <- gregexpr("[ACT]CG", paste0(m[1,],collapse=""))[[1]]+1
  } else if (context == "gch") {
    locs <- gregexpr("GC[ACT]", paste0(m[1,],collapse=""))[[1]]+1
  } else if (context == "cg") {
    locs <- gregexpr("CG", paste0(m[1,],collapse=""))[[1]]+1
  }

  m <- m[,locs]
  colnames(m) = locs
  m
}

#' Order clones by cytosine retention
#' 
#' Order clones by cytosine retention
#' 
#' This function only operates on the base matrix, but it requires the reference be intact (before subsetting).
#' 
#' @param ca CloneAlignment
#' @param context "hcg"/"gch"/"cg"
#' @param subset.c whether to subset to cytosine of the specified context, default is TRUE
#' @return a CloneAlignment (subset.c=FALSE) or a sorted base matrix (subset.c==TRUE)
#' @export
order.cytosine <- function(ca, context="hcg", subset.c=TRUE) {
  cyto.wref <- subsetCytosine(ca$m, context)
  cyto <- cyto.wref[2:nrow(cyto.wref),]
  cyto.order <- c(1, order(apply(cyto, 1, function(x) sum(x=='C')))+1)
  if (subset.c)
    cyto.wref[cyto.order,]
  else {
    ca$m <- ca$m[cyto.order,]
    ca
  }
}

#' Plot cytosine
#'
#' Plot cytosine
#'
#' @param m CloneAlignment object or base matrix of cytosine alignment
#' @param context either "hcg", "gch" or "cg"
#' @param add superimpose plot without names
#' @param draw whether to draw grobs
#' @param mar.left left margin
#' @param mar.bottom bottom margin
#' @param mar.top top margin
#' @param mar.right right margin
#' @param bar.top height of top bar
#' @param alpha transparency
#' @param radius plot radius
#' @param loc.beg start of plotting range (default to left-most cytosine)
#' @param loc.end end of plotting range (default to right-most cytosine)
#' @param text.sname.pad padding to sample name plot
#' @param show.sname whether to show sample names
#' @import grid
#' @export
plotcytosine <- function(m, context="hcg", add=FALSE, draw=TRUE,
                         mar.bottom = 0.05, mar.left = 0.05, mar.top = 0.1, mar.right = 0.05,
                         bar.top=0.05, alpha=0.9, radius=3, loc.beg=NULL, loc.end=NULL, text.sname.pad=0.02, show.sname=FALSE) {
  
  if (class(m) == "CloneAlignment")
    m <- order.cytosine(m, context, subset.c=TRUE)
  
  if (is.null(colnames(m)))
    stop("Matrix input must contain location information.")
  
  ## methylation color
  if (context == "hcg") {
    methcol <- 'black'
  } else if (context == "gch") {
    methcol <- '#40E0D0'
  } else if (context == "cg") {
    methcol <- 'black'
  }
  ##cgd.locs <- gregexpr("CG[AGT]", paste0(basemat[1,],collapse=""))[[1]]+1
  ##dgc.locs <- gregexpr("[AGT]GC", paste0(basemat[1,],collapse=""))[[1]]+1

  ## compute coordinates
  ms <- m[2:nrow(m),]
  nr <- nrow(ms)
  nc <- ncol(ms)
  locs <- as.numeric(colnames(ms))
  if (is.null(loc.beg))
    loc.beg <- min(locs)
  if (is.null(loc.end))
    loc.end <- max(locs)
  loc.range <- loc.end - loc.beg
  xlocs <- mar.left + (locs-loc.beg) * (1-mar.left-mar.right) / loc.range
  ylocs <- mar.bottom + (1:nr) * (1-mar.top-mar.bottom-bar.top) / (nr+1)
  xs <- matrix(rep(xlocs, nr), nrow=nr, byrow = T)
  ys <- matrix(rep(ylocs, nc), ncol=nc)
  
  ## plot
  p <- list()
  if (!add) {
    ## title
    p[[length(p)+1]] <- grid.text(sprintf("%d bp", loc.range), 0.5, 1-mar.top*0.7, draw=FALSE, 
                                  just = c('center','bottom'), gp=gpar(fontsize=30))
    
    ## ruler
    s150 <- 150 / (loc.end-loc.beg) * (1-mar.left-mar.right)
    p[[length(p)+1]] <- grid.rect(0.99, 0.99, s150, mar.top*0.2, draw=FALSE, 
                                  just=c('right', 'top'), gp=gpar(color='black', lwd=3))
    p[[length(p)+1]] <- grid.text("150 bp", 0.98-s150, 0.99-mar.top*0.1,
                                  draw=FALSE, just=c('right','center'))
    
    ## top bar
    p[[length(p)+1]] <- grid.lines(x=c(0,1),y=1-mar.top, draw=FALSE, gp=gpar(lwd=3))
    
    ## sample names
    if (show.sname)
      p[[length(p)+1]] <- grid.text(rownames(ms), x=mar.left - text.sname.pad, y=ylocs, draw=FALSE, just=c("right","center"))
  }

  ## circles
  p[[length(p)+1]] <- grid.segments(xlocs, 1-mar.top, xlocs, 1-mar.top-bar.top, draw=FALSE, gp=gpar(lwd=3))
  p[[length(p)+1]] <- grid.circle(xs, ys, r=radius/10/nr, draw=FALSE, 
                                  gp = gpar(
                                    alpha = alpha, lwd = 3,
                                    lty = ifelse(ms=='-', 0, 1),
                                    fill = ifelse(ms=='C', methcol, ifelse(ms=='T', 'white', 0))))
  g <- do.call(gList, p)
  if (draw) {
    grid.newpage()
    grid.draw(g)
  }
  invisible(g)
}

