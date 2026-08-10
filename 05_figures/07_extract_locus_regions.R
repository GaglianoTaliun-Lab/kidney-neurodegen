#!/usr/bin/env Rscript
###############################################################################
# 07_extract_locus_regions.R   (run on NARVAL — where the sumstats live)
#
# Extract +/-WINDOW regional association data (SNP rsID, CHR, BP, P) for the
# coloc.abf loci we want to plot, from the SAME sumstats the coloc used
# (trait_file_map in coloc_abf_common.R). Writes small per-locus tables:
#   plot_data/<label>__PD.tsv   plot_data/<label>__KID.tsv
# which the locuszoom step (08) reads. Consistent with the analysis, and it
# includes eGFRcys (for the BIN3 cross-biomarker panel).
#
# Run it in your 64 GB interactive session (it loads the big files once each).
#
# Usage:
#   Rscript 07_extract_locus_regions.R [sumstats_dir] [out_dir]
###############################################################################
suppressPackageStartupMessages(library(data.table))
gsd <- function(){ a<-commandArgs(FALSE); f<-sub("^--file=","",a[grep("^--file=",a)]); if(length(f)) normalizePath(dirname(f[1])) else getwd() }
source(file.path(gsd(), "coloc_abf_common.R"))   # trait_file_map, .get_chr_pos, .find_col, P_CANDS, LOG10P_CANDS

args <- commandArgs(TRUE)
sumstats_dir <- if (length(args) >= 1) args[1] else "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
out_dir      <- if (length(args) >= 2) args[2] else "plot_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
WINDOW_KB <- 500

# --- loci to plot: the 4 colocalized loci + BIN3 in eGFRcys + SCARB2 control ---
loci <- data.frame(
  label        = c("RAB19","MSRA_PRSS51","BIN3_eGFR","BIN3_eGFRcys","TNK2","SCARB2_distinct"),
  chr          = c(7, 8, 8, 8, 3, 4),
  center_bp    = c(140117809, 10264000, 22501830, 22501830, 195624393, 77108306),
  pd_stratum   = c("PD (meta)","PD (meta)","PD (meta)","PD (meta)","PD (female)","PD (meta)"),
  kidney_trait = c("eGFR (meta)","eGFR (meta)","eGFR (meta)","eGFRcys (meta)","eGFR (female)","eGFR (meta)"),
  stringsAsFactors = FALSE)

ID_CANDS <- c("RSID","rsID","rsid","ID","SNP","snp","variant_id","rs_id","dbSNP")
find_id_col <- function(d){ for(c in ID_CANDS){ j<-which(tolower(names(d))==tolower(c)); if(length(j)){ v<-as.character(d[[names(d)[j[1]]]][seq_len(min(500L,nrow(d)))]); if(mean(grepl("^rs[0-9]+",v),na.rm=TRUE)>0.5) return(names(d)[j[1]]) } }; NA_character_ }
resolve <- function(f) if(grepl("^(/|~)",f)) path.expand(f) else file.path(sumstats_dir,f)
read_full <- function(fpath){
  is_gz<-grepl("\\.gz$",fpath,ignore.case=TRUE); base<-if(is_gz)paste("zcat",shQuote(fpath))else paste("cat",shQuote(fpath))
  con<-if(is_gz)gzfile(fpath)else file(fpath); hdr<-tryCatch(readLines(con,1L),error=function(e)character(0),finally=close(con))
  if(length(hdr)==1L&&grepl("\t",hdr)&&grepl(" ",hdr)) fread(cmd=paste(base,"| tr -s '[:blank:]' ' '"),sep=" ",showProgress=FALSE)
  else if(is_gz) fread(cmd=base,showProgress=FALSE) else fread(fpath,showProgress=FALSE) }

cache <- new.env()
load_trait <- function(trait){
  if(exists(trait, envir=cache)) return(get(trait, envir=cache))
  f<-trait_file_map[[trait]]
  if(is.null(f)||is.na(f)){ message("!! no file map for ",trait); assign(trait,NULL,envir=cache); return(NULL) }
  fp<-resolve(unname(f)); if(!file.exists(fp)){ message("!! missing ",trait,": ",fp); assign(trait,NULL,envir=cache); return(NULL) }
  message("reading ",trait," -> ",basename(fp))
  d<-as.data.frame(read_full(fp)); cp<-.get_chr_pos(d); nm<-names(d)
  pc<-.find_col(nm,P_CANDS); lc<-.find_col(nm,LOG10P_CANDS)
  P<-if(!is.na(pc)) suppressWarnings(as.numeric(d[[pc]])) else if(!is.na(lc)) 10^(-suppressWarnings(as.numeric(d[[lc]]))) else NA_real_
  idc<-find_id_col(d); SNP<-if(!is.na(idc)) as.character(d[[idc]]) else NA_character_
  res<-data.table(SNP=SNP, CHR=cp$chr, BP=cp$pos, P=P); assign(trait,res,envir=cache); res }
region <- function(trait, chr, center){
  d<-load_trait(trait); if(is.null(d)) return(NULL)
  w<-d[CHR==chr & BP>=(center-WINDOW_KB*1000) & BP<=(center+WINDOW_KB*1000) & !is.na(P) & P>0]
  as.data.frame(w[, .(SNP,CHR,BP,P)]) }

for(i in seq_len(nrow(loci))){
  lab<-loci$label[i]; ch<-loci$chr[i]; ce<-loci$center_bp[i]
  pd<-region(loci$pd_stratum[i],ch,ce); kid<-region(loci$kidney_trait[i],ch,ce)
  if(is.null(pd)||is.null(kid)){ message("  skip ",lab); next }
  # if PD has no rsIDs, borrow them from the kidney file by position (keeps LD colouring working in 08)
  if(mean(grepl("^rs",pd$SNP),na.rm=TRUE) < 0.5){ m<-setNames(kid$SNP, kid$BP); pd$SNP<-unname(m[as.character(pd$BP)]) }
  write.table(pd,  file.path(out_dir, paste0(lab,"__PD.tsv")),  sep="\t", row.names=FALSE, quote=FALSE)
  write.table(kid, file.path(out_dir, paste0(lab,"__KID.tsv")), sep="\t", row.names=FALSE, quote=FALSE)
  message(sprintf("  %-16s  PD %5d SNPs | KID %5d SNPs", lab, nrow(pd), nrow(kid)))
}
message("\nDone -> ", normalizePath(out_dir), "  (copy this folder to wherever you run 08)")
