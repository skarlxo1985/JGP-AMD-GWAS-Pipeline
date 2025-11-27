# ==========================================================
# R/utils_gwas.R
# Helper functions for AMD GWAS Pipeline
# ==========================================================

# Constants
VALID <- c("A","C","G","T")

# 1. Parse Genotype Vector (e.g., "A|G" -> "A/G", "A/." -> NA)
parse_gt_vec <- function(x){
    if(is.null(x)) return(rep(NA_character_, length(x)))
    x <- as.character(x)
    x <- trimws(x)
    x[x %in% c("",".","-","NA","NaN")] <- NA_character_
    x <- gsub("[|\\\\]", "/", x)        # '|' or '\' to '/'
    out <- vapply(strsplit(x, "/", fixed=TRUE), function(p){
        if(length(p)!=2) return(NA_character_)
        a1 <- toupper(trimws(p[1])); a2 <- toupper(trimws(p[2]))
        if(!(a1 %in% VALID) || !(a2 %in% VALID)) return(NA_character_)
        paste(sort(c(a1,a2)), collapse="/")
    }, FUN.VALUE = character(1))
    out
}

# 2. Count Alleles from Normalized Genotype Vector
allele_counts_from_gt <- function(gt_norm){
    g <- gt_norm[!is.na(gt_norm)]
    if(length(g)==0) return(table(factor(character(0), levels=VALID)))
    alle <- unlist(strsplit(g, "/", fixed=TRUE))
    alle <- alle[alle %in% VALID]
    table(factor(alle, levels=VALID))
}

# 3. Calculate MAF safely
safe_maf <- function(gt_norm){
    ctab <- allele_counts_from_gt(gt_norm)
    n <- sum(ctab)
    if(n==0) return(tibble(N_alleles=0, MAC=NA_real_, MAF=NA_real_,
                           MinorAllele=NA_character_, MajorAllele=NA_character_))
    nz <- ctab[ctab>0]
    if(length(nz)==1){
        return(tibble(N_alleles=n, MAC=0, MAF=0,
                      MinorAllele=NA_character_, MajorAllele=names(nz)[1]))
    }
    freqs <- sort(nz/n, decreasing=TRUE)
    major <- names(freqs)[1]
    minor <- names(freqs)[length(freqs)]
    maf   <- min(freqs)
    MAC   <- round(maf*n)
    tibble(N_alleles=n, MAC=MAC, MAF=maf, MinorAllele=minor, MajorAllele=major)
}

# 4. Check if Biallelic
is_biallelic <- function(gt_norm){
    ctab <- allele_counts_from_gt(gt_norm)
    sum(ctab>0)==2
}

# 5. Check Palindromic (A/T or G/C)
palindromic_flag <- function(minor, major){
    s <- sort(unique(c(minor, major)))
    (length(s)==2 && (all(s==c("A","T")) || all(s==c("C","G"))))
}

# 6. Calculate Dosage (0, 1, 2) from Effect Allele
dosage_from_effect <- function(gt_norm, effect){
    if(is.na(effect) || !(effect %in% VALID)) return(rep(NA_real_, length(gt_norm)))
    out <- rep(NA_real_, length(gt_norm))
    ok <- !is.na(gt_norm)
    parts <- strsplit(gt_norm[ok], "/", fixed=TRUE)
    out[ok] <- vapply(parts, function(p) sum(p==effect), FUN.VALUE = numeric(1))
    out
}

# 7. HWE Exact Test (Controls)
hwe_exact_controls <- function(gt_norm){
    g <- gt_norm[!is.na(gt_norm)]
    if(length(g)==0) return(NA_real_)
    alle <- unlist(strsplit(g, "/", fixed=TRUE))
    alle <- alle[alle %in% VALID]
    alle_u <- sort(unique(alle))
    if(length(alle_u)!=2) return(NA_real_) # not biallelic
    a <- alle_u[1]; b <- alle_u[2]
    AA <- sum(g == paste0(a,"/",a))
    AB <- sum(g == paste0(a,"/",b))
    BB <- sum(g == paste0(b,"/",b))
    n  <- AA + AB + BB
    if(n==0) return(NA_real_)
    tryCatch({
        HardyWeinberg::HWExact(c(AA, AB, BB))$pval
    }, error=function(e) NA_real_)
}
