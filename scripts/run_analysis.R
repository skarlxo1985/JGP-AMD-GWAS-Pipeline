# ==========================================================
# AMD Candidate-SNP GWAS Pipeline (Execution Script)
#
# Description:
#   Performs QC, logistic regression (Firth fallback), subgroup analysis,
#   and report generation for AMD candidate SNPs.
#
# Inputs:
#   - AMD_383_with_PAQ8_R.csv (placed in project root)
#   - RSID_383_information.csv (placed in project root)
# ==========================================================

# ---------- 0) Packages ----------
req_pkgs <- c("tidyverse","readr","stringr","HardyWeinberg","logistf",
              "ggplot2","gridExtra","scales","openxlsx","png", "ggrepel")
inst <- req_pkgs[!suppressWarnings(sapply(req_pkgs, require, character.only=TRUE))]
if(length(inst)>0){
    install.packages(inst, repos="https://cloud.r-project.org", quiet=TRUE)
    sapply(inst, require, character.only=TRUE)
}

# *** Load Helper Functions ***
# Ensure the "R" folder is in your working directory
if(file.exists("R/utils_gwas.R")) {
    source("R/utils_gwas.R") 
} else {
    stop("Helper file 'R/utils_gwas.R' not found. Please check folder structure.")
}

# ---------- 1) I/O & paths ----------
input_geno <- "AMD_383_with_PAQ8_R.csv"
input_rsid <- "RSID_383_information.csv"

# Check file existence
if(!file.exists(input_geno)) stop("Input file not found: ", input_geno)
if(!file.exists(input_rsid)) stop("Input file not found: ", input_rsid)

out_dir <- "AMD_GWAS_out"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 2) Load data ----------
message("Loading data...")
geno_raw <- readr::read_csv(input_geno, 
                            show_col_types = FALSE,
                            locale = locale(encoding = "CP949"))

rs_info  <- readr::read_csv(input_rsid, 
                            show_col_types = FALSE,
                            locale = locale(encoding = "CP949"))

# phenotype processing
df <- geno_raw %>%
    mutate(
        Phenotype = case_when(
            `AMD sorting` == 1 ~ 0,
            `AMD sorting` >= 2 & `AMD sorting` <= 10 ~ 1,
            TRUE ~ NA_real_
        ),
        Age = suppressWarnings(as.numeric(연령)), 
        Sex = as.factor(성별),                    
        PAQ8 = as.factor(PAQ8) 
    )

# analysis set
df_a <- df %>% filter(!is.na(Phenotype), !is.na(Age), Age >= 50)

# SNP columns detection
rs_pool <- unique(as.character(rs_info$RSID))
snp_cols <- intersect(colnames(df_a), rs_pool)
if(length(snp_cols)==0){
    snp_cols <- grep("^rs[0-9]+", colnames(df_a), value=TRUE)
}
stopifnot("No SNP columns found" = length(snp_cols)>0)

# normalize genotypes
gt_norm <- as_tibble(lapply(df_a[, snp_cols, drop=FALSE], parse_gt_vec))

# sample / snp missingness
per_sample_missing <- rowMeans(is.na(gt_norm))
per_snp_missing <- colMeans(is.na(gt_norm))

# ---------- 3) Basic cohort summary ----------
message("Generating cohort summary...")
cohort_summary <- tibble(
    N_raw_rows = nrow(geno_raw),
    N_after_filters = nrow(df_a),
    N_controls = sum(df_a$Phenotype==0, na.rm=TRUE),
    N_cases    = sum(df_a$Phenotype==1, na.rm=TRUE),
    Age_mean = mean(df_a$Age, na.rm=TRUE),
    Age_sd   = sd(df_a$Age, na.rm=TRUE),
    Age_min  = min(df_a$Age, na.rm=TRUE),
    Age_max  = max(df_a$Age, na.rm=TRUE),
    SNP_columns_in_RS_list = sum(colnames(geno_raw) %in% rs_pool),
    SNP_columns_used = length(snp_cols)
)
readr::write_csv(cohort_summary, file.path(out_dir,"cohort_summary.csv"))

if("Sex" %in% colnames(df_a)){
    sex_dist <- df_a %>% count(Sex, name="Count")
    readr::write_csv(sex_dist, file.path(out_dir,"sex_distribution.csv"))
}

if("PAQ8" %in% colnames(df_a)){
    paq8_dist <- df_a %>% count(PAQ8, name="Count")
    readr::write_csv(paq8_dist, file.path(out_dir,"PAQ8_distribution.csv"))
}

missingness_summary <- tibble(
    Per_sample_missing_rate_mean   = mean(per_sample_missing),
    Per_sample_missing_rate_median = median(per_sample_missing),
    Per_snp_missing_rate_mean      = mean(per_snp_missing),
    Per_snp_missing_rate_median    = median(per_snp_missing),
    SNPs_missing_gt_10pct          = sum(per_snp_missing>0.10)
)
readr::write_csv(missingness_summary, file.path(out_dir,"missingness_summary.csv"))

# ---------- 4) Per-SNP QC ----------
message("Running QC...")
idx_ctrl <- which(df_a$Phenotype==0)
idx_case <- which(df_a$Phenotype==1)

qc_rows <- lapply(snp_cols, function(snp){
    g_all <- gt_norm[[snp]]
    call_rate <- 1 - mean(is.na(g_all))
    bi <- is_biallelic(g_all)
    ov <- safe_maf(g_all) %>% rename(MAF_overall=MAF, Minor_overall=MinorAllele, Major_overall=MajorAllele,
                                     N_alleles_overall=N_alleles, MAC_overall=MAC)
    g_ctrl <- g_all[idx_ctrl]; g_case <- g_all[idx_case]
    ct <- safe_maf(g_ctrl) %>% rename(MAF_controls=MAF, Minor_controls=MinorAllele, Major_controls=MajorAllele,
                                      N_alleles_controls=N_alleles, MAC_controls=MAC)
    cs <- safe_maf(g_case) %>% rename(MAF_cases=MAF, Minor_cases=MinorAllele, Major_cases=MajorAllele,
                                      N_alleles_cases=N_alleles, MAC_cases=MAC)
    palin <- palindromic_flag(ct$Minor_controls, ct$Major_controls)
    hwe_p <- tryCatch(hwe_exact_controls(g_ctrl), error=function(e) NA_real_)
    tibble(
        SNP=snp, CallRate=call_rate, Biallelic=bi, Palindromic=palin,
        MAF_overall=ov$MAF_overall, Minor_overall=ov$Minor_overall, Major_overall=ov$Major_overall,
        MAF_controls=ct$MAF_controls, Minor_controls=ct$Minor_controls, Major_controls=ct$Major_controls,
        MAF_cases=cs$MAF_cases, Minor_cases=cs$Minor_cases, Major_cases=cs$Major_cases,
        HWE_p_controls=hwe_p,
        N_alleles_counted_overall=ov$N_alleles_overall, N_alleles_counted_controls=ct$N_alleles_controls, N_alleles_counted_cases=cs$N_alleles_cases,
        MissingRate=1-call_rate
    )
})
qc_df <- bind_rows(qc_rows)
readr::write_csv(qc_df, file.path(out_dir,"per_snp_qc_summary.csv"))

qc_summary <- tibble(
    SNPs_total = nrow(qc_df),
    SNPs_biallelic = sum(qc_df$Biallelic, na.rm=TRUE),
    SNPs_monomorphic_or_nonbiallelic = sum(!qc_df$Biallelic, na.rm=TRUE),
    SNPs_missing_gt_10pct = sum(qc_df$MissingRate>0.10, na.rm=TRUE),
    SNPs_palin_flag = sum(qc_df$Palindromic, na.rm=TRUE),
    Mean_MAF_overall = mean(qc_df$MAF_overall, na.rm=TRUE)
)
readr::write_csv(qc_summary, file.path(out_dir,"SNP_QC_summary.csv"))

# QC Filtering
CALLRATE_MIN <- 0.98
MAF_MIN      <- 0.01      # controls
HWE_MIN      <- 1e-6

qc_df <- qc_df %>%
    mutate(
        PASS_CallRate = CallRate >= CALLRATE_MIN,
        PASS_Biallelic = Biallelic %in% TRUE,
        PASS_MAF = replace_na(MAF_controls, 0) >= MAF_MIN,
        PASS_HWE = replace_na(HWE_p_controls, 0) >= HWE_MIN,
        PASS_ALL = PASS_CallRate & PASS_Biallelic & PASS_MAF & PASS_HWE
    )
filtered_qc <- qc_df %>% filter(PASS_ALL)
filtered_snps <- filtered_qc$SNP

filter_summary <- tibble(
    Total_SNPs=nrow(qc_df),
    Pass_CallRate=sum(qc_df$PASS_CallRate),
    Pass_Biallelic=sum(qc_df$PASS_Biallelic),
    Pass_MAF=sum(qc_df$PASS_MAF),
    Pass_HWE=sum(qc_df$PASS_HWE),
    Pass_All_Criteria=sum(qc_df$PASS_ALL),
    Using_MAF_from="MAF_controls"
)
readr::write_csv(filter_summary, file.path(out_dir,"qc_filter_summary.csv"))
readr::write_lines(filtered_snps, file.path(out_dir,"filtered_snps_list.txt"))
readr::write_csv(filtered_qc, file.path(out_dir,"filtered_snps_qc_table.csv"))

# ---------- 6b) Final analysis cohort (post-QC) age/sex/PAQ8 comparison ----------
SAMPLE_CALLRATE_MIN <- 0.98

if(length(filtered_snps) == 0){
    warning("No SNPs passed QC; skipping final cohort age/sex comparison.")
    final_age_summary <- tibble()
    final_sex_table <- tibble()
    final_paq8_table <- tibble()
    final_tests <- tibble()
} else {
    gt_norm_filtered <- gt_norm[, filtered_snps, drop = FALSE]
    sample_callrate  <- 1 - rowMeans(is.na(gt_norm_filtered))
    
    final_keep <- which(
        sample_callrate >= SAMPLE_CALLRATE_MIN &
            !is.na(df_a$Phenotype) & !is.na(df_a$Age) & !is.na(df_a$Sex) &
            !is.na(df_a$PAQ8) # Exclude PAQ8 missing values
    )
    final_df <- df_a[final_keep, , drop=FALSE]
    
    final_df <- final_df %>%
        mutate(Group = ifelse(Phenotype==0, "Control", "Case"))
    
    final_age_summary <- final_df %>%
        group_by(Group) %>%
        summarize(
            N        = n(),
            Age_mean = mean(Age, na.rm=TRUE),
            Age_sd   = sd(Age, na.rm=TRUE),
            Age_median = median(Age, na.rm=TRUE),
            Age_min  = min(Age, na.rm=TRUE),
            Age_max  = max(Age, na.rm=TRUE),
            .groups="drop"
        )
    
    age_test <- tryCatch(t.test(Age ~ Group, data=final_df), error=function(e) NULL)
    p_age <- if(!is.null(age_test)) age_test$p.value else NA_real_
    
    sx_tab <- table(final_df$Sex, final_df$Group, useNA="no")
    sex_test_method <- NA_character_
    p_sex <- NA_real_
    if(all(dim(sx_tab) >= c(2,2))){
        expected <- suppressWarnings(chisq.test(sx_tab)$expected)
        if(any(expected < 5, na.rm=TRUE)){
            sex_test_method <- "Fisher's exact test"
            p_sex <- tryCatch(fisher.test(sx_tab)$p.value, error=function(e) NA_real_)
        } else {
            sex_test_method <- "Chi-squared test"
            p_sex <- tryCatch(chisq.test(sx_tab, correct=FALSE)$p.value, error=function(e) NA_real_)
        }
    }
    
    final_sex_table <- as_tibble(sx_tab, .name_repair="minimal")
    colnames(final_sex_table) <- c("Sex","Group","Count")
    final_sex_table <- final_sex_table %>%
        tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0) %>%
        mutate(RowTotal = rowSums(across(-Sex)))
    
    paq_tab <- table(final_df$PAQ8, final_df$Group, useNA="no")
    paq_test_method <- NA_character_
    p_paq <- NA_real_
    if(all(dim(paq_tab) >= c(2,2))){
        expected <- suppressWarnings(chisq.test(paq_tab)$expected)
        if(any(expected < 5, na.rm=TRUE)){
            paq_test_method <- "Fisher's exact test"
            p_paq <- tryCatch(fisher.test(paq_tab)$p.value, error=function(e) NA_real_)
        } else {
            paq_test_method <- "Chi-squared test"
            p_paq <- tryCatch(chisq.test(paq_tab, correct=FALSE)$p.value, error=function(e) NA_real_)
        }
    }
    
    final_paq8_table <- as_tibble(paq_tab, .name_repair="minimal")
    colnames(final_paq8_table) <- c("PAQ8","Group","Count")
    final_paq8_table <- final_paq8_table %>%
        tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0) %>%
        mutate(RowTotal = rowSums(across(-PAQ8)))
    
    final_tests <- tibble(
        Test = c("Age: Welch t-test", 
                 paste0("Sex: ", sex_test_method),
                 paste0("PAQ8: ", paq_test_method)), 
        P_value = c(p_age, p_sex, p_paq) 
    )
    
    readr::write_csv(final_age_summary, file.path(out_dir, "final_age_summary_by_group.csv"))
    readr::write_csv(final_sex_table,   file.path(out_dir, "final_sex_distribution_by_group.csv"))
    readr::write_csv(final_paq8_table,  file.path(out_dir, "final_PAQ8_distribution_by_group.csv")) 
    readr::write_csv(final_tests,       file.path(out_dir, "final_age_sex_PAQ8_pvalues.csv")) 
    
    # ==========================================================
    # 6c) Additional Request: Compare subgroups (Dry/Wet, Drusen/Pachy) among QC-passed subjects
    # ==========================================================
    
    sub_comp_dir <- file.path(out_dir, "subgroup_comparisons_age_sex")
    if(!dir.exists(sub_comp_dir)) dir.create(sub_comp_dir, recursive = TRUE)
    
    message("Running requested subgroup comparisons (Dry/Wet, Drusen/Pachy)...")
    
    # --- 1. Dry AMD (2,3,4,6,7) vs Wet AMD (5,8,9) ---
    dry_wet_df <- final_df %>%
        filter(`AMD sorting` %in% c(2,3,4,6,7, 5,8,9)) %>%
        mutate(Group = factor(case_when(
            `AMD sorting` %in% c(2,3,4,6,7) ~ "Dry_AMD",
            `AMD sorting` %in% c(5,8,9)     ~ "Wet_AMD"
        ), levels = c("Dry_AMD", "Wet_AMD"))) %>%
        filter(!is.na(Group))
    
    if(nrow(dry_wet_df) > 10 && length(unique(dry_wet_df$Group)) == 2 &&
       all(table(dry_wet_df$Group) >= 5)){ 
        
        dw_age_summary <- dry_wet_df %>%
            group_by(Group) %>%
            summarize(
                N        = n(),
                Age_mean = mean(Age, na.rm=TRUE),
                Age_sd   = sd(Age, na.rm=TRUE),
                .groups="drop"
            )
        
        dw_age_test <- tryCatch(t.test(Age ~ Group, data=dry_wet_df), error=function(e) NULL)
        p_age_dw <- if(!is.null(dw_age_test)) dw_age_test$p.value else NA_real_
        
        dw_sx_tab <- table(dry_wet_df$Sex, dry_wet_df$Group, useNA="no")
        dw_sex_method <- NA_character_
        p_sex_dw <- NA_real_
        
        if(all(dim(dw_sx_tab) >= c(2,2))){
            expected <- suppressWarnings(chisq.test(dw_sx_tab)$expected)
            if(any(expected < 5, na.rm=TRUE)){
                dw_sex_method <- "Fisher's exact test"
                p_sex_dw <- tryCatch(fisher.test(dw_sx_tab)$p.value, error=function(e) NA_real_)
            } else {
                dw_sex_method <- "Chi-squared test"
                p_sex_dw <- tryCatch(chisq.test(dw_sx_tab, correct=FALSE)$p.value, error=function(e) NA_real_)
            }
        }
        
        dw_sex_table <- as_tibble(dw_sx_tab, .name_repair="minimal")
        colnames(dw_sex_table) <- c("Sex","Group","Count")
        dw_sex_table <- dw_sex_table %>%
            tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0)
        
        dw_tests <- tibble(
            Test = c("Age: Welch t-test", paste0("Sex: ", dw_sex_method)),
            P_value = c(p_age_dw, p_sex_dw)
        )
        
        readr::write_csv(dw_age_summary, file.path(sub_comp_dir, "dry_vs_wet_age_summary.csv"))
        readr::write_csv(dw_sex_table,  file.path(sub_comp_dir, "dry_vs_wet_sex_distribution.csv"))
        readr::write_csv(dw_tests,      file.path(sub_comp_dir, "dry_vs_wet_pvalues.csv"))
        
    } else {
        warning("Skipping Dry vs Wet: insufficient samples for comparison.")
    }
    
    # --- 2. Drusen-driven (2,3,4,5) vs Pachychoroid-driven (6,7,8,9) ---
    dru_pac_df <- final_df %>%
        filter(`AMD sorting` %in% c(2,3,4,5, 6,7,8,9)) %>%
        mutate(Group = factor(case_when(
            `AMD sorting` %in% c(2,3,4,5) ~ "Drusen_driven",
            `AMD sorting` %in% c(6,7,8,9) ~ "Pachychoroid_driven"
        ), levels = c("Drusen_driven", "Pachychoroid_driven"))) %>%
        filter(!is.na(Group))
    
    if(nrow(dru_pac_df) > 10 && length(unique(dru_pac_df$Group)) == 2 &&
       all(table(dru_pac_df$Group) >= 5)){
        
        dp_age_summary <- dru_pac_df %>%
            group_by(Group) %>%
            summarize(
                N        = n(),
                Age_mean = mean(Age, na.rm=TRUE),
                Age_sd   = sd(Age, na.rm=TRUE),
                .groups="drop"
            )
        
        dp_age_test <- tryCatch(t.test(Age ~ Group, data=dru_pac_df), error=function(e) NULL)
        p_age_dp <- if(!is.null(dp_age_test)) dp_age_test$p.value else NA_real_
        
        dp_sx_tab <- table(dru_pac_df$Sex, dru_pac_df$Group, useNA="no")
        dp_sex_method <- NA_character_
        p_sex_dp <- NA_real_
        
        if(all(dim(dp_sx_tab) >= c(2,2))){
            expected <- suppressWarnings(chisq.test(dp_sx_tab)$expected)
            if(any(expected < 5, na.rm=TRUE)){
                dp_sex_method <- "Fisher's exact test"
                p_sex_dp <- tryCatch(fisher.test(dp_sx_tab)$p.value, error=function(e) NA_real_)
            } else {
                dp_sex_method <- "Chi-squared test"
                p_sex_dp <- tryCatch(chisq.test(dp_sx_tab, correct=FALSE)$p.value, error=function(e) NA_real_)
            }
        }
        
        dp_sex_table <- as_tibble(dp_sx_tab, .name_repair="minimal")
        colnames(dp_sex_table) <- c("Sex","Group","Count")
        dp_sex_table <- dp_sex_table %>%
            tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0)
        
        dp_tests <- tibble(
            Test = c("Age: Welch t-test", paste0("Sex: ", dp_sex_method)),
            P_value = c(p_age_dp, p_sex_dp)
        )
        
        readr::write_csv(dp_age_summary, file.path(sub_comp_dir, "drusen_vs_pachy_age_summary.csv"))
        readr::write_csv(dp_sex_table,  file.path(sub_comp_dir, "drusen_vs_pachy_sex_distribution.csv"))
        readr::write_csv(dp_tests,      file.path(sub_comp_dir, "drusen_vs_pachy_pvalues.csv"))
        
    } else {
        warning("Skipping Drusen vs Pachychoroid: insufficient samples for comparison.")
    }
}

# ---------- 7) Association Analysis ----------
message("Running association analysis...")
effect_lookup <- qc_df %>% select(SNP, Minor_controls) %>% deframe()

# Define association function (Closure to capture df_a, gt_norm)
run_assoc_one <- function(snp){
    g <- gt_norm[[snp]]
    eff <- effect_lookup[[snp]]
    dosage <- dosage_from_effect(g, eff)
    
    sub <- df_a %>%
        mutate(dosage=dosage) %>%
        filter(!is.na(Phenotype), !is.na(Age), !is.na(dosage), !is.na(Sex), !is.na(PAQ8)) 
    
    N_used <- nrow(sub)
    if(N_used < 20 || (length(unique(sub$dosage)) < 2)){
        return(tibble(SNP=snp, EffectAllele=eff, Beta=NA_real_, SE=NA_real_, OR_per_effect=NA_real_,
                      P=NA_real_, N_used=N_used, Method="NA", ErrorMessage="Insufficient data"))
    }
    
    fit <- tryCatch(glm(Phenotype ~ dosage + Age + Sex + PAQ8, data=sub, family=binomial()),
                    error=function(e) NULL)
    
    if(!is.null(fit)){
        co <- summary(fit)$coefficients
        if("dosage" %in% rownames(co)){
            beta <- co["dosage","Estimate"]; se <- co["dosage","Std. Error"]; p <- co["dosage","Pr(>|z|)"]
            return(tibble(SNP=snp, EffectAllele=eff, Beta=beta, SE=se, OR_per_effect=exp(beta),
                          P=p, N_used=N_used, Method="GLM", ErrorMessage=NA_character_))
        }
    }
    
    fitf <- tryCatch(logistf::logistf(Phenotype ~ dosage + Age + Sex + PAQ8, data=sub), error=function(e) NULL)
    if(!is.null(fitf) && "dosage" %in% rownames(fitf$coefficients)){
        beta <- fitf$coefficients["dosage","coef"]; se <- fitf$coefficients["dosage","se"]; p <- fitf$prob["dosage","prob"]
        return(tibble(SNP=snp, EffectAllele=eff, Beta=beta, SE=se, OR_per_effect=exp(beta),
                      P=p, N_used=N_used, Method="Firth", ErrorMessage=NA_character_))
    }
    tibble(SNP=snp, EffectAllele=eff, Beta=NA_real_, SE=NA_real_, OR_per_effect=NA_real_,
           P=NA_real_, N_used=N_used, Method="Failed", ErrorMessage="Fit error")
}

assoc_list <- lapply(filtered_snps, run_assoc_one)
assoc <- bind_rows(assoc_list)

m_eff <- sum(!is.na(assoc$P))
assoc <- assoc %>%
    mutate(FDR_BH = p.adjust(P, method="BH"),
           Bonferroni = ifelse(is.na(P), NA_real_, pmin(P * m_eff, 1.0)))

# Merge Results
keep_qc <- qc_df %>% select(SNP, CallRate, Biallelic, Palindromic,
                            MAF_overall, Minor_overall, Major_overall,
                            MAF_controls, Minor_controls, Major_controls,
                            MAF_cases, Minor_cases, Major_cases,
                            HWE_p_controls)
assoc_out <- assoc %>%
    left_join(keep_qc, by="SNP") %>%
    left_join(rs_info, by=c("SNP"="RSID")) %>%
    arrange(P)

readr::write_csv(assoc_out, file.path(out_dir,"assoc_results_filtered.csv"))
openxlsx::write.xlsx(assoc_out, file.path(out_dir,"assoc_results_filtered.xlsx"), overwrite=TRUE)

# ---------- 8) Plots ----------
message("Generating plots...")

# Volcano
assoc_out <- assoc_out %>%
    mutate(P_safe = ifelse(is.na(P) | P<=.Machine$double.xmin, .Machine$double.xmin, P),
           log10P = -log10(P_safe),
           DeltaMAF_abs = abs(MAF_cases - MAF_controls))

p_volcano <- ggplot(assoc_out, aes(x=Beta, y=log10P)) +
    geom_point(size=1.8, alpha=.8) +
    labs(title="Volcano plot", x="Beta (adj. Age+Sex+PAQ8)", y="-log10(P)") +
    theme_bw(base_size=12)
ggsave(file.path(out_dir,"volcano_beta_log10p.png"), p_volcano, width=8, height=5, dpi=150)

# Scatter
maf_for_x <- ifelse(is.na(assoc_out$MAF_overall), assoc_out$MAF_controls, assoc_out$MAF_overall)
p_maf <- ggplot(assoc_out, aes(x=maf_for_x, y=DeltaMAF_abs)) +
    geom_point(size=1.8, alpha=.8) +
    labs(title="MAF vs |DeltaMAF|", x="MAF (overall; fallback controls)", y="|DeltaMAF| (cases - controls)") +
    theme_bw(base_size=12)
ggsave(file.path(out_dir,"scatter_maf_vs_delta.png"), p_maf, width=8, height=5, dpi=150)

# Manhattan-like (Simple version)
assoc_m <- assoc_out %>% filter(!is.na(CHROM), !is.na(START))
chrom_order <- function(x){
    xx <- as.character(x); xx <- gsub("^chr", "", xx, ignore.case=TRUE)
    key <- suppressWarnings(as.integer(xx))
    ifelse(is.na(key), Inf, key)
}
chr_levels <- assoc_m %>% distinct(CHROM) %>% arrange(chrom_order(CHROM)) %>% pull(CHROM)
assoc_m <- assoc_m %>% mutate(CHROM=factor(CHROM, levels=chr_levels)) %>% arrange(CHROM, START)

offset <- 0
ticks <- c(); ticklabs <- c()
cum_pos <- numeric(nrow(assoc_m))
prev_chr <- NA
for(i in seq_len(nrow(assoc_m))){
    chr <- assoc_m$CHROM[i]
    if(!identical(chr, prev_chr)){
        if(length(cum_pos)>0){
        }
        chri <- which(levels(assoc_m$CHROM)==chr)
        if(i>1) offset <- max(cum_pos, na.rm=TRUE) + 1e6
        prev_chr <- chr
    }
    cum_pos[i] <- assoc_m$START[i] + offset
}
assoc_m$cum_pos <- cum_pos

tick_df <- assoc_m %>% group_by(CHROM) %>% summarize(mid=mean(range(cum_pos))) %>% ungroup()
bonf_alpha <- 0.05
bonf_thr   <- ifelse(m_eff > 0, bonf_alpha / m_eff, NA_real_)
bonf_y     <- ifelse(is.na(bonf_thr), NA_real_, -log10(bonf_thr))

p_manhattan <- ggplot(assoc_m, aes(x=cum_pos, y=log10P)) +
    geom_point(size=1.3, alpha=.85) +
    scale_x_continuous(breaks=tick_df$mid, labels=as.character(tick_df$CHROM)) +
    labs(title="Manhattan-like plot (candidate SNPs)", x="Chromosome", y="-log10(P)") +
    theme_bw(base_size=12)

if(!is.na(bonf_y) && is.finite(bonf_y)){
    p_manhattan <- p_manhattan + geom_hline(yintercept = bonf_y, linetype = 2)
}

ggsave(file.path(out_dir,"manhattan_candidate_snps.png"), p_manhattan, width=10, height=4.5, dpi=150)

# QQ Plot
qq_df <- assoc_out %>% filter(!is.na(P)) %>% arrange(P) %>%
    mutate(obs = -log10(pmax(P, .Machine$double.xmin)),
           exp = -log10(ppoints(n())))

lambda_gc <- {
    chisq_obs <- qchisq(1 - pmin(assoc_out$P, 1 - 1e-16), df = 1)
    median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
}

p_qq <- ggplot(qq_df, aes(x = exp, y = obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
    geom_point(size=1.6, alpha=.8) +
    labs(title = sprintf("QQ plot (lambdaGC = %.3f)", lambda_gc),
         x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_bw(base_size=12)
ggsave(file.path(out_dir, "qq_plot.png"), p_qq, width=6.5, height=6.5, dpi=150)


# ---------- 8c) Locus Zoom Plots ----------
locus_zoom_files <- c() 

if(!is.na(bonf_y) && is.finite(bonf_y)){
    padding_bp <- 250000 
    sig_snps_data <- assoc_m %>% filter(log10P > bonf_y)
    
    if(nrow(sig_snps_data) == 0) {
        message("No SNPs passed the Bonferroni threshold. Skipping Locus Zoom Plot generation.")
    } else {
        sig_chroms <- unique(sig_snps_data$CHROM)
        message("Generating Locus Zoom Plots for significant chromosomes...")
        
        for (chrom in sig_chroms) {
            chrom_data <- assoc_m %>% filter(CHROM == chrom)
            region_bounds <- sig_snps_data %>%
                filter(CHROM == chrom) %>%
                summarise(min_pos = min(START), max_pos = max(START))
            
            zoom_start <- max(0, region_bounds$min_pos - padding_bp)
            zoom_end   <- region_bounds$max_pos + padding_bp
            
            plot_data <- chrom_data %>%
                filter(START >= zoom_start, START <= zoom_end) %>%
                mutate(is_significant = (log10P > bonf_y))
            
            p_zoom <- ggplot(plot_data, aes(x = START, y = log10P)) +
                geom_point(data = . %>% filter(!is_significant), color = "grey60", alpha = 0.7) +
                geom_point(data = . %>% filter(is_significant), color = "red", size = 2) +
                geom_hline(yintercept = bonf_y, linetype = "dashed", color = "red") +
                ggrepel::geom_text_repel(
                    data = . %>% filter(is_significant),
                    aes(label = SNP), 
                    size = 3.5, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0
                ) +
                scale_x_continuous(labels = scales::comma) +
                labs(
                    title = paste("Locus Zoom Plot - Chromosome", chrom),
                    subtitle = paste0("Region: ", scales::comma(zoom_start), " - ", scales::comma(zoom_end), " bp"),
                    x = paste("Position on Chromosome", chrom, "(bp)"),
                    y = "-log10(P)"
                ) +
                theme_bw(base_size = 12) +
                theme(plot.title = element_text(face = "bold"), legend.position = "none")
            
            plot_filename <- file.path(out_dir, paste0("locus_zoom_chr", chrom, ".png"))
            ggsave(plot_filename, p_zoom, width = 10, height = 6, dpi = 150)
            locus_zoom_files <- c(locus_zoom_files, plot_filename)
        }
    }
}

# ---------- 9) Top signals ----------
top_fdr <- assoc_out %>% filter(!is.na(FDR_BH), FDR_BH < 0.05) %>%
    arrange(P) %>%
    select(SNP, CHROM, START, Beta, OR_per_effect, P, FDR_BH, Bonferroni,
           MAF_controls, MAF_cases, DeltaMAF_abs, N_used, Palindromic)
top_bonf <- assoc_out %>% filter(!is.na(Bonferroni), Bonferroni < 0.05) %>%
    arrange(P) %>%
    select(SNP, CHROM, START, Beta, OR_per_effect, P, FDR_BH, Bonferroni,
           MAF_controls, MAF_cases, DeltaMAF_abs, N_used, Palindromic)
readr::write_csv(top_fdr, file.path(out_dir,"top_signals_FDR_lt_0.05.csv"))
readr::write_csv(top_bonf, file.path(out_dir,"top_signals_Bonferroni_lt_0.05.csv"))

# ---------- 9b) LD among Bonferroni-significant SNPs ----------
sig_snps <- assoc_out %>%
    filter(!is.na(Bonferroni), Bonferroni < 0.05) %>%
    pull(SNP) %>% unique()

if(length(sig_snps) >= 2){
    # Reuse filtering objects if available
    if(!exists("gt_norm_filtered")){
        gt_norm_filtered <- gt_norm[, filtered_snps, drop = FALSE]
        sample_callrate  <- 1 - rowMeans(is.na(gt_norm_filtered))
    }
    if(!exists("final_keep")){
        SAMPLE_CALLRATE_MIN <- 0.98
        final_keep <- which(
            sample_callrate >= SAMPLE_CALLRATE_MIN &
                !is.na(df_a$Phenotype) & !is.na(df_a$Age) & !is.na(df_a$Sex) &
                !is.na(df_a$PAQ8) 
        )
    }
    sig_snps <- intersect(sig_snps, colnames(gt_norm))
    sig_snps <- intersect(sig_snps, filtered_snps)
    
    if(length(sig_snps) >= 2){
        eff_allele_lookup <- qc_df %>% select(SNP, Minor_controls) %>% deframe()
        get_dos <- function(snp){
            dosage_from_effect(gt_norm[[snp]], eff_allele_lookup[[snp]])
        }
        dos_list <- lapply(sig_snps, get_dos)
        dos_mat_all <- do.call(cbind, dos_list)
        colnames(dos_mat_all) <- sig_snps
        
        dos_mat <- dos_mat_all[final_keep, , drop=FALSE]
        
        R <- suppressWarnings(cor(dos_mat, use = "pairwise.complete.obs"))
        R2 <- R^2
        
        annot <- assoc_out %>%
            select(SNP, CHROM, START) %>%
            distinct()
        
        ld_long <- as_tibble(R2, rownames = "SNP_A") %>%
            tidyr::pivot_longer(-SNP_A, names_to = "SNP_B", values_to = "r2") %>%
            filter(SNP_A < SNP_B) %>% 
            left_join(annot, by = c("SNP_A" = "SNP")) %>% rename(CHR_A = CHROM, POS_A = START) %>%
            left_join(annot, by = c("SNP_B" = "SNP")) %>% rename(CHR_B = CHROM, POS_B = START) %>%
            mutate(
                same_chr = !is.na(CHR_A) & !is.na(CHR_B) & CHR_A == CHR_B,
                dist_bp  = ifelse(same_chr, abs(POS_A - POS_B), NA_real_)
            ) %>%
            arrange(desc(r2))
        
        readr::write_csv(ld_long, file.path(out_dir, "LD_pairs_Bonferroni_sig_r2.csv"))
        readr::write_csv(as_tibble(R2, rownames = "SNP"), file.path(out_dir, "LD_matrix_Bonferroni_sig_r2.csv"))
        
        order_df <- annot %>%
            filter(SNP %in% sig_snps) %>%
            mutate(
                CHROM_order = {
                    xx <- as.character(CHROM); xx <- gsub("^chr","",xx,ignore.case=TRUE)
                    suppressWarnings(as.integer(xx))
                }
            ) %>%
            arrange(CHROM_order, START)
        
        ord <- order_df$SNP
        ord <- ord[ord %in% colnames(R2)]
        if(length(ord) >= 2){
            R2_ord <- R2[ord, ord, drop=FALSE]
            hm_df <- as_tibble(R2_ord, rownames = "SNP_A") %>%
                tidyr::pivot_longer(-SNP_A, names_to = "SNP_B", values_to = "r2")
            
            p_ld <- ggplot(hm_df, aes(x = SNP_A, y = SNP_B, fill = r2)) +
                geom_tile() +
                scale_fill_gradient(limits=c(0,1), oob=scales::squish) +
                labs(title = "LD heatmap (r2) among Bonferroni-significant SNPs",
                     x = "", y = "", fill = expression(r^2)) +
                theme_bw(base_size=10) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
                      panel.grid = element_blank())
            
            ggsave(file.path(out_dir, "LD_heatmap_Bonferroni_sig_r2.png"),
                   p_ld, width = max(6, min(12, 0.18*length(ord))), height = max(6, min(12, 0.18*length(ord))), dpi=150)
        }
    } else {
        message("Bonferroni-significant SNPs fewer than 2 after intersection; LD skipped.")
    }
} else {
    message("No Bonferroni-significant SNPs; LD skipped.")
}

# ---------- 10) Palindromic strand checklist ----------
assoc_out <- assoc_out %>%
    mutate(MAF_for_flag = ifelse(is.na(MAF_overall), MAF_controls, MAF_overall),
           Palin_highMAF_flag = Palindromic & !is.na(MAF_for_flag) & MAF_for_flag >= 0.4 & MAF_for_flag <= 0.6)
palin_check <- assoc_out %>%
    filter(Palindromic %in% TRUE) %>%
    transmute(SNP, CHROM, START, Palindromic, MAF_for_flag, MAF_controls, MAF_cases,
              DeltaMAF_abs,
              Strand_check_recommended = TRUE,
              Reason = ifelse(Palin_highMAF_flag,
                              "Palindromic with high MAF (~0.4–0.6) — high strand-flip risk",
                              "Palindromic (A/T or C/G); verify against reference allele/AF"))
readr::write_csv(palin_check, file.path(out_dir,"palindromic_strand_checklist.csv"))

# ---------- 11) Proximity clustering (250 kb) ----------
DIST_BP <- 250000
all_with_clusters <- assoc_out %>%
    arrange(CHROM, START) %>%
    group_by(CHROM) %>%
    mutate(ClusterID_local = {
        cid <- integer(n()); cid[1] <- 1
        for(i in 2:n()){
            if(is.na(START[i-1]) || is.na(START[i]) || (START[i] - START[i-1]) > DIST_BP){
                cid[i] <- cid[i-1] + 1
            } else {
                cid[i] <- cid[i-1]
            }
        }
        cid
    }) %>%
    ungroup() %>%
    mutate(ClusterID = paste0("chr", CHROM, "_c", ClusterID_local))
cluster_leads <- all_with_clusters %>%
    arrange(P) %>%
    group_by(ClusterID) %>%
    slice_min(order_by=P, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    select(ClusterID, SNP, CHROM, START, Beta, OR_per_effect, P, FDR_BH, Bonferroni,
           MAF_controls, MAF_cases, DeltaMAF_abs, Palindromic) %>%
    mutate(Literature_Notes="")
readr::write_csv(all_with_clusters, file.path(out_dir,"all_snps_with_clusters_250kb.csv"))
readr::write_csv(cluster_leads, file.path(out_dir,"independent_signal_clusters_250kb.csv"))

lit_template <- cluster_leads %>%
    select(ClusterID, SNP, CHROM, START, P, FDR_BH, Bonferroni) %>%
    mutate(PMID_or_DOI="", Study="", Population="", Reported_EffectAllele="", Notes="")
readr::write_csv(lit_template, file.path(out_dir,"literature_crossref_template.csv"))

# ---------- 11b) QQ plot using LD-pruned set (cluster leads only) ----------
assoc_ldpruned <- assoc_out %>%
    semi_join(cluster_leads %>% select(SNP), by = "SNP") %>%
    filter(!is.na(P)) %>%
    arrange(P)

if(nrow(assoc_ldpruned) >= 5){  # Ensure at least a few points
    qq_ld <- assoc_ldpruned %>%
        mutate(
            obs = -log10(pmax(P, .Machine$double.xmin)),
            exp = -log10(ppoints(n()))
        )
    
    lambda_gc_ld <- {
        chisq_obs <- qchisq(1 - pmin(assoc_ldpruned$P, 1 - 1e-16), df = 1)
        median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df = 1)
    }
    
    p_qq_ld <- ggplot(qq_ld, aes(exp, obs)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
        geom_point(size = 1.6, alpha = .8) +
        labs(title = sprintf("QQ plot (LD-pruned, lambdaGC = %.3f)", lambda_gc_ld),
             x = "Expected -log10(P)", y = "Observed -log10(P)") +
        theme_bw(base_size = 12)
    
    ggsave(file.path(out_dir, "qq_plot_LDpruned_cluster_leads.png"),
           p_qq_ld, width = 6.5, height = 6.5, dpi = 150)
} else {
    message("Too few cluster-lead SNPs for LD-pruned QQ; skipping.")
}

# ==========================================================
# 11c) Subgroup analyses by AMD sorting
# ==========================================================

if(!exists("gt_norm_filtered") || !exists("final_keep")){
    if(length(filtered_snps) > 0){
        gt_norm_filtered <- gt_norm[, filtered_snps, drop = FALSE]
        sample_callrate  <- 1 - rowMeans(is.na(gt_norm_filtered))
        SAMPLE_CALLRATE_MIN <- get0("SAMPLE_CALLRATE_MIN", ifnotfound = 0.98)
        
        final_keep <- which(
            sample_callrate >= SAMPLE_CALLRATE_MIN &
                !is.na(df_a$Phenotype) & !is.na(df_a$Age) & !is.na(df_a$Sex) &
                !is.na(df_a$PAQ8) 
        )
    } else {
        final_keep <- integer(0)
    }
}

run_subgroup <- function(label, case_lo, case_hi){
    message("=== Subgroup: ", label, " (cases ", case_lo, "-", case_hi, ") ===")
    sub_dir <- file.path(out_dir, "subgroups", label)
    if(!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
    
    idx_ctrl <- which(df_a$`AMD sorting` == 1)
    idx_case <- which(df_a$`AMD sorting` >= case_lo & df_a$`AMD sorting` <= case_hi)
    keep_idx <- intersect(final_keep, union(idx_ctrl, idx_case))
    
    sub_pheno <- rep(NA_real_, nrow(df_a))
    sub_pheno[idx_ctrl] <- 0
    sub_pheno[idx_case] <- 1
    sub_pheno <- sub_pheno[keep_idx]
    
    if(length(keep_idx) < 20 || sum(sub_pheno==1, na.rm=TRUE) < 10 || sum(sub_pheno==0, na.rm=TRUE) < 10){
        warning("Subgroup '", label, "': insufficient samples after filters; skipping.")
        return(invisible(NULL))
    }
    
    sub_df <- df_a[keep_idx, , drop=FALSE] %>%
        mutate(Group = ifelse(sub_pheno==0, "Control(1)", paste0(label, " case(", case_lo, "-", case_hi, ")")))
    sub_age_summary <- sub_df %>%
        group_by(Group) %>%
        summarize(
            N = n(), Age_mean = mean(Age, na.rm=TRUE), Age_sd = sd(Age, na.rm=TRUE),
            Age_median = median(Age, na.rm=TRUE), Age_min = min(Age, na.rm=TRUE), Age_max = max(Age, na.rm=TRUE),
            .groups="drop"
        )
    readr::write_csv(sub_age_summary, file.path(sub_dir, "age_summary.csv"))
    
    sx_tab <- table(sub_df$Sex, sub_df$Group, useNA="no")
    sub_sex <- as_tibble(sx_tab, .name_repair="minimal")
    colnames(sub_sex) <- c("Sex","Group","Count")
    sub_sex <- sub_sex %>% tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0)
    readr::write_csv(sub_sex, file.path(sub_dir, "sex_distribution.csv"))
    
    paq_tab <- table(sub_df$PAQ8, sub_df$Group, useNA="no")
    sub_paq8 <- as_tibble(paq_tab, .name_repair="minimal")
    colnames(sub_paq8) <- c("PAQ8","Group","Count")
    sub_paq8 <- sub_paq8 %>% tidyr::pivot_wider(names_from = Group, values_from = Count, values_fill = 0)
    readr::write_csv(sub_paq8, file.path(sub_dir, "PAQ8_distribution.csv"))
    
    
    run_assoc_one_sub <- function(snp, keep_idx, pheno01){
        eff <- effect_lookup[[snp]]
        dosage_all <- dosage_from_effect(gt_norm[[snp]], eff)
        dosage <- dosage_all[keep_idx]
        Age    <- df_a$Age[keep_idx]
        Sex    <- df_a$Sex[keep_idx]
        PAQ8   <- df_a$PAQ8[keep_idx] 
        
        ok <- !is.na(pheno01) & !is.na(Age) & !is.na(Sex) & !is.na(dosage) & !is.na(PAQ8)
        
        subN <- sum(ok)
        if(subN < 20 || length(unique(dosage[ok])) < 2 || length(unique(pheno01[ok])) < 2){
            return(tibble(SNP=snp, EffectAllele=eff, Beta=NA_real_, SE=NA_real_, OR_per_effect=NA_real_,
                          P=NA_real_, N_used=subN, Method="NA", ErrorMessage="Insufficient data/variation"))
        }
        
        dsub <- tibble(Phenotype=pheno01[ok], dosage=dosage[ok], Age=Age[ok], Sex=Sex[ok],
                       PAQ8=PAQ8[ok])
        
        fit <- tryCatch(glm(Phenotype ~ dosage + Age + Sex + PAQ8, data=dsub, family=binomial()), error=function(e) NULL)
        
        if(!is.null(fit)){
            co <- summary(fit)$coefficients
            if("dosage" %in% rownames(co)){
                beta <- co["dosage","Estimate"]; se <- co["dosage","Std. Error"]; p <- co["dosage","Pr(>|z|)"]
                return(tibble(SNP=snp, EffectAllele=eff, Beta=beta, SE=se, OR_per_effect=exp(beta),
                              P=p, N_used=subN, Method="GLM", ErrorMessage=NA_character_))
            }
        }
        
        fitf <- tryCatch(logistf::logistf(Phenotype ~ dosage + Age + Sex + PAQ8, data=dsub), error=function(e) NULL)
        
        if(!is.null(fitf) && "dosage" %in% rownames(fitf$coefficients)){
            beta <- fitf$coefficients["dosage","coef"]; se <- fitf$coefficients["dosage","se"]; p <- fitf$prob["dosage","prob"]
            return(tibble(SNP=snp, EffectAllele=eff, Beta=beta, SE=se, OR_per_effect=exp(beta),
                          P=p, N_used=subN, Method="Firth", ErrorMessage=NA_character_))
        }
        tibble(SNP=snp, EffectAllele=eff, Beta=NA_real_, SE=NA_real_, OR_per_effect=NA_real_,
               P=NA_real_, N_used=subN, Method="Failed", ErrorMessage="Fit error")
    }
    
    if(length(filtered_snps) == 0){
        warning("Subgroup '", label, "': no QC-passed SNPs; skipping.")
        return(invisible(NULL))
    }
    
    assoc_list_sub <- lapply(filtered_snps, run_assoc_one_sub, keep_idx=keep_idx, pheno01=sub_pheno)
    assoc_sub <- bind_rows(assoc_list_sub)
    
    m_eff_sub <- sum(!is.na(assoc_sub$P))
    assoc_sub <- assoc_sub %>%
        mutate(
            FDR_BH      = p.adjust(P, method="BH"),
            Bonferroni = ifelse(is.na(P), NA_real_, pmin(P * m_eff_sub, 1.0))
        )
    
    keep_qc2 <- qc_df %>%
        select(SNP, CallRate, Biallelic, Palindromic,
               MAF_overall, Minor_overall, Major_overall,
               MAF_controls, Minor_controls, Major_controls,
               MAF_cases, Minor_cases, Major_cases,
               HWE_p_controls)
    assoc_sub_out <- assoc_sub %>%
        left_join(keep_qc2, by="SNP") %>%
        left_join(rs_info, by=c("SNP"="RSID")) %>%
        arrange(P)
    
    readr::write_csv(assoc_sub_out, file.path(sub_dir, paste0("assoc_results_", label, ".csv")))
    openxlsx::write.xlsx(assoc_sub_out, file.path(sub_dir, paste0("assoc_results_", label, ".xlsx")), overwrite=TRUE)
    
    assoc_sub_out <- assoc_sub_out %>%
        mutate(
            P_safe = ifelse(is.na(P) | P <= .Machine$double.xmin, .Machine$double.xmin, P),
            log10P = -log10(P_safe),
            DeltaMAF_abs = abs(MAF_cases - MAF_controls)
        )
    
    p_vol <- ggplot(assoc_sub_out, aes(Beta, log10P)) +
        geom_point(size=1.8, alpha=.8) +
        labs(title=paste0("Volcano plot (", label, ")"), x="Beta (per effect allele, adj. Age+Sex+PAQ8)", y="-log10(P)") +
        theme_bw(base_size=12)
    ggsave(file.path(sub_dir, paste0("volcano_", label, ".png")), p_vol, width=8, height=5, dpi=150)
    
    maf_for_x <- ifelse(is.na(assoc_sub_out$MAF_overall), assoc_sub_out$MAF_controls, assoc_sub_out$MAF_overall)
    p_mf <- ggplot(assoc_sub_out, aes(maf_for_x, DeltaMAF_abs)) +
        geom_point(size=1.8, alpha=.8) +
        labs(title=paste0("MAF vs |DeltaMAF| (", label, ")"), x="MAF (overall; fallback controls)", y="|DeltaMAF|") +
        theme_bw(base_size=12)
    ggsave(file.path(sub_dir, paste0("scatter_maf_vs_delta_", label, ".png")), p_mf, width=8, height=5, dpi=150)
    
    assoc_m_sub <- assoc_sub_out %>% filter(!is.na(CHROM), !is.na(START))
    if(nrow(assoc_m_sub) > 0){
        chrom_order <- function(x){
            xx <- as.character(x); xx <- gsub("^chr", "", xx, ignore.case=TRUE)
            key <- suppressWarnings(as.integer(xx))
            ifelse(is.na(key), Inf, key)
        }
        chr_levels <- assoc_m_sub %>% distinct(CHROM) %>% arrange(chrom_order(CHROM)) %>% pull(CHROM)
        assoc_m_sub <- assoc_m_sub %>% mutate(CHROM=factor(CHROM, levels=chr_levels)) %>% arrange(CHROM, START)
        
        offset <- 0; cum_pos <- numeric(nrow(assoc_m_sub)); prev_chr <- NA
        for(i in seq_len(nrow(assoc_m_sub))){
            chr <- assoc_m_sub$CHROM[i]
            if(!identical(chr, prev_chr)){
                if(i>1) offset <- max(cum_pos, na.rm=TRUE) + 1e6
                prev_chr <- chr
            }
            cum_pos[i] <- assoc_m_sub$START[i] + offset
        }
        assoc_m_sub$cum_pos <- cum_pos
        tick_df <- assoc_m_sub %>% group_by(CHROM) %>% summarize(mid=mean(range(cum_pos)), .groups="drop")
        
        bonf_alpha <- 0.05
        bonf_thr   <- ifelse(m_eff_sub > 0, bonf_alpha / m_eff_sub, NA_real_)
        bonf_y     <- ifelse(is.na(bonf_thr), NA_real_, -log10(bonf_thr))
        
        p_manh <- ggplot(assoc_m_sub, aes(cum_pos, log10P)) +
            geom_point(size=1.3, alpha=.85) +
            scale_x_continuous(breaks=tick_df$mid, labels=as.character(tick_df$CHROM)) +
            labs(title=paste0("Manhattan-like (", label, ")"), x="Chromosome", y="-log10(P)") +
            theme_bw(base_size=12)
        if(!is.na(bonf_y) && is.finite(bonf_y)){
            p_manh <- p_manh + geom_hline(yintercept = bonf_y, linetype = 2)
        }
        ggsave(file.path(sub_dir, paste0("manhattan_candidate_snps_", label, ".png")),
               p_manh, width=10, height=4.5, dpi=150)
    }
    
    qq_df_sub <- assoc_sub_out %>% filter(!is.na(P)) %>% arrange(P) %>%
        mutate(obs = -log10(pmax(P, .Machine$double.xmin)),
               exp = -log10(ppoints(n())))
    lambda_gc_sub <- {
        chisq_obs <- qchisq(1 - pmin(assoc_sub_out$P, 1 - 1e-16), df = 1)
        median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
    }
    p_qq_sub <- ggplot(qq_df_sub, aes(exp, obs)) +
        geom_abline(slope=1, intercept=0, linetype=2, size=0.5) +
        geom_point(size=1.6, alpha=.8) +
        labs(title=sprintf("QQ plot (%s, lambdaGC=%.3f)", label, lambda_gc_sub),
             x="Expected -log10(P)", y="Observed -log10(P)") +
        theme_bw(base_size=12)
    ggsave(file.path(sub_dir, paste0("qq_plot_", label, ".png")), p_qq_sub, width=6.5, height=6.5, dpi=150)
    
    top_fdr_sub <- assoc_sub_out %>% filter(!is.na(FDR_BH), FDR_BH < 0.05) %>%
        arrange(P) %>%
        select(SNP, CHROM, START, Beta, OR_per_effect, P, FDR_BH, Bonferroni,
               MAF_controls, MAF_cases, DeltaMAF_abs, N_used, Palindromic)
    top_bonf_sub <- assoc_sub_out %>% filter(!is.na(Bonferroni), Bonferroni < 0.05) %>%
        arrange(P) %>%
        select(SNP, CHROM, START, Beta, OR_per_effect, P, FDR_BH, Bonferroni,
               MAF_controls, MAF_cases, DeltaMAF_abs, N_used, Palindromic)
    readr::write_csv(top_fdr_sub,  file.path(sub_dir, paste0("top_signals_FDR_lt_0.05_", label, ".csv")))
    readr::write_csv(top_bonf_sub, file.path(sub_dir, paste0("top_signals_Bonferroni_lt_0.05_", label, ".csv")))
    
    cat("\n[Subgroup:", label, "] N_ctrl=", sum(sub_pheno==0, na.rm=TRUE),
        " N_case=", sum(sub_pheno==1, na.rm=TRUE),
        " m_eff=", m_eff_sub,
        "\n  ->", file.path(sub_dir, paste0("assoc_results_", label, ".csv")), "\n")
}

run_subgroup("Drusen", 2, 5)
run_subgroup("Pachychoroid", 6, 9)

# ==========================================================
# 11d) Head-to-Head: Drusen-driven (2–5) vs Pachychoroid-driven (6–9)
# ==========================================================

h2h_dir <- file.path(out_dir, "subgroups", "HeadToHead_2to5_vs_6to9")
if(!dir.exists(h2h_dir)) dir.create(h2h_dir, recursive = TRUE)

if(!exists("gt_norm_filtered") || !exists("final_keep")){
    if(length(filtered_snps) > 0){
        gt_norm_filtered <- gt_norm[, filtered_snps, drop = FALSE]
        sample_callrate  <- 1 - rowMeans(is.na(gt_norm_filtered))
        SAMPLE_CALLRATE_MIN <- get0("SAMPLE_CALLRATE_MIN", ifnotfound = 0.98)
        
        final_keep <- which(
            sample_callrate >= SAMPLE_CALLRATE_MIN &
                !is.na(df_a$Phenotype) & !is.na(df_a$Age) & !is.na(df_a$Sex) &
                !is.na(df_a$PAQ8) 
        )
    } else {
        final_keep <- integer(0)
    }
}

idx_dru  <- which(df_a$`AMD sorting` >= 2 & df_a$`AMD sorting` <= 5)
idx_pach <- which(df_a$`AMD sorting` >= 6 & df_a$`AMD sorting` <= 9)

keep_idx_h2h <- intersect(final_keep, union(idx_dru, idx_pach))

Subtype01_all <- rep(NA_real_, nrow(df_a))
Subtype01_all[idx_dru]  <- 0
Subtype01_all[idx_pach] <- 1
Subtype01 <- Subtype01_all[keep_idx_h2h]

if(length(keep_idx_h2h) < 20 ||
   sum(Subtype01==0, na.rm=TRUE) < 10 ||
   sum(Subtype01==1, na.rm=TRUE) < 10){
    warning("Head-to-Head: insufficient samples after filters; skipping.")
} else {
    
    h2h_df <- df_a[keep_idx_h2h, , drop=FALSE] %>%
        mutate(Subtype = factor(ifelse(Subtype01==0, "Drusen(2–5)", "Pachy(6–9)"),
                                levels=c("Drusen(2–5)","Pachy(6–9)")))
    
    h2h_age <- h2h_df %>%
        group_by(Subtype) %>%
        summarize(
            N = n(),
            Age_mean = mean(Age, na.rm=TRUE),
            Age_sd   = sd(Age, na.rm=TRUE),
            Age_median = median(Age, na.rm=TRUE),
            Age_min  = min(Age, na.rm=TRUE),
            Age_max  = max(Age, na.rm=TRUE),
            .groups="drop"
        )
    readr::write_csv(h2h_age, file.path(h2h_dir, "age_summary.csv"))
    
    sx_tab <- table(h2h_df$Sex, h2h_df$Subtype, useNA="no")
    h2h_sex <- as_tibble(sx_tab, .name_repair="minimal")
    colnames(h2h_sex) <- c("Sex","Subtype","Count")
    h2h_sex <- h2h_sex %>% tidyr::pivot_wider(names_from=Subtype, values_from=Count, values_fill=0)
    readr::write_csv(h2h_sex, file.path(h2h_dir, "sex_distribution.csv"))
    
    paq_tab <- table(h2h_df$PAQ8, h2h_df$Subtype, useNA="no")
    h2h_paq8 <- as_tibble(paq_tab, .name_repair="minimal")
    colnames(h2h_paq8) <- c("PAQ8","Subtype","Count")
    h2h_paq8 <- h2h_paq8 %>% tidyr::pivot_wider(names_from=Subtype, values_from=Count, values_fill=0)
    readr::write_csv(h2h_paq8, file.path(h2h_dir, "PAQ8_distribution.csv"))
    
    
    run_h2h_one <- function(snp){
        eff <- effect_lookup[[snp]]     
        gt_all <- gt_norm[[snp]]
        
        g_dru  <- gt_all[intersect(keep_idx_h2h, idx_dru)]
        g_pach <- gt_all[intersect(keep_idx_h2h, idx_pach)]
        
        maf_dru  <- safe_maf(g_dru)  %>% rename(MAF=MAF, Minor=MinorAllele, Major=MajorAllele,
                                                N_alleles=N_alleles, MAC=MAC)
        maf_pach <- safe_maf(g_pach) %>% rename(MAF=MAF, Minor=MinorAllele, Major=MajorAllele,
                                                N_alleles=N_alleles, MAC=MAC)
        
        dosage_all <- dosage_from_effect(gt_all, eff)
        Age    <- df_a$Age
        Sex    <- df_a$Sex
        PAQ8   <- df_a$PAQ8 
        
        ok <- keep_idx_h2h[!is.na(Subtype01) &
                               !is.na(dosage_all[keep_idx_h2h]) &
                               !is.na(Age[keep_idx_h2h]) &
                               !is.na(Sex[keep_idx_h2h]) &
                               !is.na(PAQ8[keep_idx_h2h])] 
        
        if(length(ok) < 20 || length(unique(dosage_all[ok])) < 2 || length(unique(Subtype01[ok])) < 2){
            beta <- se <- p_glm <- OR <- NA_real_; method <- "NA"; err <- "Insufficient data/variation"
        } else {
            
            dsub <- tibble(
                Subtype01 = Subtype01_all[ok],
                dosage = dosage_all[ok],
                Age = Age[ok],
                Sex = Sex[ok],
                PAQ8 = PAQ8[ok] 
            )
            
            fit <- tryCatch(glm(Subtype01 ~ dosage + Age + Sex + PAQ8, data=dsub, family=binomial()),
                            error=function(e) NULL)
            
            if(!is.null(fit)){
                co <- summary(fit)$coefficients
                if("dosage" %in% rownames(co)){
                    beta <- co["dosage","Estimate"]; se <- co["dosage","Std. Error"]
                    p_glm <- co["dosage","Pr(>|z|)"]; OR <- exp(beta); method <- "GLM"; err <- NA_character_
                } else {
                    beta <- se <- p_glm <- OR <- NA_real_; method <- "GLM_failed"; err <- "No dosage coef"
                }
            } else {
                
                fitf <- tryCatch(logistf::logistf(Subtype01 ~ dosage + Age + Sex + PAQ8, data=dsub),
                                 error=function(e) NULL)
                
                if(!is.null(fitf) && "dosage" %in% rownames(fitf$coefficients)){
                    beta <- fitf$coefficients["dosage","coef"]; se <- fitf$coefficients["dosage","se"]
                    p_glm <- fitf$prob["dosage","prob"]; OR <- exp(beta); method <- "Firth"; err <- NA_character_
                } else {
                    beta <- se <- p_glm <- OR <- NA_real_; method <- "Failed"; err <- "Fit error"
                }
            }
        }
        
        alle_tab <- function(gv, eff){
            a <- allele_counts_from_gt(gv)
            eff_n <- unname(a[eff])
            noneff_n <- sum(a, na.rm=TRUE) - eff_n
            c(eff_n, noneff_n)
        }
        a_dru  <- alle_tab(g_dru,  eff);  a_pach <- alle_tab(g_pach, eff)
        m <- rbind(Drusen=a_dru, Pachy=a_pach)
        colnames(m) <- c("EffectAllele","NonEffect")
        p_fisher <- tryCatch(fisher.test(m)$p.value, error=function(e) NA_real_)
        OR_unadj <- tryCatch((m[2,1]*m[1,2]) / (m[1,1]*m[2,2]), error=function(e) NA_real_)  
        
        tibble(
            SNP=snp, EffectAllele=eff,
            Beta_adj=beta, SE_adj=se, OR_per_effect_adj=OR, P_adj=p_glm, Method_adj=method, ErrorMessage=err,
            Allele_OR_unadj=OR_unadj, P_fisher_allele=p_fisher,
            MAF_Drusen=maf_dru$MAF, Minor_Drusen=maf_dru$Minor, Major_Drusen=maf_dru$Major,
            N_alleles_Drusen=maf_dru$N_alleles, MAC_Drusen=maf_dru$MAC,
            MAF_Pachy=maf_pach$MAF, Minor_Pachy=maf_pach$Minor, Major_Pachy=maf_pach$Major,
            N_alleles_Pachy=maf_pach$N_alleles, MAC_Pachy=maf_pach$MAC
        )
    }
    
    if(length(filtered_snps) == 0){
        warning("Head-to-Head: no QC-passed SNPs; skipping.")
    } else {
        h2h_list <- lapply(filtered_snps, run_h2h_one)
        h2h_assoc <- bind_rows(h2h_list)
        
        m_eff_h2h <- sum(!is.na(h2h_assoc$P_adj))
        h2h_assoc <- h2h_assoc %>%
            mutate(
                FDR_BH_adj      = p.adjust(P_adj, method="BH"),
                Bonferroni_adj = ifelse(is.na(P_adj), NA_real_, pmin(P_adj * m_eff_h2h, 1.0))
            )
        
        keep_qc3 <- qc_df %>%
            select(SNP, CallRate, Biallelic, Palindromic,
                   MAF_overall, Minor_overall, Major_overall,
                   MAF_controls, Minor_controls, Major_controls,
                   MAF_cases, Minor_cases, Major_cases,
                   HWE_p_controls)
        h2h_out <- h2h_assoc %>%
            left_join(keep_qc3, by="SNP") %>%
            left_join(rs_info, by=c("SNP"="RSID")) %>%
            arrange(P_adj)
        
        readr::write_csv(h2h_out, file.path(h2h_dir, "assoc_results_head2head.csv"))
        openxlsx::write.xlsx(h2h_out, file.path(h2h_dir, "assoc_results_head2head.xlsx"), overwrite=TRUE)
        
        h2h_out <- h2h_out %>%
            mutate(
                P_safe = ifelse(is.na(P_adj) | P_adj <= .Machine$double.xmin, .Machine$double.xmin, P_adj),
                log10P = -log10(P_safe),
                DeltaMAF_cases = abs(MAF_Pachy - MAF_Drusen)
            )
        
        p_vol_h2h <- ggplot(h2h_out, aes(x=Beta_adj, y=log10P)) +
            geom_point(size=1.8, alpha=.8) +
            labs(title="Volcano: Pachychoroid(6–9) vs Drusen(2–5)",
                 x="Beta (per controls' minor, adj. Age+Sex+PAQ8)", y="-log10(P, adjusted model)") +
            theme_bw(base_size=12)
        ggsave(file.path(h2h_dir, "volcano_head2head.png"), p_vol_h2h, width=8, height=5, dpi=150)
        
        p_dmaf_h2h <- ggplot(h2h_out, aes(x=MAF_overall, y=DeltaMAF_cases)) +
            geom_point(size=1.8, alpha=.8) +
            labs(title="MAF vs |DeltaMAF| between case subtypes",
                 x="MAF (overall; fallback controls)", y="|DeltaMAF| (Pachy - Drusen)") +
            theme_bw(base_size=12)
        ggsave(file.path(h2h_dir, "scatter_maf_vs_delta_head2head.png"), p_dmaf_h2h, width=8, height=5, dpi=150)
        
        h2h_m <- h2h_out %>% filter(!is.na(CHROM), !is.na(START))
        if(nrow(h2h_m) > 0){
            chrom_order <- function(x){
                xx <- as.character(x); xx <- gsub("^chr","",xx,ignore.case=TRUE)
                key <- suppressWarnings(as.integer(xx))
                ifelse(is.na(key), Inf, key)
            }
            chr_levels <- h2h_m %>% distinct(CHROM) %>% arrange(chrom_order(CHROM)) %>% pull(CHROM)
            h2h_m <- h2h_m %>% mutate(CHROM=factor(CHROM, levels=chr_levels)) %>% arrange(CHROM, START)
            
            offset <- 0; cum_pos <- numeric(nrow(h2h_m)); prev_chr <- NA
            for(i in seq_len(nrow(h2h_m))){
                chr <- h2h_m$CHROM[i]
                if(!identical(chr, prev_chr)){
                    if(i>1) offset <- max(cum_pos, na.rm=TRUE) + 1e6
                    prev_chr <- chr
                }
                cum_pos[i] <- h2h_m$START[i] + offset
            }
            h2h_m$cum_pos <- cum_pos
            tick_df <- h2h_m %>% group_by(CHROM) %>% summarize(mid=mean(range(cum_pos)), .groups="drop")
            
            bonf_alpha <- 0.05
            bonf_thr   <- ifelse(m_eff_h2h > 0, bonf_alpha / m_eff_h2h, NA_real_)
            bonf_y     <- ifelse(is.na(bonf_thr), NA_real_, -log10(bonf_thr))
            
            p_manh_h2h <- ggplot(h2h_m, aes(cum_pos, log10P)) +
                geom_point(size=1.3, alpha=.85) +
                scale_x_continuous(breaks=tick_df$mid, labels=as.character(tick_df$CHROM)) +
                labs(title="Manhattan-like: Head-to-Head (6–9 vs 2–5)", x="Chromosome", y="-log10(P)") +
                theme_bw(base_size=12)
            if(!is.na(bonf_y) && is.finite(bonf_y)){
                p_manh_h2h <- p_manh_h2h + geom_hline(yintercept = bonf_y, linetype = 2)
            }
            ggsave(file.path(h2h_dir, "manhattan_head2head.png"), p_manh_h2h, width=10, height=4.5, dpi=150)
        }
        
        qq_h2h <- h2h_out %>% filter(!is.na(P_adj)) %>% arrange(P_adj) %>%
            mutate(obs = -log10(pmax(P_adj, .Machine$double.xmin)),
                   exp = -log10(ppoints(n())))
        lambda_gc_h2h <- {
            chisq_obs <- qchisq(1 - pmin(h2h_out$P_adj, 1 - 1e-16), df = 1)
            median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
        }
        p_qq_h2h <- ggplot(qq_h2h, aes(exp, obs)) +
            geom_abline(slope=1, intercept=0, linetype=2, size=0.5) +
            geom_point(size=1.6, alpha=.8) +
            labs(title=sprintf("QQ plot (Head-to-Head, lambdaGC = %.3f)", lambda_gc_h2h),
                 x="Expected -log10(P)", y="Observed -log10(P)") +
            theme_bw(base_size=12)
        ggsave(file.path(h2h_dir, "qq_plot_head2head.png"), p_qq_h2h, width=6.5, height=6.5, dpi=150)
        
        top_fdr_h2h <- h2h_out %>% filter(!is.na(FDR_BH_adj), FDR_BH_adj < 0.05) %>%
            arrange(P_adj) %>%
            select(SNP, CHROM, START, Beta_adj, OR_per_effect_adj, P_adj, FDR_BH_adj, Bonferroni_adj,
                   MAF_Drusen, MAF_Pachy, DeltaMAF_cases, Allele_OR_unadj, P_fisher_allele,
                   Palindromic, CallRate)
        top_bonf_h2h <- h2h_out %>% filter(!is.na(Bonferroni_adj), Bonferroni_adj < 0.05) %>%
            arrange(P_adj) %>%
            select(SNP, CHROM, START, Beta_adj, OR_per_effect_adj, P_adj, FDR_BH_adj, Bonferroni_adj,
                   MAF_Drusen, MAF_Pachy, DeltaMAF_cases, Allele_OR_unadj, P_fisher_allele,
                   Palindromic, CallRate)
        readr::write_csv(top_fdr_h2h,  file.path(h2h_dir, "top_signals_FDR_lt_0.05_head2head.csv"))
        readr::write_csv(top_bonf_h2h, file.path(h2h_dir, "top_signals_Bonferroni_lt_0.05_head2head.csv"))
        
        cat("\n[Head-to-Head 2–5 vs 6–9]  N_Drusen=", sum(Subtype01==0, na.rm=TRUE),
            " N_Pachy=", sum(Subtype01==1, na.rm=TRUE),
            " m_eff=", m_eff_h2h,
            "\n ->", file.path(h2h_dir, "assoc_results_head2head.csv"), "\n")
    }
}

# ---------- 12) PDF report ----------
pdf_path <- file.path(out_dir,"AMD_candidate_GWAS_report.pdf")
pdf(pdf_path, width=11, height=8.5)
grid::grid.newpage()
grid::grid.text("AMD Candidate-SNP Association Report", x=.5, y=.8, gp=grid::gpar(cex=1.6))
grid::grid.text("Age >= 50; Logistic regression adjusted for Age, Sex & PAQ8\nQC: CallRate>=0.98, MAF>=0.01 (controls), HWE p(controls)>=1e-6, biallelic",
                x=.5, y=.65, gp=grid::gpar(cex=1.0))

grid::grid.newpage()
gridExtra::grid.table(round(cohort_summary,3))
grid::grid.newpage()
gridExtra::grid.table(filter_summary)

if(exists("final_age_summary") && nrow(final_age_summary) > 0){
    grid::grid.newpage(); gridExtra::grid.table(final_age_summary)
}
if(exists("final_sex_table") && nrow(final_sex_table) > 0){
    grid::grid.newpage(); gridExtra::grid.table(final_sex_table)
}
if(exists("final_paq8_table") && nrow(final_paq8_table) > 0){
    grid::grid.newpage(); gridExtra::grid.table(final_paq8_table)
}
if(exists("final_tests") && nrow(final_tests) > 0){
    grid::grid.newpage(); gridExtra::grid.table(final_tests)
}


grid::grid.newpage(); grid::grid.raster(as.raster(png::readPNG(file.path(out_dir,"manhattan_candidate_snps.png"))))
grid::grid.newpage(); grid::grid.raster(as.raster(png::readPNG(file.path(out_dir,"volcano_beta_log10p.png"))))
grid::grid.newpage(); grid::grid.raster(as.raster(png::readPNG(file.path(out_dir,"scatter_maf_vs_delta.png"))))
grid::grid.newpage(); grid::grid.raster(as.raster(png::readPNG(file.path(out_dir,"qq_plot.png"))))
if(file.exists(file.path(out_dir,"LD_heatmap_Bonferroni_sig_r2.png"))){
    grid::grid.newpage(); grid::grid.raster(as.raster(png::readPNG(file.path(out_dir,"LD_heatmap_Bonferroni_sig_r2.png"))))
}


if(nrow(top_fdr)>0){
    grid::grid.newpage(); gridExtra::grid.table(head(top_fdr, 25))
}
if(nrow(top_bonf)>0){
    grid::grid.newpage(); gridExtra::grid.table(head(top_bonf, 25))
}
if(nrow(palin_check)>0){
    grid::grid.newpage(); gridExtra::grid.table(head(palin_check, 30))
}
if(nrow(cluster_leads)>0){
    grid::grid.newpage(); gridExtra::grid.table(head(cluster_leads, 30))
}

if(file.exists(file.path(out_dir,"qq_plot_LDpruned_cluster_leads.png"))){
    grid::grid.newpage(); grid::grid.raster(as.raster(
        png::readPNG(file.path(out_dir,"qq_plot_LDpruned_cluster_leads.png"))
    ))
}
dev.off() 

# ---------- 13) Console summary ----------
cat("=== DONE ===\nOutput dir:", normalizePath(out_dir), "\n",
    "Assoc results  :", file.path(out_dir,"assoc_results_filtered.csv"), "\n",
    "Manhattan plot :", file.path(out_dir,"manhattan_candidate_snps.png"), "\n",
    "Report PDF      :", pdf_path, "\n", sep="")