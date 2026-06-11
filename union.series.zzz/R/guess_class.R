
# ==========================================================================
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.general_prefix_of_class <- function() {
  prefix <- c(srn = "job_seurat5n",
    metadata = "data.frame",
    sr = "job_seurat",
    mr = "job_mr",
    cc = "job_cellchat",
    ssg = "job_ssgsea",
    lm = "job_limma",
    ml = "job_mlearn",
    des = "job_deseq2",
    ven = "job_venn",
    sce = "job_scenic",
    vn = "job_vina",
    mn = "job_monocle2",
    rt = "job_reactome",
    rn = "job_regNet",
    gb = "job_gBan",
    mb = "job_mebocost",
    sc = "job_scenic",
    hd = "job_hdwgcna",
    hdw = "job_hdwgcna",
    wgc = "jop_wgcna",
    ssr = "job_scissor",
    au = "job_aucell",
    rms = "job_rms",
    en = "job_enrich",
    sct = "job_scteni",
    geo = "job_geo",
    mdb = "job_msigdb",
    fea = "feature",
    corgsea = "job_corgsea",
    iobr = "job_iobr",
    gn = "job_genecard",
    mti = "job_metaInte",
    mdiff = "job_metaboDiff"
  )
  if (any(duplicated(names(prefix)))) {
    stop('any(duplicated(names(prefix))), in `.general_prefix_of_class`')
  }
  prefix
}

