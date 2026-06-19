# ==========================================================================
# workflow of vina
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (!exists("vinaFuns")) {
  vinaFuns <- new.env(parent = emptyenv())
}

.job_vina <- setClass("job_vina", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorials: https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"),
    cite = "[@AutodockVina1Eberha2021; @AutogridfrImpZhang2019; @AutodockCrankpZhang2019; @AutositeAnAuRavind2016; @AutodockfrAdvRavind2015]",
    method = "The CLI tools of `AutoDock vina` and `ADFR` software used for auto molecular docking",
    tag = "dock:vina",
    analysis = "AutoDock vina 分子对接"
    ))

job_vina <- function(cids, hgnc_symbols, .layout = NULL)
{
  if (missing(cids) || missing(hgnc_symbols)) {
    message("Get from `.layout` first 3 columns: cpd names, hgnc_symbols, cids")
    cids <- nl(.layout[[ 1 ]], .layout[[ 3 ]], as.list = FALSE)
    cids <- cids[!duplicated(cids)]
    hgnc_symbols <- .layout[[ 2 ]]
    hgnc_symbols <- hgnc_symbols[ !duplicated(hgnc_symbols) ]
  }
  x <- .job_vina(object = namel(cids, hgnc_symbols))
  x$.layout <- .layout
  x
}

.select_pdb <- setClass("select_pdb",
  representation = representation(
    id = "character", chain = "character", resi = "integer"
  ),
  prototype = NULL)

select_pdb <- function(id, chain, resi) {
  .select_pdb(id = id, chain = chain, resi = resi)
}

setGeneric("asjob_vina",
  function(x, ...) standardGeneric("asjob_vina"))

setMethod("asjob_vina", signature = c(x = "job_stringdb"),
  function(x, cids, job_herb = NULL, compounds = NULL, hubs = 10)
  {
    hgnc_symbols <- head(x@tables$step1$hub_genes$hgnc_symbol, n = hubs)
    if (!is.null(job_herb)) {
      compounds_targets <- dplyr::filter(job_herb@tables$step2$compounds_targets,
        Target.name %in% dplyr::all_of(hgnc_symbols))
      compounds <- dplyr::filter(object(job_herb)$component,
        Ingredient_id %in% compounds_targets$Ingredient_id)
      cids <- compounds$PubChem_id
      cids <- cids[ !is.na(cids) ]
      from_job_herb <- TRUE
    } else {
      from_job_herb <- NULL
    }
    x <- job_vina(cids, hgnc_symbols)
    x$compounds <- compounds
    x$from_job_herb <- from_job_herb
    return(x)
  })

setMethod("step0", signature = c(x = "job_vina"),
  function(x){
    step_message("Prepare your data with function `job_vina`. ")
  })

setMethod("step1", signature = c(x = "job_vina"),
  function(x, order = TRUE, each_target = 1, custom_pdbs = NULL, 
    bdb_file = .prefix("BindingDB_All_202401.tsv", "db"), 
    forceAF = FALSE, exclude_pdb = NULL, recode = NULL)
  {
    step_message("Prepare Docking Combination.")
    x <- .find_proper_pdb(
      x, order, each_target, custom_pdbs, exclude_pdb, 
      forceAF = forceAF, recode = recode
    )
    x$dock_layout <- sapply(object(x)$cids, function(cid) x$used_pdbs, simplify = FALSE)
    names(x$dock_layout) <- object(x)$cids
    return(x)
  })

.find_proper_pdb <- function(x, order = TRUE, each_target = 1, 
  custom_pdbs = NULL, exclude_pdb, forceAF = FALSE, recode = NULL)
{
  if (!is(x, "job")) {
    stop('!is(x, "job").')
  }
  if (!is(object(x)$hgnc_symbols, "character")) {
    stop('!is(object(x)$hgnc_symbols, "character").')
  }
  if (is.null(x$mart)) {
    mart <- new_biomart()
  } else {
    mart <- x$mart
  }
  if (!forceAF && is.null(x$targets_annotation)) {
    genes <- object(x)$hgnc_symbols
    if (!is.null(recode)) {
      genes <- dplyr::recode(genes, !!!recode)
    }
    x$targets_annotation <- filter_biomart(
      mart, c("hgnc_symbol", "pdb"), "hgnc_symbol",
      genes, distinct = FALSE
    )
    if (!is.null(recode)) {
      x$targets_annotation <- dplyr::mutate(
        x$targets_annotation, hgnc_symbol = dplyr::recode(
          hgnc_symbol, !!!setNames(names(recode), unname(recode))
        )
      )
      genes <- object(x)$hgnc_symbols
    }
    if (nrow(x$targets_annotation)) {
      x <- methodAdd(x, "以 R 包 `biomaRt` ⟦pkgInfo('biomaRt')⟧ {cite_show('MappingIdentifDurinc2009')} 获取基因 Symbol 对应的蛋白结构 PDB (<https://www.rcsb.org/>) 数据库 ID。")
      # x <- snapAdd(x, "以 `biomaRt` 获取基因 Symbol 对应的蛋白结构 (PDB，详见方法章节)。")
      x$targets_annotation <- dplyr::filter(x$targets_annotation, pdb != "")
      x$targets_annotation <- dplyr::distinct(x$targets_annotation, pdb, .keep_all = TRUE)
    }
  }
  if (!is.null(custom_pdbs)) {
    if (is.null(x$targets_annotation)) {
      x$targets_annotation <- data.frame()
    }
    custom_pdbs <- list(hgnc_symbol = names(custom_pdbs), pdb = unname(custom_pdbs))
    if (!is.character(x$targets_annotation$pdb)) {
      x$targets_annotation <- dplyr::mutate(x$targets_annotation, pdb = as.character(pdb))
    }
    x$targets_annotation <- tibble::add_row(x$targets_annotation, !!!custom_pdbs)
  }
  if (any(isThat <- !object(x)$hgnc_symbols %in% x$targets_annotation$hgnc_symbol)) {
    x$pdb_notGot <- unique(object(x)$hgnc_symbols[ isThat ])
    message("PDB not found:\n\t", paste0(x$pdb_notGot, collapse = ", "))
    if (!sureThat("Continue? (will retrieve from AlphaFold in step3 by `.get_pdb_files`)"))
    {
      stop("Consider other ways to found PDB files for docking.")
    }
  } else {
    message("Got PDB for all `hgnc_symbol`.")
  }
  if (!forceAF && nrow(x$targets_annotation)) {
    x$annoPdbs <- as_tibble(e(bio3d::pdb.annotate(x$targets_annotation$pdb)))
    x <- methodAdd(
      x, "以 R 包 `bio3d` ⟦pkgInfo('bio3d')⟧ 获取 PDB ID 对应的注释 (蛋白结构分辨率, resolution) 。"
    )
    if (order) {
      annoPdbs <- dplyr::distinct(x$annoPdbs, structureId, resolution)
      x$targets_annotation <- map(
        x$targets_annotation, "pdb", annoPdbs, "structureId", 
        "resolution", col = "resolution"
      )
      x$targets_annotation <- dplyr::arrange(x$targets_annotation, resolution)
      x <- methodAdd(x, "首要以 resolution 选取用于分子对接的蛋白结构 (resolution 越小，分辨率越高) 。")
      # x <- snapAdd(x, "选取分辨率最高 (即，resolution 值最小) 的 PDB 作为分子对接的蛋白结构。")
    }
    if (!is.null(exclude_pdb)) {
      message("Exclude pdb: ", bind(exclude_pdb))
      x$targets_annotation <- dplyr::filter(x$targets_annotation, !pdb %in% !!exclude_pdb)
    }
    used_pdbs <- dplyr::distinct(
      x$targets_annotation, hgnc_symbol, .keep_all = TRUE
    )
    isMultiChains <- vapply(used_pdbs$pdb,
      function(pdb) {
        length(which(x$annoPdbs$structureId == pdb)) > 1
      }, logical(1))
    if (any(isMultiChains)) {
      x$pdb_MultiChains <- used_pdbs$pdb[ isMultiChains ]
      message(glue::glue("Some 'pdb' has multiple chains ({bind(x$pdb_MultiChains)})"))
    }
    used_pdbs <- setLegend(used_pdbs, "基因 Symbol 所用的蛋白结构，以及对应的 PDB ID 和分辨率。")
    x <- tablesAdd(x, t.proteins_used_PDB = used_pdbs)
    x$used_pdbs <- used_pdbs <- nl(
      used_pdbs$hgnc_symbol, used_pdbs$pdb, FALSE
    )
    if (!length(used_pdbs)) {
      stop('!length(used_pdbs), no any protein need docking?')
    }
  } else {
    x$used_pdbs <- NULL
  }
  return(x)
}

setMethod("step2", signature = c(x = "job_vina"),
  function(x, try_cluster_random = FALSE, nGroup = 100, 
    nMember = 3, cl = 5, sdf.3d = NULL, dir_save = paste0(x@sig, "_cpd"),
    conda_env = "base", use_pubchem_3d = TRUE, use_obgen = TRUE,
    rdkit_fallback = TRUE, obabel_fallback = TRUE,
    strip_salts = TRUE, neutralize_ligand = FALSE,
    ligand_forcefield = c("MMFF", "UFF"), ligand_minimize = TRUE,
    ligand_per_molecule = TRUE, obabel_partialcharge = "gasteiger",
    overwrite_sdf = FALSE, overwrite_pubchem_3d = FALSE,
    overwrite_obgen = FALSE)
  {
    step_message("Download sdf files and convert as pdbqt for ligands.")
    ligand_forcefield <- match.arg(ligand_forcefield)
    if (!is.null(conda_env)) {
      activate_env(conda_env)
    }
    input_cids <- unique(as.character(names(x$dock_layout)))
    dir_sdf <- paste0(dir_save, "_SDF")
    file_sdf <- file.path(
      dir_sdf, add_filename_suffix("all_compounds.sdf", x@sig)
    )
    if (!overwrite_sdf && file.exists(file_sdf) && file.size(file_sdf) > 0L) {
      message("Use existing 2D SDF: ", file_sdf)
      sdfFile <- file_sdf
    } else {
      sdfFile <- query_sdfs(
        input_cids,
        dir_sdf,
        curl_cl = cl, filename = add_filename_suffix("all_compounds.sdf", x@sig)
      )
    }
    x$sdf_2d_file <- sdfFile
    x <- methodAdd(x, "以 PubChem API
      (<https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>) 获取化合物 SDF
      结构文件。"
    )
    Show_filter <- FALSE
    if (try_cluster_random && length(object(x)$cids) > nGroup) {
      message("To reduce docking candidates, clustering the molecules, and random sample each group to get `n`")
      set.seed(100)
      sdfset <- e(ChemmineR::read.SDFset(sdfFile))
      message("Read ", length(sdfset), " molecules.")
      apset <- e(ChemmineR::sdf2ap(sdfset))
      cluster <- e(ChemmineR::cmp.cluster(db = apset, cutoff = seq(0.9, 0.4, by = -0.1)))
      x <- methodAdd(x, "以 R 包 `ChemmineR` ⟦pkgInfo('ChemmineR')⟧ {cite_show('ChemminerACoCaoY2008')} 计算化合物结构相似度。")
      cluster <- dplyr::select(cluster, 1, dplyr::starts_with("CLID"))
      nCl <- apply(cluster[, -1], 2, function(x) length(unique(x)))
      useWhich <- which( nCl < 30 )[1] + 1
      if (is.na(useWhich)) {
        useWhich <- menu(paste0(names(nCl), "__", unname(nCl)), title = "Need custom select for cutoff.")
        useWhich <- useWhich + 1
      }
      message("Use cutoff: ", colnames(cluster)[useWhich])
      groups <- split(seq_along(sdfset), cluster[useWhich])
      message("Now, get group number: ", length(groups))
      groups <- lapply(groups,
        function(x) {
          if (length(x) > nMember) {
            sample(x, nMember)
          } else x
        })
      sdfset <- sdfset[unlist(groups, use.names = FALSE)]
      MaybeError <- vapply(ChemmineR::cid(sdfset), FUN.VALUE = logical(1),
        function(x) {
          dim(sdfset[[x]][[2]])[2] < 3
        })
      sdfset <- sdfset[which(!MaybeError)]
      sdfFile <- add_filename_suffix(sdfFile, "random")
      e(ChemmineR::write.SDF(sdfset, sdfFile))
      input_cids <- intersect(input_cids, as.character(ChemmineR::cid(sdfset)))
      Show_filter <- TRUE
      .add_internal_job(.job(method = "R package `ChemmineR` used for similar chemical compounds clustering",
          cite = "[@ChemminerACoCaoY2008]"
          ))
      x$chem_lich <- new_lich(list(
          "Clustering method:" = "Binning Clustering",
          "Use Cut-off:" = colnames(cluster)[useWhich],
          "Cluster number:" = length(groups),
          "Each sampling for next step (docking):" = nMember
          )
      )
    }
    if (Show_filter) {
      message("Filter out (Due to Chemmine bug): ", length(which(MaybeError)))
      message("Now, Total docking molecules: ", length(sdfset))
    }
    if (!is.null(sdf.3d) && !file.exists(sdf.3d)) {
      stop("file.exists(sdf.3d) == FALSE")
    }
    res.pdbqt <- vinaFuns$prepare_ligand_pdbqt_recovery(
      cids = input_cids,
      sdf_2d = sdfFile,
      sdf_3d = sdf.3d,
      dir_save = dir_save,
      mkdir.pdbqt = paste0(dir_save, "_pdbqt"),
      use_pubchem_3d = use_pubchem_3d,
      use_obgen = use_obgen,
      rdkit_fallback = rdkit_fallback,
      obabel_fallback = obabel_fallback,
      strip_salts = strip_salts,
      neutralize_ligand = neutralize_ligand,
      ligand_forcefield = ligand_forcefield,
      ligand_minimize = ligand_minimize,
      ligand_per_molecule = ligand_per_molecule,
      obabel_partialcharge = obabel_partialcharge,
      overwrite_pubchem_3d = overwrite_pubchem_3d,
      overwrite_obgen = overwrite_obgen,
      cl = cl
    )
    x <- methodAdd(x, "配体结构准备采用多级策略：优先使用 PubChem 已计算的 3D conformer；对于无可用 3D 记录或未能成功转化的化合物，使用 `openbabel` 生成 3D 构象；对于仍未成功的结构，进一步使用 RDKit 进行结构标准化、去除游离盐或反离子并重新生成 3D 构象。随后优先以 Python `meeko` 包 (`mk_prepare_ligand.py`) 转化得到配体 PDBQT 文件；若 `meeko` 对复杂盐型、多片段或电荷结构未能完成转化，则使用 `openbabel` 作为备选 PDBQT 转换工具，并记录每个配体的结构来源和转换工具。")
    x$ligand_prepare_stat <- res.pdbqt$stat
    x$ligand_prepare_status <- res.pdbqt$status
    x$ligand_prepare_sdf_files <- res.pdbqt$sdf_files
    x$res.ligand <- nl(res.pdbqt$pdbqt.cid, res.pdbqt$pdbqt)
    message("Got ligand PDBQT: ", length(x$res.ligand))
    alls <- unique(as.character(input_cids))
    notGot <- alls[!alls %in% names(x$res.ligand)]
    message("Not got: ", paste0(notGot, collapse = ", "))
    if (length(notGot)) {
      x <- methodAdd(x, "SDF 输入分子数量为 {length(alls)}，经多级配体结构准备后输出 PDBQT 的分子数量为 {length(x$res.ligand)}；其余 {length(notGot)} 个化合物因结构标准化、3D 构象生成或 PDBQT 转换未成功而未进入后续分子对接。")
    } else {
      x <- methodAdd(x, "SDF 输入分子数量为 {length(alls)}，经多级配体结构准备后均成功获得 PDBQT 文件。")
    }
    x$ligand_notGot <- notGot
    message("Filter the `x$dock_layout`")
    x$dock_layout <- x$dock_layout[ names(x$dock_layout) %in% names(x$res.ligand) ]
    return(x)
  })

setMethod("step3", signature = c(x = "job_vina"),
  function(x, cl = 10, pattern = NULL,
    # extra_pdb.files: hgnc_symbol = **pdb ID** file
    extra_pdb.files = NULL, extra_layouts = NULL,
    filter = TRUE, use_complex = TRUE, select = .select_pdb(), exclude_nonStd = c(
      "NAG", "BMA", "FUL"
      ), tryAF = TRUE, forceAF = x$forceAF %||% FALSE, split_chain = FALSE,
    tool_prepare = c("prepare_receptor", "mk_prepare_receptor.py"),
    dir_save = paste0(x@sig, "_protein_pdb"),
    path_adfr = getOption("path_adfr"))
  {
    step_message("Dowload pdb files for Receptors.")
    if (!nchar(Sys.which("prepare_receptor")) && !is.null(path_adfr)) {
      message(glue::glue("Add path to system PATH: {path_adfr}"))
      path <- normalizePath(path_adfr)
      Sys.setenv(PATH = paste(path, Sys.getenv("PATH"), sep = ":"))
    }
    res <- .get_pdb_files(
      x, cl = cl, pattern = pattern, extra_pdb.files = extra_pdb.files, 
      extra_layouts = extra_layouts,
      tryAF = tryAF, forceAF = forceAF, split_chain = split_chain, dir_save = dir_save
    )
    x <- res$x
    pdb.files <- res$pdb.files
    used_pdbs <- x$used_pdbs
    dir.create(cpdir <- paste0(dir_save, "_clean"), FALSE)
    pdb.files <- lapply(pdb.files, 
      function(file) {
        newfile <- file.path(cpdir, basename(file))
        .pymol_select_polymer.protein(file, newfile)
        newfile
      })
    # x <- snapAdd(x, "随后，以 `pymol` 仅保留蛋白结构 (polymer.protein)。")
    x <- methodAdd(x, "以 `pymol` 仅保留蛋白结构 (polymer.protein) (去除了原 PDB 中的配体等其他结构)。")
    tool_prepare <- match.arg(tool_prepare)
    x$res.receptor <- prepare_receptor(
      pdb.files, paste0(dir_save, "qt"), use = tool_prepare
    )
    if (tool_prepare == "prepare_receptor") {
      x <- methodAdd(x, "以 `ADFR` {cite_show('AutogridfrImpZhang2019')} 工具组的准备受体蛋白的 PDBQT 文件 (以 `prepare_receptor` 添加氢原子，移除对对接而言不必要的分子水、配体、辅因子和离子等，并转化为 PDBQT 文件) 。请参考 <https://autodock-vina.readthedocs.io/en/latest/docking_basic.html>。")
    } else {
      x <- methodAdd(x, "以 Python `meeko` 包 (`mk_prepare_receptor.py`) 准备受体蛋白的 PDBQT 文件 (添加氢原子，移除对对接而言不必要的分子水、配体、辅因子和离子等，并转化为 PDBQT 文件) 。 ")
    }
    # x <- snapAdd(x, "以 `ADFR` 工具给受体添加氢原子，转化为 PDBQT 文件。")
    names <- names(x$res.receptor)
    gotSymbols <- names(x$used_pdbs)[ match(tolower(names(x$res.receptor)), tolower(x$used_pdbs)) ]
    x$res.receptor.symbol <- gotSymbols
    fun <- function(x) x[ !x %in% gotSymbols ]
    message("Not got: ", paste0(fun(unique(object(x)$hgnc_symbol)), collapse = ", "))
    if (filter) {
      if (!is.null(x$.layout)) {
        message("Customize using `x$.layout` columns: ", paste0(colnames(x$.layout)[1:2], collapse = ", "))
        x <- filter(x, x$.layout[[ 1 ]], x$.layout[[ 2 ]])
      }
    }
    if (any(lengths(x$dock_layout) > 1)) {
      names <- rep(names(x$dock_layout), lengths(x$dock_layout))
      dock_layout <- unlist(x$dock_layout, use.names = FALSE)
      names(dock_layout) <- names
      x$dock_layout <- dock_layout
    }
    if (use_complex && any(isThats <- grpl(names(used_pdbs), "\\+"))) {
      isUsedChains <- vapply(x$dock_layout, 
        function(x) {
          grpl(x, "_")
        }, logical(1))
      cpdsUsedChains <- names(x$dock_layout)[isUsedChains]
      dat <- dplyr::distinct(
        data.frame(
          cpd = cpdsUsedChains, pdb = gs(
            unlist(x$dock_layout[isUsedChains]), "_[a-z]$", ""
          )
        )
      )
      x$dock_layout <- c(x$dock_layout, nl(dat$cpd, dat$pdb))
    }
    if (!length(x$res.receptor) || !length(x$res.ligand)) {
      stop(
        '!length(x$res.receptor) || !length(x$res.ligand), no any ligand or receptor?'
      )
    }
    return(x)
  })

.get_pdb_files <- function(x, cl = 10, pattern = NULL, 
  extra_pdb.files = NULL, extra_layouts = NULL,
  tryAF = TRUE, forceAF = FALSE, split_chain = FALSE, dir_save = "protein_pdb")
{
  if (!is(x, "job")) {
    stop('!is(x, "job").')
  }
  if (is.null(x$dock_layout)) {
    stop('is.null(x$dock_layout).')
  }
  ids <- rm.no(unlist(x$dock_layout, use.names = FALSE))
  if (length(ids) && !forceAF) {
    # The returned files is got by list.files
    pdb.files <- get_pdb(ids, cl = cl, mkdir.pdb = dir_save)
    pdb.files <- pdb.files[ names(pdb.files) %in% tolower(x$used_pdbs) ]
    x <- methodAdd(x, "以 RCSB API  (<https://www.rcsb.org/docs/programmatic-access/web-apis-overview>) 获取蛋白 PDB 文件。")
    # x <- snapAdd(x, "从 RCSB PDB 获取 PDB 文件。")
  } else {
    pdb.files <- NULL
  }
  if (forceAF || (length(x$pdb_notGot) && tryAF)) {
    if (forceAF) {
      dir.create(dir_save, FALSE)
      genes_touch_AF <- unique(x@object$hgnc_symbols)
      x$used_pdbs <- NULL
    } else {
      genes_touch_AF <- x$pdb_notGot
    }
    res_af <- get_pdb_from_alphaFold(genes_touch_AF, dir_save)
    # x <- methodAdd(x, "以 R 包 `UniProt.ws` ⟦pkgInfo('UniProt.ws')⟧ 获取基因 (symbol) 的 `UniProtKB-Swiss-Prot` ID (Entry ID)。")
    if (forceAF || !is.null(extra_pdb.files)) {
      getFromAF <- unique(c(genes_touch_AF, names(extra_pdb.files)))
      x <- methodAdd(
        x, "从数据库 `AlphaFold` (<https://alphafold.ebi.ac.uk/>) 获取 {bind(unique(getFromAF))} 蛋白结构 (已有诸多文献报道使用 alphaFold 数据库提供的蛋白用于虚拟筛选 {cite_show('Bioinformatics_Salama_2025')}，{cite_show('Virtual_Screeni_Wang_2025')} {cite_show('Accurate_struct_Abrams_2024')})。"
      )
    }
    # x <- snapAdd(x, "{if (forceAF) '' else '对于未从 `PDB` 数据库找到结构文件的，'}从数据库 `AlphaFold` 获取 {less(genes_touch_AF)} 预测的蛋白结构 (根据 `UniProtKB-Swiss-Prot` ID，详见方法章节)。")
    x$pdb_notGot_uniprot <- res_af$info
    # extra_pdb.files: hgnc_symbol = file
    if (!is.null(extra_pdb.files)) {
      # symbol = pdb_ID
      customInput <- basename(tools::file_path_sans_ext(extra_pdb.files))
      names(customInput) <- names(extra_pdb.files)
      customFiles <- extra_pdb.files
      # pdb_ID = files
      names(customFiles) <- unname(customInput)
    } else {
      customInput <- NULL
      customFiles <- NULL
    }
    extra_pdb.files <- c(
      res_af$files, customFiles
    )
    used_pdbs_extra <- c(res_af$used_pdbs, customInput)
    # need revision: x$dock_layout, pdb.files, x$used_pdbs
  }
  if (!is.null(extra_pdb.files)) {
    pdb.files <- c(pdb.files, extra_pdb.files)
    x$used_pdbs <- c(used_pdbs_extra, x$used_pdbs)
    layoutNeedRevise <- TRUE
  } else {
    layoutNeedRevise <- FALSE
  }
  if (!is.null(extra_layouts)) {
    if (is(extra_layouts, "character")) {
      extra_layouts <- as.list(extra_layouts)
    }
    x$dock_layout <- c(x$dock_layout, extra_layouts)
  } else if (layoutNeedRevise) {
    x$dock_layout <- lapply(
      x$dock_layout, function(x) c(x, used_pdbs_extra)
    )
  }
  if (!is.null(x$pdb_MultiChains) && !forceAF && split_chain) {
    pdb_MultiChains <- tolower(unique(x$pdb_MultiChains))
    fun_extract_cpdName <- function(lines) {
      chains <- stringr::str_extract(lines, "(?<=CHAIN: ).*(?=;)")
      names <- stringr::str_extract(lines, "(?<=MOLECULE: ).*(?=;)")
      names <- names[ !is.na(names) ]
      names <- if (length(names)) names else "__Omit__"
      names(names) <- chains[ !is.na(chains) ]
      return(names)
    }
    used_pdbs <- x$used_pdbs
    pdb_files_MultiChains <- lapply(pdb_MultiChains,
      function(pdb) {
        file <- pdb.files[[ pdb ]]
        contentPdb <- readLines(file)
        contentPdb <- contentPdb[ grpl(contentPdb, "^COMPND") ]
        contentPdb <- sep_list(contentPdb, "^COMPND.*MOL_ID", 0)
        if (length(contentPdb) > 1) {
          chainNames <- unlist(lapply(unname(contentPdb), fun_extract_cpdName))
          message(glue::glue("PDB: {pdb.files[[ pdb ]]} has multiple compounds: {bind(chainNames)}."))
          if (any(nchar(names(chainNames)) > 1)) {
            message(
              glue::glue('Chain names not in expection: {bind(names(chainNames), co = " | ")}')
            )
            return(file)
          }
          chains <- names(chainNames)
          newfiles <- add_filename_suffix(file, tolower(chains))
          .pymol_select_chains(file, chains, newfiles)
          gene_pdbChains <- nl(chainNames, tools::file_path_sans_ext(basename(newfiles)), FALSE)
          # match gene in chain name
          genesPattern <- make.names(names(used_pdbs))
          hasWhich <- vapply(
            genesPattern, function(p) any(grpl(names(gene_pdbChains), p)), logical(1)
          )
          if (any(hasWhich)) {
            genes <- names(used_pdbs)[hasWhich]
            isTheChain <- grpl(names(gene_pdbChains), genesPattern[hasWhich])
            message(
              glue::glue("Detected input gene ({bind(genes)}) in chain: {names(gene_pdbChains)[isTheChain]}.")
            )
            used_pdbs[ hasWhich ] <- gene_pdbChains[ isTheChain ]
            complexPdb <- pdb
            names(complexPdb) <- paste0(
              formal_name(names(gene_pdbChains)), collapse = "___"
            )
            used_pdbs <<- c(used_pdbs, complexPdb)
          }
          c(file, newfiles)
        } else {
          file
        }
      })
    lessUsedChain <- head(
      used_pdbs[ grpl(used_pdbs, "_[a-e]$") ], n = 2
    )
    if (length(lessUsedChain)) {
      x <- methodAdd(x, "对于复合体 PDB (文件中包含支链分子信息)，将使用对应的支链 (以 `pymol` 获取支链) 进行分子对接 (例如 {bind(names(lessUsedChain))}，使用 {bind(lessUsedChain)}) 。")
    }
    x$used_pdbs <- used_pdbs
    x$dock_layout <- lapply(x$dock_layout,
      function(vec) {
        vec[ names(vec) %in% names(used_pdbs) ] <- used_pdbs[ match(names(vec), names(used_pdbs)) ]
        vec
      })
    pdb_files_MultiChains <- unlist(
      pdb_files_MultiChains, use.names = FALSE
    )
    names(pdb_files_MultiChains) <- tools::file_path_sans_ext(basename(pdb_files_MultiChains))
    x$pdb_files_MultiChains <- pdb_files_MultiChains
    pdb.files <- c(pdb.files[!names(pdb.files) %in% pdb_MultiChains], pdb_files_MultiChains)
  }
  if (!is.null(pattern)) {
    pdb.files <- filter_pdbs(pdb.files, pattern)
  }
  namel(x, pdb.files)
}

get_pdb_from_alphaFold <- function(symbols, dir = "protein_pdb")
{
  message("Try get from alphaFold database (taxId: 9606, Human)")
  dir.create(dir, FALSE)
  info <- UniProt.ws::mapUniProt(
    from = "Gene_Name",
    to = "UniProtKB-Swiss-Prot",
    columns = c("accession", "id"),
    query = list(taxId = 9606, ids = symbols)
  )
  if (any(which <- !symbols %in% info$From)) {
    message(glue::glue("Not Got in alphaFold: {bind(symbols[which])}"))
  }
  if (any(duplicated(info$From))) {
    message('any(duplicated(info$From)), distinct herein.')
    info <- dplyr::distinct(info, From, .keep_all = TRUE)
  }
  files <- apply(info, 1, 
    function(vec) {
      id <- vec[[ "Entry" ]]
      if (!is.na(id)) {
        url <- glue::glue("https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v4.pdb")
        save <- file.path(dir, paste0(id, ".pdb"))
        if (!file.exists(save) || (file.exists(save) && !file.size(save))) {
          # res <- try(download.file(url, save))
          # if (inherits(res, "try-error")) {
          #   message(glue::glue("Download Failed of {id}, skip."))
          #   return(NULL)
          # }
          return(NULL)
        }
        nl(vec[[ "Entry" ]], save, FALSE)
      } else {
        NULL
      }
    }, simplify = FALSE)
  files <- unlist(files)
  info <- dplyr::filter(info, Entry %in% names(files))
  lst <- list(
    files = files, 
    info = info, used_pdbs = nl(info$From, info$Entry, FALSE)
  )
  lst
}

setMethod("filter", signature = c(x = "job_vina"),
  function(x, cpd, symbol, cid = NULL)
  {
    message("Custom specified docking.")
    if (x@step != 3L) {
      stop("x@step != 3L")
    }
    if (is.null(cid)) {
      cid <- unname(object(x)$cids[ match(cpd, names(object(x)$cids)) ])
    }
    if (is.null(x$used_pdbs)) {
      stop('is.null(x$used_pdbs).')
    }
    layout <- pdb <- x$used_pdbs[ match(symbol, names(x$used_pdbs)) ]
    names(layout) <- unlist(cid)
    x$dock_layout <- as.list(layout)
    return(x)
  })

setMethod("step4", signature = c(x = "job_vina"),
  function(x, time = 3600 * 2, savedir = paste0(x@sig, "_vina_space"),
    log = "~/vina.log", save.object = "vn3.rds", scoring = c("vina", "ad4"),
    exhaustiveness = 32,
    path_autodock_scripts = getOption("path_autodock_scripts"), ...)
  {
    step_message("Run vina ...")
    if (!nchar(Sys.which("prepare_gpf.py")) && !is.null(path_autodock_scripts)) {
      message(glue::glue("Add path to system PATH: {path_autodock_scripts}"))
      path <- normalizePath(path_autodock_scripts)
      Sys.setenv(PATH = paste(path, Sys.getenv("PATH"), sep = ":"))
    }
    runs <- tibble::tibble(
      Ligand = rep(names(x$dock_layout), lengths(x$dock_layout)),
      Receptor = tolower(unlist(x$dock_layout, use.names = FALSE))
    )
    x$show_layout <- runs
    runs <- apply(runs, 1, unname, simplify = FALSE)
    x$runs <- runs
    x$savedir <- savedir
    n <- 0
    saveRDS(x, save.object)
    scoring <- match.arg(scoring)
    message(glue::glue("Save this 'job_vina' as: {save.object}."))
    if (scoring == "ad4") {
      x <- methodAdd(x, "以 `AutoDock-Vina` 提供的工具 (`prepare_gpf.py`) (<https://github.com/ccsb-scripps/AutoDock-Vina>) 创建 GPF (grid parameter file)。")
      x <- methodAdd(x, "以 `ADFR` {cite_show('AutogridfrImpZhang2019')} 工具 `autogrid4` 计算亲和图谱 (Affinity Maps)。")
    }
    # x <- snapAdd(x, "以 `ADFR` 创建 Affinity Maps (详见方法章节) 。")
    if (any(grpl(names(x$res.receptor), "[A-Z]"))) {
      message('any(grpl(names(x$res.receptor), "[A-Z]")), convert to lower case.')
      names(x$res.receptor) <- tolower(names(x$res.receptor) )
    }
    if (is.remote(x)) {
      dir.create(savedir, FALSE)
      .script <- file.path(savedir, "script_remote.sh")
      scriptPrefix <- scriptPrefix(x)
      if (is.null(scriptPrefix)) {
        cat("", file = .script)
      } else {
        writeLines(scriptPrefix, .script)
      }
    } else {
      .script <- NULL
    }
    res <- pbapply::pblapply(x$runs,
      function(v) {
        lig <- x$res.ligand[[ v[1] ]]
        recep <- x$res.receptor[[ tolower(v[2]) ]]
        n <<- n + 1
        if (!is.null(lig) & !is.null(recep)) {
          vina_limit(
            lig, recep, time, dir = savedir, x = x, 
            .script = .script, stout = log, scoring = scoring, 
            exhaustiveness, ...
          )
        } else {
          NULL
        }
      }
    )
    if (is.remote(x)) {
      x$remote_operation <- dplyr::bind_rows(res)
    }
    x <- methodAdd(
      x, "运行 AutoDock-Vina {cite_show('AutodockVina1Eberha2021')} (参数，计分方式 score 设定为 {scoring}; 穷尽性 exhaustiveness 设定为 {exhaustiveness})。"
    )
    # x <- snapAdd(x, "以 `Autodock-Vina` 进行自动分子对接。")
    if (.Platform$OS.type == "unix" && Sys.which("notify-send") != "") {
      system("notify-send 'AutoDock vina' 'All job complete'")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_vina"),
  function(x, compounds, by.y, axis = "Ingredient_name", excludes = NULL, top = NULL,
    cutoff.af = NULL, sig.af = -5, maxShow = 20)
  {
    step_message("Summary and visualization for results.")
    x$summary_vina <- summary_vina(x$savedir)
    x$summary_vina <- dplyr::filter(x$summary_vina,
      PubChem_id %in% !!names(x$dock_layout),
      PDB_ID %in% tolower(unlist(x$dock_layout, use.names = FALSE))
    )
    res_dock <- dplyr::mutate(
      x$summary_vina, PubChem_id = as.integer(PubChem_id),
      hgnc_symbol = names(x$used_pdbs)[ match(PDB_ID, tolower(x$used_pdbs)) ]
    )
    if (!is.null(x$from_job_herb)) {
      res_dock <- map(res_dock, "PubChem_id",
        x$compounds, "PubChem_id", "Ingredient_name", col = "Ingredient_name")
      axis <- "Ingredient_name"
    } else {
      if (missing(compounds)) {
        compounds <- data.frame(PubChem_id = as.integer(object(x)$cids))
        if (is.null(names(object(x)$cids))) {
          compounds$Ingredient_name <- paste0("CID:", object(x)$cids)
        } else {
          compounds$Ingredient_name <- names(object(x)$cids)
        }
        by.y <- "PubChem_id"
        axis <- "Ingredient_name"
      }
      res_dock <- map(res_dock, "PubChem_id",
        compounds, "PubChem_id", "Ingredient_name", col = "Ingredient_name")
    }
    res_dock <- dplyr::arrange(res_dock, Affinity)
    res_dock <- dplyr::filter(res_dock, !is.na(!!rlang::sym(axis)))
    data <- dplyr::distinct(res_dock, PubChem_id, hgnc_symbol, .keep_all = TRUE)
    if (!is.null(excludes)) {
      data <- dplyr::filter(data, !hgnc_symbol %in% !!excludes)
    }
    trunc <- FALSE
    if (nrow(data) > maxShow) {
      trunc <- TRUE
      allProtein <- unique(data$hgnc_symbol)
      nProtein <- length(allProtein)
      each <- maxShow %/% nProtein
      data <- split_lapply_rbind(data, ~ PDB_ID, head, n = each)
    }
    data <- dplyr::mutate(
      data, receptor = paste0(hgnc_symbol, " (", s(PDB_ID, "-f1-model.*", ""), ")"),
      label = stringr::str_trunc(!!rlang::sym(axis), 30),
      label = paste0(label, " (CID:", PubChem_id, ")")
    )
    p.res_vina <- ggplot(data) + 
      geom_col(
        aes(x = reorder(label, Affinity, decreasing = TRUE),
          y = Affinity, fill = Affinity), width = .7
        ) +
      geom_text(data = dplyr::filter(data, Affinity <= 0),
        aes(x = label, y = Affinity - .5, label = round(Affinity, 1)), hjust = 1) +
      labs(x = "", y = "Affinity (kcal/mol)") +
      coord_flip() +
      ylim(zoRange(c(-1, data$Affinity, 1), 1.4)) +
      facet_wrap(~ receptor, ncol = 1, scales = "free_y") +
      theme()
    p.res_vina <- wrap_scale(p.res_vina, 15, nrow(data), h.size = .2)
    p.res_vina <- set_lab_legend(
      p.res_vina,
      glue::glue("{x@sig} Overall combining Affinity"),
      glue::glue("分子对接亲和度|||分子对接能量越低，代表亲和度越高 (图中对接的配体已注释化合物名称、PubChem ID，对接的受体已注释对应的基因名，以及对应的 PDB ID 或者 AlphaFold 数据库 ID)。")
    )
    res_dock <- set_lab_legend(
      res_dock,
      glue::glue("{x@sig} All combining Affinity data"),
      glue::glue("分子对接得分 (亲和度) 附表。")
    )
    t.sigData <- dplyr::filter(res_dock, Affinity < sig.af)
    t.showData <- dplyr::select(
      data, Compound = Ingredient_name, Protein = hgnc_symbol, 
      Affinity, Compound_CID = PubChem_id, Protein_Structure_ID = PDB_ID
    )
    x$dataInOverall <- data
    t.showData <- set_lab_legend(
      t.showData,
      glue::glue("{x@sig} Overall combining Affinity data"),
      glue::glue("分子对接亲和度概览。")
    )
    x <- snapAdd(x, "一共进行了 {nrow(res_dock)} 次对接。其中，有 {nrow(t.sigData)} 对配体受体组合的结合能 &lt; {sig.af} kcal/mol (具有良好的亲和力)。")
    if (trunc) {
      snap_each <- .stat_target_in_table(
        t.sigData, allProtein, "hgnc_symbol"
      )
      x <- snapAdd(x, "其中，{snap_each}。")
      x <- snapAdd(x, "如图{aref(p.res_vina)}展示了各个蛋白对应的 Top {each} 的对接亲和能。")
      snap_table <- glue::glue("表格{aref(t.showData)}为对应的蛋白与分子数据 (其余数据可见文件夹中数据表)。")
    } else {
      x <- snapAdd(x, "如图{aref(p.res_vina)}展示了蛋白与分子的对接亲和能。")
      snap_table <- glue::glue("表格{aref(t.showData)}为对应的蛋白与分子数据。")
    }
    x <- snapAdd(x, "{snap_table}(表格中，如蛋白结构 ID 以 AF 开头，则该蛋白结构获取于 AlphaFold 数据库，ID 为对应的数据库 ID) 。")
    x <- tablesAdd(
      x, t.sigData, t.showData, res_dock, unique_tops = data
    )
    x <- plotsAdd(x, p.res_vina)
    return(x)
  })

.stat_target_in_table <- function(data, targets, col) {
  snaps <- vapply(targets, FUN.VALUE = character(1), 
    function(x) {
      data <- dplyr::filter(data, !!rlang::sym(col) == !!x)
      glue::glue("{x} 包含 {nrow(data)} 例")
    })
  bind(snaps)
}

setMethod("step6", signature = c(x = "job_vina"),
  function(x, time = 30, top = 1, save = TRUE, unique = FALSE, 
    symbol = NULL, cpd = NULL, rerun = FALSE)
  {
    step_message("Use pymol for all visualization.")
    data <- x$dataInOverall
    if (!is.null(symbol)) {
      data <- dplyr::filter(data, hgnc_symbol %in% !!symbol)
    }
    if (!is.null(cpd)) {
      data <- dplyr::filter(data, Ingredient_name %in% !!cpd)
    }
    if (unique) {
      data <- dplyr::distinct(data, hgnc_symbol, .keep_all = TRUE)
    }
    if (is.null(top) && nrow(data) > 20) {
      message(
        glue::glue("Too many data for manually drawingl, switch to 10.")
      )
      top <- 10L
    }
    if (!is.null(top)) {
      data <- split_lapply_rbind(data, ~ hgnc_symbol, head, n = top)
    }
    fun_draw <- function(data) {
      figs <- pbapply::pbapply(data, 1, simplify = FALSE,
        function(v) {
          vinaShow(
            v[[ "Combn" ]], v[[ "PDB_ID" ]], timeLimit = time, save = save, dir = x$savedir
          )
        }
      )
    }
    figs <- expect_local_data(
      "tmp", "vina_pymol", fun_draw, list(data = data), rerun = rerun
    )
    figs <- unlist(figs, recursive = FALSE)
    figs <- lapply(figs, 
      function(x) {
        if (is(x, "file_fig")) {
          as_data_binary(x)
        } else {
          x
        }
      })
    names(figs) <- paste0(
      "Top", "_", data$hgnc_symbol, "_", data$PubChem_id
    )
    figs <- set_lab_legend(
      figs,
      glue::glue("{x@sig} {data$hgnc_symbol} {data$PubChem_id} docking visualization"),
      glue::glue("Top 亲和度分子对接结果|||蛋白(Symbol: {data$hgnc_symbol}) (Protein Structure ID: {data$PDB_ID}) 与化合物 (PubChem CID: {data$PubChem_id}) (name: {data$Ingredient_name})，亲和度为 {data$Affinity} kcal/mol。")
    )
    # x <- snapAdd(x, "使用 `pymol` 将分子对接结果可视化。")
    x$data_selectVis <- data
    x <- plotsAdd(x, figs)
    return(x)
  })

setMethod("step7", signature = c(x = "job_vina"),
  function(x, save = TRUE, rerun = FALSE){
    step_message("Show docking results in deep and in detail.")
    data <- x$data_selectVis
    fun_draw <- function(data) {
      figs <- pbapply::pbapply(data, 1, simplify = FALSE,
        function(v) {
          vinaShow(
            v[[ "Combn" ]], v[[ "PDB_ID" ]], save = save, detail = TRUE, dir = x$savedir
          )
        }
      )
    }
    figs <- expect_local_data(
      "tmp", "vina_pymol_detail", fun_draw, list(data = data), rerun = rerun
    )
    figs <- unlist(figs, recursive = FALSE)
    figs <- lapply(figs, 
      function(x) {
        if (is(x, "file_fig")) {
          as_data_binary(x)
        } else {
          x
        }
      })
    names(figs) <- paste0(
      "Top_", data$hgnc_symbol, "_", data$PubChem_id
    )
    figs <- set_lab_legend(
      figs,
      glue::glue("{x@sig} {data$hgnc_symbol} {data$PubChem_id} docking interaction details"),
      glue::glue("Top 亲和度分子对接局部细节|||蛋白(Symbol: {data$hgnc_symbol}) (Protein Structure ID: {data$PDB_ID}) 与化合物 (PubChem CID: {data$PubChem_id}) (name: {data$Ingredient_name})，亲和度为 {data$Affinity} kcal/mol。图中蛋白与分子之间的虚线为可能存在的氢键结合。")
    )
    x <- plotsAdd(x, figs)
    return(x)
  })

setMethod("step8", signature = c(x = "job_vina"),
  function(x){
    step_message("Merge as pdb for molecular dynamics simulation")
    x@tables$step5$res_dock
    x$res_dock_merge <- .merge_ligand_recepter_as_pdb(x@tables$step5$res_dock, x$savedir)
    return(x)
  })

.merge_ligand_recepter_as_pdb <- function(res_dock, dir, overwrite = FALSE) {
  if (!nchar(Sys.which("obabel"))) {
    stop('!nchar(Sys.which("obabel")).')
  }
  fun_convert <- function(file) {
    real <- tools::file_path_sans_ext(file)
    file_new <- paste0(real, ".pdb")
    if (!file.exists(file_new) || overwrite) {
      system(glue::glue("obabel {file} -O {file_new}"))
    }
    return(file_new)
  }
  res_dock$pdb_merge <- apply(res_dock, 1,
    function(v) {
      subdir <- Combn <- v[[ "Combn" ]]
      recep <- paste0(v[[ "PDB_ID" ]], ".pdbqt")
      path <- file.path(dir, subdir)
      maybethat <- list.files(
        path, recep, ignore.case = TRUE, full.names = TRUE
      )
      if (length(maybethat)) {
        file_recep <- maybethat[1]
      }
      file_ligand <- v[[ "file" ]]
      files <- lapply(c(file_recep, file_ligand),
        function(file) {
          message(glue::glue("Convert: {file}"))
          fun_convert(file)
        })
      cmd_load <- glue::glue(
        "load {files[[1]]}, recep; load {files[[2]]}, ligand; create complex, recep or ligand;"
      )
      file_complex <- file.path(path, paste0(Combn, ".pdb"))
      cmd <- glue::glue("{pg('pymol')} -c -Q -d '{cmd_load} save {file_complex}, complex; quit'")
      system(cmd)
      return(file_complex)
    })
  res_dock
}

setMethod("upload", signature = c(x = "job_vina"),
  function(x, ..., testFinish = TRUE){
    if (!is.remote(x)) {
      stop('!is.remote(x).')
    }
    if (x@step < 4L) {
      stop('x@step < 4L.')
    }
    if (!is.null(x$.upload) && x$.upload) {
      stop('!is.null(x$.upload). The job has been upload?')
    }
    if (is.null(x$remote_operation)) {
      stop('is.null(x$remote_operation).')
    }
    data <- x$remote_operation
    script <- unique(data$script)
    if (length(script) != 1) {
      stop('length(script) != 1.')
    }
    .run_script_in_remote(script, x$wd, remote = x$remote, ...)
    if (testFinish) {
      expectFile <- paste0(
        tail(data$remote_wd, n = 1), "/", tail(data$output, n = 1)
      )
      testRem_file.exists(x, expectFile, x$wait)
    }
    x$.upload <- TRUE
    return(x)
  })

setMethod("pull", signature = c(x = "job_vina"),
  function(x, force = FALSE){
    if (!is.remote(x)) {
      stop('!is.remote(x).')
    }
    if (x@step < 4L) {
      stop('x@step < 4L.')
    }
    if (is.null(x$.upload) || !x$.upload) {
      stop('is.null(x$.upload) || !x$.upload.')
    }
    if (is.null(x$remote_operation)) {
      stop('is.null(x$remote_operation).')
    }
    data <- x$remote_operation
    continue <- TRUE
    gotThat <- pbapply::pbapply(data, 1, simplify = FALSE,
      function(args) {
        if ((notfirst <- !file.exists(args[[ "to" ]])) && continue) {
          try(get_file_from_remote(
            args[[ "output" ]], args[[ "remote_wd" ]], args[[ "to" ]], remote = x$remote
          ))
        }
        second <- file.exists(args[[ "to" ]])
        if (notfirst && !second && !force) {
          continue <<- FALSE
        }
        second
      })
    gotThat <- unlist(gotThat)
    message(glue::glue("Got results: {length(which(gotThat))} / {length(gotThat)}"))
    return(x)
  })

pretty_docking <- function(protein, ligand, path,
  save = "annotation.png",
  script = paste0(.expath, "/pretty_docking.pymol"), detail = TRUE)
{
  script <- readLines(script)
  script <- gs(script, "{{protein.pdbqt}}", protein, fixed = TRUE)
  script <- gs(script, "{{ligand.pdbqt}}", ligand, fixed = TRUE)
  if (detail) {
    # script <- c(script, "center Ligand", "zoom Ligand, 4", "clip slab, 10")
  } else {
    script <- c(
      script, "center Protein", "zoom Protein, 4", "", "label all, ''"
    )
  }
  temp <- tempfile("Pymol", fileext = ".pml")
  writeLines(script, temp)
  cli::cli_alert_info(paste0("Pymol run script: ", temp))
  output <- file.path(path, save)
  message("Save png to ", output)
  expr <- paste0(" png ", save, ",2500,2000,dpi=300")
  gett(expr)
  cdRun(pg("pymol"), " ",
    " -d \"run ", temp, "\"",
    " -d \"ray; ", expr, "\"",
    path = path)
  return(output)
}

vina_limit <- function(lig, recep, timeLimit = 120, dir = "vina_space", ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit, dir = dir), TRUE)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoring = c("vina", "ad4"),
  stout = "~/vina.log", timeLimit = 60,
  excludes.atom = c("G0", "CG0"), .script, x, remoteTest = FALSE)
{
  remote <- FALSE
  if (!missing(x)) {
    if (is.remote(x)) {
      remote <- TRUE
    }
  }
  if (!file.exists(dir)) {
    dir.create(dir, FALSE)
  }
  scoring <- match.arg(scoring)
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
  if (!file.exists(paste0(wd, "/", subdir, "_out.pdbqt"))) {
    dir.create(wd, FALSE)
    if (!remoteTest) {
      file.copy(c(recep, lig), wd)
    }
    .message_info("Generating affinity maps", subdir)
    .cdRun <- function(...) {
      if (!remoteTest) {
        cdRun(..., path = wd)
      }
    }
    files <- basename(c(lig, recep))
    if (scoring == "ad4") {
      .cdRun(pg("prepare_gpf.py"), " -l ", files[1], " -r ", files[2], " -y")
      if (!is.null(excludes.atom)) {
        message("Excludes atom type: ", paste0(excludes.atom, collapse = ", "))
        ## ligand type and map file
        .cdRun("sed -i",
          " -e '/", paste0(paste0("^map.*", excludes.atom), collapse = "\\|"), "/d'",
          " -e 's/", paste0(paste0("\\b", excludes.atom, "\\b"), collapse = "\\|"), "//g'",
          " ", reals[2], ".gpf")
      }
      .cdRun(pg("autogrid4"), " -p ", reals[2], ".gpf ", " -l ", reals[2], ".glg")
      config <- paste0(" --maps ", reals[2])
    } else if (scoring == "vina") {
      box <- glue::glue("center_x = 0\ncenter_y = 0\ncenter_z = 0\nsize_x = 20.0\nsize_y = 20.0\nsize_z = 20.0")
      writeLines(box, file.path(wd, "receptor.box.txt"))
      config <- glue::glue(" --receptor {files[2]} --config receptor.box.txt ")
    }
    if (remote) {
      message("Run in remote server.")
      if (!remoteTest) {
        cdRun("scp -r ", subdir, " ", x$remote, ":", x$wd, path = dir)
      }
      remoteWD <- x$wd
      x$wd <- paste0(x$wd, "/", subdir)
      output <- paste0(subdir, "_out.pdbqt")
      rem_run("timeout ", timeLimit, " ",
        pg("vina", TRUE), " --ligand ", files[1],
        " --scoring ", scoring,
        config,
        " --exhaustiveness ", exhaustiveness,
        " --out ", output, 
        " >> ", stout, .script = .script, .append = TRUE)
      return(
        namel(script = .script, output, remote_wd = x$wd, to = file.path(wd, output))
      )
    } else {
      cat("\n$$$$\n", date(), "\n", subdir, "\n\n", file = stout, append = TRUE)
      try(.cdRun("timeout ", timeLimit, 
          " ", pg("vina"), "  --ligand ", files[1],
          " --scoring ", scoring,
          config,
          " --exhaustiveness ", exhaustiveness,
          " --out ", subdir, "_out.pdbqt",
          " >> ", stout), TRUE)
    }
  }
}

vinaShow <- function(Combn, recep, subdir = Combn, dir = "vina_space",
  timeLimit = 3, backup = NULL, save = TRUE, detail = FALSE)
{
  if (!file.exists(path <- file.path(dir, subdir))) {
    stop('!file.exists(path <- file.path(dir, subdir)).')
  } 
  wd <- path
  out <- paste0(Combn, "_out.pdbqt")
  recep <- paste0(recep, ".pdbqt")
  if (!file.exists(recep)) {
    maybethat <- list.files(wd, recep, ignore.case = TRUE)
    if (length(maybethat)) {
      recep <- maybethat[1]
    }
  }
  res <- paste0(Combn, ".png")
  if (detail) {
    res <- paste0("detail_", res)
  }
  .cdRun <- function(...) cdRun(..., path = wd)
  img <- file.path(wd, res)
  if (file.exists(img)) {
    file.remove(img)
  }
  if (TRUE) {
    pretty_docking(recep, out, wd, save = res, detail = detail)
  } else {
    expr <- paste0(" png ", res, ",2500,2000,dpi=300")
    gett(expr)
    if (!save) {
      expr <- ""
    }
    try(.cdRun("timeout ", timeLimit, 
        " ", pg("pymol"), " ",
        " -d \"load ", out, ";",
        " load ", recep, ";",
        " ray; zoom; bg white; color grey70, Protein; ", expr, "\" "), TRUE)
  }
  if (is.character(backup)) {
    dir.create(backup, FALSE)
    file.copy(img, backup, TRUE)
  }
  fig <- as_data_binary(.file_fig(img))
  lab(fig) <- paste0("docking ", Combn, if (!detail) NULL else " detail")
  return(fig)
}

summary_vina <- function(space = "vina_space", pattern = "_out\\.pdbqt$")
{
  files <- list.files(space, pattern, recursive = TRUE, full.names = TRUE)
  res_dock <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (length(lines) >= 1) {
        name <- gsub("_out\\.pdbqt", "", basename(file))
        top <- stringr::str_extract(lines[2], "[\\-0-9.]{1,}")
        top <- as.double(top)
        names(top) <- name
        top
      }
    })
  res_dock <- unlist(res_dock)
  res_dock <- tibble::tibble(
    Combn = names(res_dock), Affinity = unname(res_dock)
  )
  res_dock <- dplyr::mutate(
    res_dock, PubChem_id = gs(Combn, "(.*)_into_.*", "\\1"),
    PDB_ID = tolower(gs(Combn, ".*_into_(.*)", "\\1")),
    dir = paste0(space, "/", Combn),
    file = paste0(dir, "/", Combn, "_out.pdbqt")
  )
  dplyr::select(res_dock, PubChem_id, PDB_ID, Affinity, dir, file, Combn)
}

smiles_as_sdfs.obabel <- function(smiles) {
  lst.sdf <- pbapply::pbsapply(smiles, simplify = FALSE,
    function(smi) {
      ChemmineOB::convertFormat("SMI", "SDF", source = test)
    })
  lst.sdf
}

tbmerge <- function(x, y, ...) {
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)
  tibble::as_tibble(data.table::merge.data.table(x, y, ...))
}

ld_cols <- function(file, sep = "\t") {
  line <- readLines(file, n = 1)
  strsplit(line, sep)[[ 1 ]]
}

lst_clear0 <- function(lst, len = 0) {
  lst[ vapply(lst, function(v) if (length(v) > len) TRUE else FALSE, logical(1)) ]
}

ld_cutRead <- function(file, cols, abnum = TRUE, sep = "\t", tmp = "/tmp/ldtmp.txt") {
  if (is.character(cols)) {
    names <- ld_cols(file, sep)
    message("Find columns of: ", paste0(names[ names %in% cols ], collapse = ", "))
    pos <- which( names %in% cols )
    if (abnum) {
      pos <- head(pos, n = length(cols))
    }
  } else {
    pos <- cols
  }
  if (!file.exists(tmp)) {
    if (tools::file_ext(file) == "gz") {
      cdRun(glue::glue("zcat {file} | cut -f {paste0(pos, collapse = ',')} > {tmp}"))
    } else {
      cdRun(glue::glue("cut -f {paste0(pos, collapse = ',')} {file} > {tmp}"))
    }
  }
  ftibble(tmp)
}

mk_prepare_ligand.sdf <- function(sdf_file, mkdir.pdbqt = "pdbqt",
  check = FALSE, per_molecule = TRUE)
{
  dir.create(mkdir.pdbqt, FALSE, recursive = TRUE)
  check_sdf_validity <- function(file) {
    lst <- sep_list(readLines(file), "^\\${4,}$")
    lst <- lst[ - length(lst) ]
    osum <- length(lst)
    lst <- lapply(lst,
      function(line) {
        pos <- grep("^\\s*-OEChem|M\\s*END", line)
        if (length(pos) >= 2L && pos[2] - pos[1] > 5) {
          line
        }
      })
    lst <- lst[ !vapply(lst, is.null, logical(1)) ]
    nsum <- length(lst)
    list(data = lst, osum = osum, nsum = nsum, dec = osum - nsum)
  }
  if (check) {
    lst <- check_sdf_validity(sdf_file)
    lines <- unlist(lst$data)
    lst <- lst[ -1 ]
    writeLines(lines, sdf_file <- gsub("\\.sdf", "_modified.sdf", sdf_file))
  } else {
    lst <- list(file = sdf_file)
  }
  if (per_molecule) {
    records <- vinaFuns$split_sdf_records(sdf_file)
    split_dir <- file.path(
      dirname(sdf_file),
      paste0(tools::file_path_sans_ext(basename(sdf_file)), "_split_for_meeko")
    )
    dir.create(split_dir, FALSE, recursive = TRUE)
    status <- lapply(seq_along(records),
      function(i) {
        cid <- vinaFuns$sdf_record_id(records[[ i ]], fallback = as.character(i))
        cid_file <- gsub("[^A-Za-z0-9_.-]", "_", cid)
        one_sdf <- file.path(split_dir, paste0(cid_file, ".sdf"))
        log_file <- file.path(split_dir, paste0(cid_file, ".meeko.log"))
        writeLines(records[[ i ]], one_sdf)
        cmd <- paste(
          pg("mk_prepare_ligand.py"),
          "-i", shQuote(one_sdf),
          "--multimol_outdir", shQuote(mkdir.pdbqt)
        )
        res <- try(cdRun(cmd), TRUE)
        data.frame(
          CID = cid,
          input_sdf = one_sdf,
          meeko_ok = !inherits(res, "try-error"),
          log_file = log_file,
          stringsAsFactors = FALSE
        )
      })
    lst$meeko_status <- if (length(status)) {
      tibble::as_tibble(do.call(rbind, status))
    } else {
      tibble::tibble()
    }
  } else {
    cdRun(pg("mk_prepare_ligand.py"), " -i ", sdf_file,
      " --multimol_outdir ", mkdir.pdbqt)
    lst$meeko_status <- tibble::tibble()
  }
  got <- vinaFuns$collect_pdbqt(mkdir.pdbqt)
  lst$file <- sdf_file
  lst$pdbqt <- got$pdbqt
  lst$pdbqt.num <- length(lst$pdbqt)
  lst$pdbqt.cid <- got$pdbqt.cid
  return(lst)
}

vinaFuns$split_sdf_records <- function(file)
{
  lines <- readLines(file, warn = FALSE)
  if (!length(lines)) {
    return(list())
  }
  end <- grep("^\\$\\$\\$\\$", lines)
  if (!length(end)) {
    return(list(lines))
  }
  start <- c(1L, end[-length(end)] + 1L)
  records <- Map(function(i, j) lines[i:j], start, end)
  records <- records[lengths(records) > 1L]
  records
}

vinaFuns$sdf_record_id <- function(record, fallback = NA_character_)
{
  if (!length(record)) {
    return(fallback)
  }
  title <- trimws(record[[ 1L ]])
  if (grepl("^[0-9]+$", title)) {
    return(title)
  }
  pos <- grep("^> *<PUBCHEM_COMPOUND_CID>", record)
  if (length(pos) && length(record) >= pos[[ 1L ]] + 1L) {
    cid <- trimws(record[[ pos[[ 1L ]] + 1L ]])
    if (grepl("^[0-9]+$", cid)) {
      return(cid)
    }
  }
  cid <- stringr::str_extract(title, "[0-9]+")
  if (!is.na(cid)) {
    return(cid)
  }
  fallback
}

vinaFuns$sdf_ids <- function(file)
{
  records <- vinaFuns$split_sdf_records(file)
  ids <- vapply(records, vinaFuns$sdf_record_id, character(1), fallback = NA_character_)
  ids[!is.na(ids)]
}

vinaFuns$filter_sdf_records <- function(file, ids, output)
{
  ids <- as.character(ids)
  records <- vinaFuns$split_sdf_records(file)
  keep <- vapply(records,
    function(record) {
      vinaFuns$sdf_record_id(record) %in% ids
    }, logical(1))
  records <- records[keep]
  if (!length(records)) {
    return(NULL)
  }
  writeLines(unlist(records, use.names = FALSE), output)
  output
}

vinaFuns$sdf_file_has_records <- function(file)
{
  file.exists(file) && file.size(file) > 0L &&
    any(grepl("^\\$\\$\\$\\$", readLines(file, warn = FALSE)))
}

vinaFuns$collect_pdbqt <- function(mkdir.pdbqt)
{
  files <- list.files(mkdir.pdbqt, "\\.pdbqt$", full.names = TRUE)
  cid <- tools::file_path_sans_ext(basename(files))
  cid <- stringr::str_extract(cid, "[0-9]+")
  keep <- !is.na(cid)
  list(pdbqt = files[keep], pdbqt.cid = cid[keep])
}

vinaFuns$query_pubchem_3d_sdfs <- function(cids, dir_save = "pubchem_3d_sdf",
  filename = "pubchem_3d.sdf", record_type = c("3d", "2d"),
  overwrite = FALSE, download_method = "libcurl")
{
  record_type <- match.arg(record_type)
  cids <- unique(as.character(cids))
  dir.create(dir_save, FALSE, recursive = TRUE)
  output <- file.path(dir_save, filename)
  empty_status <- tibble::tibble(
    CID = character(), record_type = character(),
    sdf_file = character(), got_sdf = logical(), from_cache = logical(),
    reason = character()
  )
  if (!length(cids)) {
    return(list(file = NULL, status = empty_status))
  }
  if (!overwrite && vinaFuns$sdf_file_has_records(output)) {
    ids_output <- unique(as.character(vinaFuns$sdf_ids(output)))
    if (all(cids %in% ids_output)) {
      message("Use existing PubChem ", record_type, " SDF: ", output)
      status <- tibble::tibble(
        CID = cids, record_type = record_type,
        sdf_file = output, got_sdf = TRUE, from_cache = TRUE,
        reason = "combined_sdf_cache"
      )
      return(list(file = output, status = status))
    }
  }
  files <- character()
  status <- lapply(cids,
    function(cid) {
      file <- file.path(dir_save, paste0(cid, "_", record_type, ".sdf"))
      url <- glue::glue(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type={record_type}"
      )
      ok <- FALSE
      from_cache <- FALSE
      reason <- ""
      if (!overwrite && file.exists(file) && file.size(file) > 0L) {
        ok <- vinaFuns$sdf_file_has_records(file)
        from_cache <- ok
        reason <- if (ok) "cid_sdf_cache" else "invalid_cached_sdf"
      }
      if (!ok) {
        message("Download PubChem ", record_type, " SDF: CID ", cid)
        res <- try(utils::download.file(
          url, file, quiet = FALSE, mode = "wb", method = download_method
        ), TRUE)
        ok <- !inherits(res, "try-error") && vinaFuns$sdf_file_has_records(file)
        reason <- if (ok) "downloaded" else "download_or_sdf_validation_failed"
      }
      if (ok) {
        files <<- c(files, file)
      } else {
        warning("Failed to obtain PubChem ", record_type, " SDF for CID ", cid,
          call. = FALSE)
      }
      data.frame(
        CID = as.character(cid), record_type = record_type,
        sdf_file = if (ok) file else NA_character_,
        got_sdf = ok, from_cache = from_cache, reason = reason,
        stringsAsFactors = FALSE
      )
    })
  status <- tibble::as_tibble(do.call(rbind, status))
  if (length(files)) {
    writeLines(unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE), output)
  } else {
    output <- NULL
  }
  list(file = output, status = status)
}

vinaFuns$standardize_ligand_sdf_rdkit <- function(sdf_file, output,
  strip_salts = TRUE, neutralize_ligand = FALSE, ligand_minimize = TRUE,
  ligand_forcefield = c("MMFF", "UFF"), seed = 614L)
{
  ligand_forcefield <- match.arg(ligand_forcefield)
  py <- tempfile("vina_rdkit_ligand_", fileext = ".py")
  code <- c(
    "import sys, csv",
    "from rdkit import Chem",
    "from rdkit.Chem import AllChem",
    "try:",
    "    from rdkit.Chem.MolStandardize import rdMolStandardize",
    "except Exception:",
    "    rdMolStandardize = None",
    "infile, outfile = sys.argv[1], sys.argv[2]",
    "strip_salts = bool(int(sys.argv[3]))",
    "neutralize = bool(int(sys.argv[4]))",
    "minimize = bool(int(sys.argv[5]))",
    "forcefield = sys.argv[6]",
    "seed = int(sys.argv[7])",
    "suppl = Chem.SDMolSupplier(infile, sanitize=False, removeHs=False)",
    "writer = Chem.SDWriter(outfile)",
    "rows = []",
    "for i, mol in enumerate(suppl):",
    "    name = str(i + 1)",
    "    try:",
    "        if mol is None:",
    "            raise ValueError('RDKit returned None')",
    "        if mol.HasProp('_Name') and mol.GetProp('_Name').strip():",
    "            name = mol.GetProp('_Name').strip()",
    "        try:",
    "            Chem.SanitizeMol(mol)",
    "        except Exception:",
    "            mol.UpdatePropertyCache(strict=False)",
    "            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)",
    "        fragment_count = len(Chem.GetMolFrags(mol))",
    "        if strip_salts and fragment_count > 1 and rdMolStandardize is not None:",
    "            mol = rdMolStandardize.LargestFragmentChooser().choose(mol)",
    "        if neutralize and rdMolStandardize is not None:",
    "            mol = rdMolStandardize.Uncharger().uncharge(mol)",
    "        mol = Chem.AddHs(mol, addCoords=True)",
    "        params = AllChem.ETKDGv3()",
    "        params.randomSeed = seed + i",
    "        params.useRandomCoords = True",
    "        emb = AllChem.EmbedMolecule(mol, params)",
    "        if emb != 0:",
    "            emb = AllChem.EmbedMolecule(mol, randomSeed=seed + i, useRandomCoords=True)",
    "        if emb != 0:",
    "            raise ValueError('3D embedding failed')",
    "        if minimize:",
    "            try:",
    "                if forcefield.upper() == 'MMFF' and AllChem.MMFFHasAllMoleculeParams(mol):",
    "                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)",
    "                else:",
    "                    AllChem.UFFOptimizeMolecule(mol, maxIters=500)",
    "            except Exception:",
    "                pass",
    "        mol.SetProp('_Name', name)",
    "        writer.write(mol)",
    "        rows.append([name, 'TRUE', str(fragment_count), ''])",
    "    except Exception as e:",
    "        rows.append([name, 'FALSE', '', str(e)])",
    "writer.close()",
    "with open(outfile + '.status.tsv', 'w', newline='') as f:",
    "    w = csv.writer(f, delimiter='\\t')",
    "    w.writerow(['CID', 'rdkit_ok', 'fragment_count_before', 'reason'])",
    "    w.writerows(rows)"
  )
  writeLines(code, py)
  cmd <- paste(
    pg("docking_python"), shQuote(py), shQuote(sdf_file), shQuote(output),
    as.integer(strip_salts), as.integer(neutralize_ligand),
    as.integer(ligand_minimize), ligand_forcefield, as.integer(seed)
  )
  res <- try(cdRun(cmd), TRUE)
  status_file <- paste0(output, ".status.tsv")
  status <- if (file.exists(status_file)) {
    tibble::as_tibble(utils::read.delim(status_file, stringsAsFactors = FALSE))
  } else {
    tibble::tibble()
  }
  list(
    file = if (!inherits(res, "try-error") && file.exists(output) && file.size(output)) output else NULL,
    status = status,
    error = inherits(res, "try-error")
  )
}

vinaFuns$prepare_ligand_pdbqt_obabel <- function(sdf_file, mkdir.pdbqt = "pdbqt",
  ids = NULL, partialcharge = "gasteiger")
{
  dir.create(mkdir.pdbqt, FALSE, recursive = TRUE)
  records <- vinaFuns$split_sdf_records(sdf_file)
  if (!is.null(ids)) {
    ids <- as.character(ids)
    keep <- vapply(records,
      function(record) vinaFuns$sdf_record_id(record) %in% ids,
      logical(1)
    )
    records <- records[keep]
  }
  if (!length(records)) {
    return(list(
      pdbqt = character(), pdbqt.cid = character(),
      status = tibble::tibble(
        CID = character(), input_sdf = character(), pdbqt_file = character(),
        obabel_ok = logical(), reason = character()
      )
    ))
  }
  split_dir <- file.path(dirname(sdf_file), paste0(
    tools::file_path_sans_ext(basename(sdf_file)), "_split_for_obabel"
  ))
  dir.create(split_dir, FALSE, recursive = TRUE)
  status <- lapply(seq_along(records), function(i) {
    cid <- vinaFuns$sdf_record_id(records[[ i ]], fallback = as.character(i))
    cid_file <- gsub("[^A-Za-z0-9_.-]", "_", cid)
    one_sdf <- file.path(split_dir, paste0(cid_file, ".sdf"))
    out <- file.path(mkdir.pdbqt, paste0(cid_file, ".pdbqt"))
    writeLines(records[[ i ]], one_sdf)
    base_cmd <- paste(
      pg("obabel"),
      "-isdf", shQuote(one_sdf),
      "-opdbqt", "-O", shQuote(out),
      "-r", "-h"
    )
    cmd <- if (!is.null(partialcharge) && nchar(partialcharge)) {
      paste(base_cmd, "--partialcharge", shQuote(partialcharge))
    } else {
      base_cmd
    }
    res <- try(cdRun(cmd), TRUE)
    ok <- !inherits(res, "try-error") && file.exists(out) && file.size(out) > 0L
    reason <- if (ok) "" else as.character(res)[1]
    if (!ok && !is.null(partialcharge) && nchar(partialcharge)) {
      message("Open Babel retry without partial charge for CID ", cid)
      res2 <- try(cdRun(base_cmd), TRUE)
      ok <- !inherits(res2, "try-error") && file.exists(out) && file.size(out) > 0L
      reason <- if (ok) "partialcharge_failed_retry_without_partialcharge" else as.character(res2)[1]
    }
    data.frame(
      CID = as.character(cid), input_sdf = one_sdf,
      pdbqt_file = if (ok) out else NA_character_,
      obabel_ok = ok, reason = reason, stringsAsFactors = FALSE
    )
  })
  status <- tibble::as_tibble(do.call(rbind, status))
  status_ok <- status[status$obabel_ok %in% TRUE, , drop = FALSE]
  list(
    pdbqt = status_ok$pdbqt_file,
    pdbqt.cid = as.character(status_ok$CID),
    status = status
  )
}

vinaFuns$prepare_ligand_pdbqt_recovery <- function(cids, sdf_2d, sdf_3d = NULL,
  dir_save = "ligand", mkdir.pdbqt = "pdbqt", use_pubchem_3d = TRUE,
  use_obgen = TRUE, rdkit_fallback = TRUE, obabel_fallback = TRUE,
  strip_salts = TRUE, neutralize_ligand = FALSE,
  ligand_forcefield = c("MMFF", "UFF"), ligand_minimize = TRUE,
  ligand_per_molecule = TRUE, obabel_partialcharge = "gasteiger",
  overwrite_pubchem_3d = FALSE, overwrite_obgen = FALSE, cl = NULL)
{
  ligand_forcefield <- match.arg(ligand_forcefield)
  cids <- unique(as.character(cids))
  dir.create(dir_save, FALSE, recursive = TRUE)
  dir.create(mkdir.pdbqt, FALSE, recursive = TRUE)
  got_files <- character()
  got_cids <- character()
  lst_status <- list()
  lst_sdf <- list()

  register_pdbqt <- function(res, label) {
    if (is.null(res$pdbqt) || !length(res$pdbqt)) {
      return(invisible(NULL))
    }
    keep <- !as.character(res$pdbqt.cid) %in% got_cids
    if (any(keep)) {
      got_files <<- c(got_files, res$pdbqt[keep])
      got_cids <<- c(got_cids, as.character(res$pdbqt.cid[keep]))
    }
    invisible(NULL)
  }

  run_meeko <- function(label, sdf_file, ids, suffix = "meeko") {
    if (!length(ids)) {
      return(invisible(NULL))
    }
    filtered <- file.path(dir_save, paste0(label, "_", suffix, "_input.sdf"))
    filtered <- vinaFuns$filter_sdf_records(sdf_file, ids, filtered)
    if (is.null(filtered)) {
      return(invisible(NULL))
    }
    outdir <- file.path(mkdir.pdbqt, paste0(label, "_", suffix))
    res <- mk_prepare_ligand.sdf(filtered, outdir, per_molecule = ligand_per_molecule)
    if (!is.null(res$meeko_status) && nrow(res$meeko_status)) {
      lst_status[[ paste0(label, "_", suffix) ]] <<- dplyr::mutate(
        res$meeko_status,
        strategy = label, tool = "meeko", input_stage = suffix
      )
    }
    register_pdbqt(res, label)
    invisible(NULL)
  }

  run_rdkit_meeko <- function(label, sdf_file, ids) {
    if (!rdkit_fallback || !length(ids)) {
      return(invisible(NULL))
    }
    filtered <- file.path(dir_save, paste0(label, "_rdkit_input.sdf"))
    filtered <- vinaFuns$filter_sdf_records(sdf_file, ids, filtered)
    if (is.null(filtered)) {
      return(invisible(NULL))
    }
    rdkit_file <- file.path(dir_save, paste0(label, "_rdkit3D.sdf"))
    rdkit_res <- vinaFuns$standardize_ligand_sdf_rdkit(
      filtered, rdkit_file,
      strip_salts = strip_salts,
      neutralize_ligand = neutralize_ligand,
      ligand_minimize = ligand_minimize,
      ligand_forcefield = ligand_forcefield
    )
    if (nrow(rdkit_res$status)) {
      lst_status[[ paste0(label, "_rdkit") ]] <<- dplyr::mutate(
        rdkit_res$status,
        strategy = label, tool = "rdkit", input_stage = "standardize"
      )
    }
    if (!is.null(rdkit_res$file)) {
      lst_sdf[[ paste0(label, "_rdkit") ]] <<- rdkit_res$file
      run_meeko(label, rdkit_res$file, setdiff(ids, got_cids), suffix = "rdkit_meeko")
    }
    invisible(NULL)
  }

  run_obabel <- function(label, sdf_file, ids) {
    if (!obabel_fallback || !length(ids)) {
      return(invisible(NULL))
    }
    obabel_dir <- file.path(mkdir.pdbqt, paste0(label, "_obabel"))
    obabel_res <- vinaFuns$prepare_ligand_pdbqt_obabel(
      sdf_file,
      mkdir.pdbqt = obabel_dir,
      ids = ids,
      partialcharge = obabel_partialcharge
    )
    if (!is.null(obabel_res$status) && nrow(obabel_res$status)) {
      lst_status[[ paste0(label, "_obabel") ]] <<- dplyr::mutate(
        obabel_res$status,
        strategy = label, tool = "obabel", input_stage = "fallback"
      )
    }
    register_pdbqt(obabel_res, label)
    invisible(NULL)
  }

  run_source <- function(label, sdf_file) {
    if (is.null(sdf_file) || !file.exists(sdf_file) || !file.size(sdf_file)) {
      return(invisible(NULL))
    }
    remaining <- setdiff(cids, got_cids)
    if (!length(remaining)) {
      return(invisible(NULL))
    }
    message("Ligand preparation source: ", label,
      " (remaining CID: ", length(remaining), ")")
    lst_sdf[[ label ]] <<- sdf_file
    run_meeko(label, sdf_file, remaining, suffix = "raw_meeko")
    remaining <- setdiff(cids, got_cids)
    run_rdkit_meeko(label, sdf_file, remaining)
    remaining <- setdiff(cids, got_cids)
    run_obabel(label, sdf_file, remaining)
    invisible(NULL)
  }

  if (!is.null(sdf_3d)) {
    run_source("user_3d_sdf", sdf_3d)
  }

  if (use_pubchem_3d && length(setdiff(cids, got_cids))) {
    q3d <- vinaFuns$query_pubchem_3d_sdfs(
      setdiff(cids, got_cids),
      dir_save = paste0(dir_save, "_PubChem3D"),
      filename = add_filename_suffix("pubchem_3d.sdf", basename(dir_save)),
      record_type = "3d",
      overwrite = overwrite_pubchem_3d
    )
    lst_status[[ "pubchem_3d_download" ]] <- q3d$status
    run_source("pubchem_3d", q3d$file)
  }

  if (use_obgen && length(setdiff(cids, got_cids))) {
    message("PubChem 3D source did not yield PDBQT for all ligands; start obgen fallback.")
    sdf_obgen <- cal_3d_sdf(sdf_2d, cl = cl, overwrite = overwrite_obgen)
    run_source("obgen_3d", sdf_obgen)
  }

  if (rdkit_fallback && length(setdiff(cids, got_cids))) {
    message("Use RDKit from 2D SDF fallback for remaining ligands.")
    run_source("rdkit_from_2d_sdf", sdf_2d)
  }

  status <- if (length(lst_status)) {
    lst_status <- lapply(lst_status, function(dat) {
      dat <- tibble::as_tibble(dat)
      if ("CID" %in% colnames(dat)) {
        dat$CID <- as.character(dat$CID)
      }
      dat
    })
    dplyr::bind_rows(lst_status)
  } else {
    tibble::tibble()
  }
  stat <- tibble::tibble(
    Strategy = names(lst_sdf),
    SDF = unname(unlist(lst_sdf)),
    Prepared_total = vapply(names(lst_sdf),
      function(nm) {
        as.integer(sum(got_cids %in% vinaFuns$sdf_ids(lst_sdf[[ nm ]])))
      }, integer(1))
  )
  list(
    pdbqt = got_files,
    pdbqt.cid = got_cids,
    status = status,
    stat = stat,
    sdf_files = lst_sdf
  )
}

select_files_by_grep <- function(files, pattern){
  vapply(files,
    function(file) {
      any(grepl(pattern, readLines(file), TRUE))
    }, logical(1))
}

filter_pdbs <- function(files, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS") {
  files[ select_files_by_grep(files, pattern) ]
}

prepare_receptor <- function(files, mkdir.pdbqt = "protein_pdbqt", 
  use = c("prepare_receptor", "mk_prepare_receptor.py"))
{
  dir.create(mkdir.pdbqt, FALSE)
  use <- match.arg(use)
  file <- lapply(files,
    function(file) {
      if (!is.null(file)) {
        newfile <- file.path(mkdir.pdbqt, paste0(basename(file), "qt"))
        if (use == "prepare_receptor") {
          cdRun(
            paste0(pg("prepare_receptor"), " -r ", file, " -o ", newfile, " -A checkhydrogens -e True -U nphs_lps_waters_nonstdres ")
          )
        } else if (use == "mk_prepare_receptor.py") {
          name <- tools::file_path_sans_ext(basename(file))
          extra <- " -v --box_size 20 20 20 --box_center 0 0 0 "
          cdRun(glue::glue("{pg('mk_prepare_receptor.py')} -i {file} -o {name} --delete_residues -p {extra}"))
        }
        return(newfile)
      }
    })
  file <- file[ !vapply(file, is.null, logical(1)) ]
  file[ vapply(file, file.exists, logical(1)) ]
}

setMethod("set_remote", signature = c(x = "job_vina"),
  function(x, wd = paste0("~/", x@sig, "_vina_space"))
  {
    if (!grpl(wd, "^[~/]")) {
      stop('!grpl(wd, "^[~/]").')
    }
    rem_dir.create(wd)
    x$wd <- wd
    return(x)
  })

cal_3d_sdf <- function(sdf, group = 10, cl = NULL, overwrite = FALSE) {
  output <- add_filename_suffix(sdf, "3D")
  if (!overwrite && file.exists(output) && file.size(output) > 0L) {
    message("Use existing obgen 3D SDF: ", output)
    return(output)
  }
  db <- sep_list(readLines(sdf), sep = "^\\$\\$\\$\\$")
  valid <- TRUE
  db <- lapply(db,
    function(lines) {
      if (lines[1] == "") {
        valid <<- FALSE
        for (i in seq_along(lines)) {
          if (lines[i] != "") {
            break
          }
        }
        lines <- lines[-seq_len(i - 1)]
      }
      lines
    })
  if (!valid) {
    writeLines(unlist(db), sdf)
  }
  if (length(db) <= group) {
    cdRun(pg("obgen"), " ", sdf, " -ff UFF > ", output)
  } else {
    groupSeqs <- grouping_vec2list(seq_along(db), group, TRUE)
    file_groups <- add_filename_suffix(
      sdf, seq_along(groupSeqs)
    )
    N <- 0L
    lapply(groupSeqs,
      function(ns) {
        N <<- N + 1L
        writeLines(unlist(db[ns]), file_groups[N])
      })
    output_groups <- add_filename_suffix(output, seq_len(N))
    pg <- pg('obgen')
    pbapply::pblapply(seq_len(N), cl = cl,
      function(n) {
        cdRun(glue::glue("{pg} {file_groups[n]} -ff UFF > {output_groups[n]}"))
      })
    lines <- unlist(lapply(output_groups, readLines))
    writeLines(lines, output)
  }
  output
}

.pymol_select_polymer.protein <- function(file, newfile) {
  cmd <- glue::glue(
    "select protein, polymer.protein; save {newfile}, protein"
  )
  cdRun(glue::glue("{pg('pymol')} -c -Q -d 'load {file}; {cmd}; quit'"))
}

.pymol_select_chains <- function(file, chains, newfiles = add_filename_suffix(file, chains))
{
  cmd <- glue::glue(
    "select chain{chains}, chain {chains}; save {newfiles}, chain{chains}"
  )
  cmd <- paste0(cmd, collapse = "; ")
  cdRun(glue::glue("{pg('pymol')} -c -Q -d 'load {file}; {cmd}; quit'"))
}

setMethod("res", signature = c(x = "job_vina"),
  function(x, meta,
    use = "Ingredient_name", target = "Ingredient.name", get = "Herb_pinyin_name")
  {
    data <- dplyr::select(x@tables$step5$res_dock, -dir, -file)
    data <- dplyr::relocate(data, hgnc_symbol, Ingredient_name, Affinity)
    if (!missing(meta)) {
      meta <- dplyr::filter(meta, !!rlang::sym(target) %in% !!data[[ use ]])
      meta <- dplyr::distinct(meta, !!rlang::sym(target), !!rlang::sym(get))
      meta <- dplyr::group_by(meta, !!rlang::sym(target))
      fun <- function(x) paste0(x, collapse = "; ")
      meta <- dplyr::reframe(meta, .get = fun(!!rlang::sym(get)))
      data <- map(data, use, meta, target, ".get", col = get)
    }
    data
  })

# deprecated:
# if (bdb) {
#   if (file.exists(bdb_file)) {
#     message("Database `BindingDB` used for pre-filtering of docking candidates.")
#     bdb <- ld_cutRead(bdb_file, c("PubChem CID", "PDB ID(s) of Target Chain"))
#     bdb <- dplyr::filter(bdb, `PubChem CID` %in% object(x)$cids)
#     colnames(bdb) <- c("PubChem_id", "pdb_id")
#     x@params$bdb_compounds_targets <- bdb
#     bdb <- nl(bdb$PubChem_id, strsplit(bdb$pdb_id, ","))
#     bdb <- lst_clear0(bdb)
#     bdb <- lapply(bdb, function(v) v[ v %in% x$targets_annotation$pdb ])
#     bdb <- lst_clear0(bdb)
#     if (!length(bdb)) {
#       stop("No candidates found in BindingDB.")
#     }
#     x$dock_layout <- bdb
#     .add_internal_job(.job(method = "Database `BindingDB` used for pre-filtering of docking candidates.",
#         cite = "[@BindingdbIn20Gilson2016]"))
#   }
# }
