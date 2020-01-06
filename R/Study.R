#' Study Factory
#'
#' Intanciates an object that extends the epimedtools internal class Study_abstract.
#'
#' This function intanciates an object that extends the epimedtools internal abstract class Study_abstract. You can learn more about Study_abstract class by typing ?epimedtools::Study_abstract
#' 
#' @param cache_filename A character string that describes the study cache file name.
#' @param ... Parameters that are given to the Study contructor
#' @examples
#' # create study and load data from GSE57831_series_matrix.txt.gz contained in the package
#' study = create_study()
#' study$series_matrix_filename = system.file(
#'    "extdata/GSE26471", 
#'    "GSE26471_series_matrix.txt.gz", 
#'    package = "epimedtools"
#'  )  
#' # str(study$get_gset(dest_dir="data"))
#' 
#' # load data from GEO
#' study = create_study()
#' study$gse = "GSE26471" # 1 sample, expression
#' # study$gse = "GSE42707" # 1 sample
#' # study$gse = "GSE56382" # 2 samples
#' # study$gse = "GSE49585" # 2 samples
#' # str(study$get_gset(dest_dir="data"))
#' 
#' # cache study
#' study$save("/tmp/tmp_cached_study.rds")
#' # load a cached study
#' study = create_study("/tmp/tmp_cached_study.rds")
#' @return It returns an object that extends the Study_abstract class.
#' @export
create_study = function(cache_filename, ...) {
  study = Study_abstract()
  if (!missing(cache_filename)) {
    s = readRDS(cache_filename)
    study$plaform = s$plaform
    study$data    = s$data   
    study$exp_grp = s$exp_grp
  }
  return(study)
}


#' A Reference Class to represent a multi-omic study.
#'
#' This class allow to manipulate useful concepts involved in a multi-omics
#' study. It an abstract class.
#' 
#' @field plaform The Table version of the GEO platform description.
#' @field platform_filenames A character string that describes the platform file name.
#' @field cache_filename A character string that describes the study cache file name.
#' @field gset The GEO set of data as it is return by getGEO.
#' @field gse A character string that describes the GEO accession number.
#' @field series_matrix_filename A character string that describes the GEO serie matrix family file name.
#' @field data A matris of experiemental data. Each column is a sample, named by it's sample name, and each row is a probes.
#' @field platform_name A character string that describes the platform used to perform.
#' @field exp_grp A data.frame that describe the experimental grouping. Each line is a sample.
#' @field cel_filedirs A character string that describes the directory taht contain associated cel files.
#' @export
#' @importFrom Biobase exprs
#' @importFrom affy justRMA
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @importFrom beanplot beanplot
#' @import methods
Study_abstract = setRefClass(
  "Study_abstract",
  fields = list(
    "dest_dir" = "character",
    cel_filedirs="character",
    cel_files="character",
    "data" = "ANY",
    "exp_grp" = "ANY",
    "platform" = "ANY",
    "platform_filename" = "character",
    "platform_name" = "character",
    "cache_filename" = "character",
    "gset" = "ANY",
    "gse" = "character",
    "series_matrix_filename" = "character", 
    "stuffs" = "list"
  ),
  methods = list(
    get_data = function(CACHE=FALSE, ...) {
      if (is.null(dim(.self$data))) {
        if (length(.self$gse) != 0) {
          .self$get_gset(...)
        } else {
          .self$set_data(...)          
        }
        if (CACHE) {
          .self$save(FORCE=TRUE)
        }
      }
      return(data)
    },
    get_platform_name = function(...) {
      if (length(.self$platform_name) == 0) {
        if (is.null(dim(.self$gset))) {
          if (length(.self$gse) != 0) {
            .self$get_gset(...)
          }
        }
      }
      return(platform_name)
    },
    get_exp_grp = function(...) {
      if (is.null(dim(.self$exp_grp))) {
        if (is.null(dim(.self$gset))) {
          if (length(.self$gse) != 0) {
            .self$get_gset(...)
          }
        }
      }
      return(exp_grp)
    },
    initialize = function(cfn, MEMOISATION=FALSE) {
      "Constructor."
      if (!missing(cfn)) {
        .self$cache_filename = cfn
        if (!file.exists(.self$cache_filename)) {
          cache_dir = paste(rev(rev(
            strsplit(.self$cache_filename, "/")[[1]]
          )[-1]), collapse = "/")
          dir.create(cache_dir, showWarnings = FALSE, recursive =
                       TRUE)
          .self$save(FORCE=TRUE)
        } else {
          if (MEMOISATION) {
            readRDS_funcname = "mreadRDS"            
          } else {
            readRDS_funcname = "readRDS"
          }
          print("reifying study...")
          s2 = get(readRDS_funcname)(.self$cache_filename)
          # print(.self)
          for (f in names(s2$getRefClass()$fields())) {
            .self[[f]] = s2[[f]]
          }
          .self$cache_filename = cfn
          print("done.")
        }
      }
    },
    save = function(cache_filename, FORCE=FALSE) {
      "Writes the study on the disk."
      if (!FORCE) {
        check_study(.self)        
      }
      if (!missing(cache_filename)) {
        .self$cache_filename = cache_filename
      }
      if (length(.self$cache_filename) == 1) {
        print("caching study...")
        tmp_list = list()        
        for (f in names(.self$getRefClass()$fields())) {
           tmp_list[[f]] = .self[[f]]
        }
        saveRDS(tmp_list, .self$cache_filename)
        print("done.")
      }
    },
    get_platform = function(CACHE = FALSE, MEMOISE = TRUE, dest_dir = "data/platforms") {
      "Computes if not and returns the platform field."
      if (is.null(dim(.self$platform))) {
        if (length(.self$platform_filename) == 1) {
          if (file.exists(.self$platform_filename)) {
            print(paste("Launching platform informations from ", .self$platform_filename, "...", sep = ""))
            if (MEMOISE) {
              .self$platform = mtgetGEO(filename = .self$platform_filename, getGPL = FALSE)
            } else {
              .self$platform = Table(getGEO(
                filename = .self$platform_filename, getGPL = FALSE
              ))
            }
            print("done.")
          } else {
            stop(paste("No", .self$platform_filename, "file."))
          }
        } else {
          print(paste(
            "Launching ", .self$get_platform_name(), " from GEO...", sep =
              ""
          ))
          dir.create(dest_dir, showWarnings = FALSE, recursive =TRUE)
          if (MEMOISE) {
            .self$platform = mtgetGEO(.self$get_platform_name(), getGPL = FALSE, destdir =
                                        dest_dir)
          } else {
            .self$platform = Table(getGEO(
              .self$get_platform_name(), getGPL = FALSE, destdir = dest_dir
            ))
            .self$platform_filename = paste(dest_dir, "/", .self$get_platform_name(), ".soft", sep="")
          }
        }
        .self$platform$ID = as.character(.self$platform$ID)
        rownames(.self$platform) = .self$platform$ID
        if (CACHE) {
          .self$save(FORCE=TRUE)
        }
      }
      return(.self$platform)
    },
    set_data = function(hook, method_to_read_cel_files, ...) {
            "Set the data field and update experiment grouping."
      if (!missing(hook)) {
        get(hook)(...)
      }
      if (missing(method_to_read_cel_files)) {
        method_to_read_cel_files = justRMA        
      }
      if (length(.self$cel_files) == 0) {
        tmp_self_cel_files = lapply(.self$cel_filedirs, get_cel_filenames)
        orig = unlist(sapply(1:length(.self$cel_filedirs), function(i) {
          split1 = unique(unlist(strsplit(.self$cel_filedirs[i], "/")))  
          split2 = unique(unlist(strsplit(.self$cel_filedirs[-i], "/")))  
          n = paste(split1[!(split1 %in% split2)], collapse="_")
          n = gsub("[.]", "_", n)
          n = gsub("__*", "_", n)
          n = gsub("^_", "", n)
          rep(n, length(tmp_self_cel_files[[i]]))
        }))
        .self$cel_files = unlist(tmp_self_cel_files)        
        cel_files_short = unlist(lapply(.self$cel_filedirs, get_cel_filenames, full.names=FALSE))
      } else {
        .self$cel_files = path.expand(.self$cel_files)
        orig = .self$cel_files
        cel_files_short = sapply(.self$cel_files, function(cel_file) {
          rev(strsplit(cel_file, "/")[[1]])[1]
        })
      }
      .self$cel_files = sapply(.self$cel_files, function(cel_file) {
        if (substr(cel_file, 1, 1) != "/") {
          return(paste(getwd(), "/", cel_file, sep=""))
        } else {
          return(cel_file)
        }
      })
      .self$cel_files = .self$cel_files[!duplicated(cel_files_short)]
      orig = orig[!duplicated(cel_files_short)]
      sample_names = cel_files_short[!duplicated(cel_files_short)]
      tmp_exp_grp = data.frame(orig=orig)
      rownames(tmp_exp_grp) = simplify_sample_names(sample_names) 
      if (!is.null(dim(.self$get_exp_grp()))) {
        tmp_exp_grp = fuse_exp_grp(.self$get_exp_grp(), tmp_exp_grp)
      }
      .self$exp_grp = tmp_exp_grp 
      .self$data = exprs(method_to_read_cel_files(filenames=.self$cel_files, celfile.path=""))
      colnames(.self$data) = simplify_sample_names(colnames(.self$data))
    },
    # Analysis method dealing with Study_abstract class intances...
    plot_qc = function(method = "boxplot", ...) {
      "Plot the quality control of the study."
      get(method)(t(na.omit(.self$get_data())) ~ colnames(.self$get_data()), las =
                    2, ...)
    },
    get_cel_files = function(dest_dir="data", ...) {
      "Retrieve .CEL.gz from NCBI GEO eweb site"
      gsms = as.character(.self$get_gset(dest_dir=dest_dir, ...)@phenoData@data$geo_accession)
      basedir = paste(dest_dir, "/", .self$gse, "/raw", sep="")
      return(get_gsm(gsms, basedir))
    },
    get_gset = function(CACHE=TRUE, MEMOISE=FALSE, dest_dir="data") {
      "Computes (if not yet done) and returns the gset field."
      if (is.null(dim(.self$gset))) {
        if (length(.self$series_matrix_filename) == 1) {
          if (file.exists(.self$series_matrix_filename)) {
            print(paste(
              "Launching data from ", .self$series_matrix_filename, "...", sep = ""
            ))
            if (MEMOISE) {
              .self$gset = mgetGEO(filename = .self$series_matrix_filename, getGPL = FALSE)
            } else {
              .self$gset = getGEO(filename = .self$series_matrix_filename, getGPL = FALSE)
            }
            print("done.")
          } else {
            stop(paste("No", .self$series_matrix_filename, "file."))
          }
        } else if (length(.self$gse) == 1) {
          print(paste("Launching ", .self$gse, " from GEO...", sep = ""))
          dest_dir_gse = paste(dest_dir, "/", .self$gse, sep="")
          l = list.files(dest_dir_gse)
          l_grep = grep("series_matrix.txt.gz", l)
          if (length(l_grep) > 0) {
            .self$series_matrix_filename = paste(dest_dir_gse, "/", l[l_grep[1]], sep="")
            return(.self$get_gset(CACHE=CACHE, MEMOISE=MEMOISE, dest_dir=dest_dir)) 
          }          
          dir.create(dest_dir_gse, showWarnings = FALSE, recursive = TRUE)
          if (MEMOISE) {
            tmp_gset = mgetGEO(.self$gse, getGPL = FALSE, destdir=dest_dir_gse)
          } else {
            tmp_gset = getGEO(.self$gse, getGPL = FALSE, destdir=dest_dir_gse)
          }
          if (length(tmp_gset) == 1) {
            .self$gset = tmp_gset[[1]]
            .self$series_matrix_filename = paste(dest_dir_gse, "/", names(tmp_gset)[1], sep = "")
          } else {
            if (length(.self$platform_name) != 0) {
              grep_pf = grep(.self$platform_name, names(tmp_gset))
              if (length(grep_pf) == 1) {
                warning(paste(.self$gse, " is a SuperSeries.", sep = ""))
                .self$gset = tmp_gset[[grep_pf]]
                tmp_series_matrix_filename = names(tmp_gset)[grep_pf]
                .self$series_matrix_filename = paste(dest_dir_gse, "/", tmp_series_matrix_filename, sep = "")                
              } else {
                stop(paste(.self$gse, " is a SuperSeries with a format that is not yet supported.", sep = ""))                
              }
            } else {
              stop(paste(.self$gse, " is a SuperSeries, platform_name needs to be define to avoid any ambiguity.", sep = ""))
            }
          }
          print(paste("done. File locally cached here: ", .self$series_matrix_filename, sep=""))
        } else {
          stop("You need to define gse (length 1) field to get it from GEO")
        }
        .self$platform_name = .self$gset@annotation
        .self$exp_grp = .self$gset@phenoData@data
        .self$data = exprs(.self$gset)
        if (CACHE) {
          .self$save(FORCE=TRUE)
        }
      }
      return(.self$gset)
    },
    get_ratio = function(ctrl_sample_names, method = "mean", ...) {
      "Normalize data accross probes using method *method* according to a reference poulation describe by  *ctrl_sample_names* its sample names vector."
      if (missing(ctrl_sample_names)) {
        d = .self$get_data(...)
      } else {
        d = .self$get_data(...)[, ctrl_sample_names]
      }
      ctrl_op = apply(d, 1, get(method))        
      ratio = .self$get_data(...) / ctrl_op
      return(ratio)
    },
    gs_to_probe = function(gene, ALL_PF_COL=FALSE, pf_col_name="Gene Symbol", col_delimiter="[ /]+", ...) {
      "Return probes name for a given gene."
      pf = .self$get_platform(...)
      if (ALL_PF_COL) {
        probes = unique(unlist(sapply(colnames(pf), function(cn) {
          rownames(pf[grep(gene, pf[[cn]]), ])
        })))
      } else {
        probes = rownames(pf[grep(gene, pf[[pf_col_name]]), ])    
        probes = sapply(probes, function(probe) {
          entry = as.character(pf[probe, pf_col_name])
          # print(entry)
          if (gene %in% unlist(strsplit(entry, col_delimiter))) {
            return(probe)
          } else {
            # print(NULL)
            return(NULL)
          }
        })
        probes = unlist(probes)
        # print(pf[probes, pf_col_name])
      }
      return(probes)
    },
    get_probe_gene_tab = function(genes, DEEP_SEARCH=FALSE, pf_col_name="Gene Symbol", col_delimiter="[ /]+", ...) {
      "Build a gene / probe table from a gene list."
      pf = .self$get_platform(...)
      probe_gene_tab = lapply(genes, function(gene) {
        probe = gs_to_probe(gene, pf_col_name=pf_col_name, col_delimiter=col_delimiter, ...)
        if (length(probe) == 0 & DEEP_SEARCH) {
            probe = gs_to_probe(gene, ALL_PF_COL=TRUE, col_delimiter=col_delimiter, ...)
        }
        if (length(probe) == 0) {
          return(NULL)
        }      
        df = data.frame(probe=probe)
        df$gene = gene
        return(df)
      })
      probe_gene_tab = do.call(rbind,probe_gene_tab)
      probe_gene_tab$probe = as.character(probe_gene_tab$probe)
      probe_gene_tab$gene = as.character(probe_gene_tab$gene)
      probe_gene_tab = data.frame(probe_gene_tab, stringsAsFactors = FALSE)
      return(probe_gene_tab)
    },
    # do_sw = function(sample_names, probe_names) {
    #   "Performs the Shapiro-Wilk test of normality over for each probe names.Perform shaanova test for a given `probe_name` and a given `exp_grp_key`."
    #   # data = .self$get_data()
    #   if (missing(probe_names)) {
    #     probe_names = rownames(data)
    #   }
    #   bar = data[probe_names, sample_names]
    #   apply()
    #
    #   foo = msapply(probe_names, function(probe_name) {
    #     s = shapiro.test(data[probe_name, sample_names])
    #     s_pval = -log10(s$p.value)
    #     return(s_pval)
    #   })
    #   idx_sample = rownames(.self$get_exp_grp())
    #   # Dealing with ratio data, reduce data to interesting probes and samples
    #   filtred_bp_data = .self$get_data()[probe_name, idx_sample]
    #   m = aov(filtred_bp_data~.self$get_exp_grp()[,exp_grp_key])
    #   # tests
    #   s = shapiro.test(residuals(m))
    #   shap = -log10(s$p.value)
    #   f_kc = m$coefficients[2]
    #   pval = -log10(summary(m)[[1]][["Pr(>F)"]][1])
    #   ret = list(shap=shap, pval=pval)
    #   return(ret)
    # },
    plot_pca_eig = function(pca_res, xlim, ...) {
      if (missing(xlim)) {
        xlim=c(0, min(length(pca_res$var), 50))
      }
      barplot(pca_res$pvar, xlab="PC", ylab="% of var", xlim=xlim, ...)
      barplot(pca_res$bs[,2], add=TRUE, col=adjustcolor(2, alpha.f=0.1))
      legend("topright", col=c(1,2), legend=c("eigenvalues", "broken-stick"), pch=15)
      # abline(h=pca_res$kaiser_crit)
    },
    plot_pca_pc = function(pca_res, pc1=1, pc2=2, exp_grp_key, col, LEGEND=TRUE, ...) {
      if (missing(exp_grp_key)) {
        col = 1
      }
      if (missing(col)) {
        col = as.factor(.self$get_exp_grp()[rownames(pca_res$x),][[exp_grp_key]])        
      }
      plot(pca_res$x[,c(pc1, pc2)], col=col, pch=16, 
        xlab = paste("PC", pc1, " (", round(pca_res$pvar[pc1]*100,2),"%)", sep=""), 
        ylab = paste("PC", pc2, " (", round(pca_res$pvar[pc2]*100,2),"%)", sep=""), 
        ...
      )
      if (LEGEND) {
        legend("topright", legend=sort(unique(col)), col=sort(unique(col)), pch=16)        
      }
    },
    do_pca = function(probe_names, sample_names, ...) {
      if (missing(probe_names)) {
        probe_names = rownames(.self$get_data())
      }
      if (missing(sample_names)) {
        sample_names = colnames(.self$get_data())
      }
      pca_data = .self$get_data()[probe_names, sample_names]
      # pca_res = FactoMineR::PCA(t(pca_data), graph=FALSE, ncp=8)
      pca_data[apply(t(pca_data), 2, var, na.rm=TRUE) != 0, ]
      pca_res = prcomp(t(pca_data[apply(pca_data, 1, var, na.rm=TRUE) != 0, ]), scale. = TRUE, ...)
      pca_res$var = pca_res$sdev * pca_res$sdev
      pca_res$pvar = pca_res$var/ sum(pca_res$var)
      pca_res$bs = broken_stick(length(pca_res$pvar))
      pca_res$bs_crit = min(which(pca_res$pvar < pca_res$bs[,2])) - 1
      pca_res$kaiser_crit = 1/length(pca_res$var)      
      return(pca_res)
    },
    pretreat_before_a_test = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, two_grp_test_func, 
        for_each_line_process_groups_func="for_each_line_process_groups_two_by_two_func", ...
      ) {
      " Perform a pretreatment on data, according to `probe_names` and `exp_grp` keys and factors, before processing the dataset." 
      # Check exp_grp keys
      if (missing(ctrl_key)) {
        stop("You need to specify a control column key.")
      }
      if (missing(case_key)) {
        case_key = ctrl_key
      }
      if (!(ctrl_key %in% colnames(.self$get_exp_grp()))) {
        stop(paste(ctrl_key, " is not a column name of exp_grp.", sep=""))
      }
      if (!(case_key %in% colnames(.self$get_exp_grp()))) {
        stop(paste(case_key, " is not a column name of exp_grp.", sep=""))
      }
      # Check exp_grp factors
      if (missing(ctrl_fctr)) {
        ctrl_fctr = na.omit(unique(.self$get_exp_grp()[[ctrl_key]]))
      }
      if (missing(case_fctr)) {
        case_fctr = na.omit(unique(.self$get_exp_grp()[[case_key]]))
      }
      if (case_key == ctrl_key) {
        case_fctr = case_fctr[which(!case_fctr %in% ctrl_fctr)]
      }
      # preprocess groups
      ctrl_list = lapply(ctrl_fctr, function(ctrl_f) {
        if (!(ctrl_f %in% unique(.self$get_exp_grp()[[ctrl_key]]))) {
          stop(paste(ctrl_f, " is not a factor of ", ctrl_key, " exp_grp column.", sep=""))
        }
        rownames(.self$get_exp_grp())[which(.self$get_exp_grp()[[ctrl_key]] == ctrl_f)]
      })
      case_list = lapply(case_fctr, function(case_f) {
        if (!(case_f %in% unique(.self$get_exp_grp()[[case_key]]))) {
          stop(paste(case_f, " is not a factor of ", case_key, " exp_grp column.", sep=""))
        }
        rownames(.self$get_exp_grp())[which(.self$get_exp_grp()[[case_key]] == case_f)]      
      })
      if (ctrl_key != case_key) {
        names(ctrl_list) = paste(ctrl_key, ctrl_fctr, sep="_")
        names(case_list) = paste(case_key, case_fctr, sep="_")        
      } else {
        names(ctrl_list) = paste(ctrl_fctr, sep="_")
        names(case_list) = paste(case_fctr, sep="_")        
      }
      # names(ctrl_list) = paste(ctrl_key, ctrl_fctr, sep="_")
      # names(case_list) = paste(case_key, case_fctr, sep="_")
      # if (length(ctrl_fctr) > 1) {
      #   names(ctrl_list) = simplify_factor_names(names(ctrl_list))
      # }
      # if (length(case_fctr) > 1) {
      #   names(case_list) = simplify_factor_names(names(case_list))
      # }
      # print(ctrl_fctr)
      # print(ctrl_list)
      # print(case_fctr)
      # print(case_list)
      # Filter data
      tmp_data = .self$get_data()
      if (missing(probe_names)) {
        probe_names = rownames(tmp_data)
      }
      tmp_data = tmp_data[probe_names, unique(unlist(c(ctrl_list, case_list)))]
      if (length(probe_names) == 1) {
        tmp_data = t(as.matrix(tmp_data))
      }
      # print(tmp_data)
      # Go!
      for_each_line_process_groups_two_by_two_func = function(line, ctrl_list, case_list, ...) {
        ctrl_ret = lapply(ctrl_list, function(ctrl_samples) {
          case_ret = lapply(case_list, function(case_samples, ctrl_samples, two_grp_test_func, ...) {
            ctrl = line[ctrl_samples]
            case = line[case_samples]
            # print("____________________________")
            # print(ctrl)
            # print("            ____            ")
            # print(case)
            # print("____________________________")
            return(two_grp_test_func(ctrl, case, ...))
          }, ctrl_samples, two_grp_test_func=two_grp_test_func, ...)
          case_ret = unlist(case_ret, recursive=FALSE)
          
          # print("______________________ ")
          # print(length(case_list))
          # print(rep(paste(names(case_list)), length(case_list))), "vs",
          # print(case_ret)
          return(case_ret)
        })
        ctrl_ret = unlist(ctrl_ret, recursive=FALSE)
        # print("______________________ ")
        # print(nb_fields)
        # nb_fields = length(ctrl_ret) / (length(case_list) * length(ctrl_list))
        # colnames = sapply(names(ctrl_list), function(ctrl_name) {
        #   sapply(names(case_list), function(case_name) {
        #     paste(rep(case_name, nb_fields), "vs", ctrl_name, sep="_")
        #   })
        # })
        # print(colnames)
        # print(ctrl_ret)
        return(ctrl_ret)
      }
      for_each_line_process_all_groups_func = function(line, ctrl_list, case_list, ...) {
        unified_list = c(ctrl_list, case_list)
        # print(unified_list)
        samples = unlist(unified_list)
        filtred_bp_data = line[samples]
        filtred_bp_factors = as.factor(unlist(lapply(names(unified_list),function(n){rep(n, length(unified_list[[n]]))})))
        # print("-------------------------------------------------")
        # print(names(unified_list))
        # print(length(samples))
        # print(length(filtred_bp_data))
        # print(filtred_bp_data)
        # print(length(filtred_bp_factors))
        # print(filtred_bp_factors)
        two_grp_test_func(filtred_bp_data, filtred_bp_factors, "prb", "gn", ...)
        # filtred_bp_factors, probe_name, gene_name
        return(NULL)
      }
      test_res = monitored_apply(tmp_data, 1, function(line, ctrl_list, case_list, ...) {
        get(for_each_line_process_groups_func)(line=line, ctrl_list=ctrl_list, case_list=case_list, ...)
      }, ctrl_list=ctrl_list, case_list=case_list, ...)
      # print(test_res)
      if (!is.null(test_res)) {
        test_res = do.call(rbind, test_res)
        test_res = data.frame(lapply(data.frame(test_res, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
      }
      return(test_res)
    },
    do_mw_test = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, ...) {
      " Perform a Mann-Whitney test to detect shifted (`two.sided`, `less` or `greater`) expression groups for a given `probe_name` and a given `exp_grp_key`." 
      mw_func = function(ctrl, case, alternative=NULL, PLOT=FALSE)  {
        ctrl_median = median(ctrl)
        ctrl_mean = mean(ctrl)
        case_median = median(case)
        case_mean = mean(case)
        s = sign(case_mean - ctrl_mean)
        # print(s)
        lr = s * abs(case_mean - ctrl_mean)
        fc = logratio2foldchange(lr)
        # fc = s * max(case_mean/ctrl_mean, ctrl_mean/case_mean)
        # d_med = case_median - ctrl_median
        if (is.null(alternative)) {
          if (s >= 0) {            
            alternative="less"
          } else {
            alternative="greater"            
          }
        }
        mw = suppressWarnings(wilcox.test(ctrl, case, alternative=alternative))
        if (PLOT) {
          beanplot(ctrl, case, col=(mw$p.value < 0.05) + 1, main = mw$p.value, log="")
        }
        return(list(logratio=lr, foldchange=fc, mw_pval = mw$p.value))#, d_med = d_med))
      }
      mw_res = .self$pretreat_before_a_test(probe_names=probe_names, ctrl_key=ctrl_key, case_key=case_key, ctrl_fctr=ctrl_fctr, case_fctr=case_fctr, two_grp_test_func=mw_func, ...)
      # colnames(mw_res) = simplify_column_names(colnames(mw_res))
      # colnames(mw_res) = simplify_factor_names(colnames(mw_res), "_")
      mw_pval_colnames = colnames(mw_res)[grep("mw_pval", colnames(mw_res))]
      adj_pvals = apply(t(mw_res[,mw_pval_colnames]), 1, p.adjust, method="BH")
      colnames(adj_pvals) = paste(mw_pval_colnames, "_adj", sep="")
      mw_res = cbind(mw_res, adj_pvals)
      return(mw_res)
    },
    do_fast_fc = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr) {
      " Perform a Mann-Whitney test to detect shifted (`two.sided` `less` or `greater`) expression groups for a given `probe_name` and a given `exp_grp_key`." 
      fast_fc_func = function(ctrl, case)  {
        ctrl_median = median(ctrl)
        ctrl_mean = mean(ctrl)
        case_median = median(case)
        case_mean = mean(case)
        s = sign(case_mean - ctrl_mean)
        # print(s)
        fc = s * max(case_mean/ctrl_mean, ctrl_mean/case_mean)
        return(list(fc=fc))
      }
      mw_res = .self$pretreat_before_a_test(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, two_grp_test_func=fast_fc_func)
      colnames(mw_res) = simplify_column_names(colnames(mw_res))
      colnames(mw_res) = simplify_factor_names(colnames(mw_res), "_")
      return(mw_res)
    },
    do_gm2sd_analysis = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, ctrl_thres_func=m2sd, case_value_func=mean, comp_func=get("<"), nb_perm=100, MONITORED=FALSE) {
      "Performs permutation test to detect right shifted expression groups for two given groups."
      gm2sd_func = function(ctrl, case, ctrl_thres_func, case_value_func, comp_func, nb_perm, MONITORED)  {
        ctrl_thres = ctrl_thres_func(ctrl)
        freq = sum(comp_func(ctrl_thres, case)) / length(case)
        if (MONITORED) {
          apply_function_name = "monitored_apply"
        } else {
          apply_function_name = "apply"
        }
        perm_data = get(apply_function_name)(t(0:nb_perm), 2, function(perm){
          if (perm > 0) {
            set.seed(perm)
            pooled_values = sample(c(ctrl, case))
            ctrl = pooled_values[1:length(ctrl)]
            case = pooled_values[(length(ctrl)+1):length(pooled_values)]
          }
          # ctrl
          ctrl_thres = ctrl_thres_func(ctrl)
          case_value = case_value_func(case)
          comp = comp_func(ctrl_thres, case_value)
          return(comp)
        })
        # print(perm_data)
        idx = perm_data[1]
        if (nb_perm > 0) {
          pval = sum(perm_data[2:(nb_perm+1)])/nb_perm
          pval[pval == 0] = 1 / (10*nb_perm)
        } else {
          pval = rep(1, length(idx))
        }
        return(list(idx=idx, pval=pval, freq=freq))
        # return(list(freq=freq))
      }
      # process_group
      ret = .self$pretreat_before_a_test(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, two_grp_test_func=gm2sd_func, ctrl_thres_func=ctrl_thres_func, case_value_func=case_value_func, comp_func=comp_func, nb_perm=nb_perm, MONITORED=MONITORED)
      return(ret)
    },
    # do_mw_test = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, alternative, PLOT=FALSE) {
    plot_boxplot2 = function(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, ylim, las=2, col="grey", border="black", bp_function_name="boxplot") {
      "Draw the box plot for a given `probe_name` and a given `exp_grp_key`."
      boxplot_func = function(filtred_bp_data, filtred_bp_factors, probe_name, gene_name, ylim=ylim, las=2, col="grey", border="black", bp_function_name="boxplot", ...)  {
        # Box plots
        if (missing(ylim)) {
          ylim = c(min(filtred_bp_data), max(filtred_bp_data))
        }
        if (missing(gene_name)) {
          main = probe_name
        } else {
          main = paste(gene_name, probe_name, sep="@")
        }
        # boxplot_filename <- paste(study_dirname, "/", gene, "_", probe_name, "_", exp_grp_key, ".pdf", sep="")
        # pdf(file=boxplot_filename, height=10, width=10)# open jpeg device with specified dimensions
        bp_function = get(bp_function_name)
        bp_function(filtred_bp_data~filtred_bp_factors, # open boxplot
        # what=c(1,1,1,0),
          ylim = ylim,                     # fixed vertical scale
          col = col, border = border,  # colors of inside and border of box
          las = las,                         # written vertically
          # xlab = exp_grp_key,
          main = main, # title
          ...
        )
        # dev.off()
      }
      ret = .self$pretreat_before_a_test(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, two_grp_test_func=boxplot_func, for_each_line_process_groups_func="for_each_line_process_all_groups_func")
      return(ret)
    },
    do_anova = function(probe_names, samples_names, exp_grp_key) {
      "Performs anova test for a given `probe_name` and a given `exp_grp_key`."
      # Check data
      tmp_data = .self$get_data()
      if (missing(probe_names)) {
        probe_names = rownames(tmp_data)
      }
      if (missing(samples_names)) {
        samples_names = colnames(tmp_data)
      }
      tmp_data = tmp_data[probe_names, samples_names]
      histo = .self$get_exp_grp()[samples_names,exp_grp_key]
      anova_res = apply(tmp_data, 1, function(line){
        m = aov(line~histo)
        # tests
        s = shapiro.test(residuals(m))
        shap_pval = s$p.value
        f_kc = m$coefficients[2]
        pval = -log10(summary(m)[[1]][["Pr(>F)"]][1])
        ret = list(shap=shap, pval=pval)
        return(ret)        
      })
    },
    
    # plot_m2s_analysis = function(m2s, histo, label_col_names="gene", xlim, ...) {
    #   h = histo
    #   nsd = m2s[[paste(h, "nsd", sep="_")]]
    #   pval = m2s[[paste(h, "p_nsd", sep="_")]]
    #   pval[pval==0] = 1/(1/min(pval[pval!=0]) + 1)
    #   # idx = nsd > 0
    #   # nsd = nsd[idx]
    #   # pval = pval[idx]
    #   # if (label_col_names %in% colnames(m2s)) {
    #   #   col = as.factor(m2s[[label_col_names]])
    #   #   col = col[idx]
    #   # } else {
    #   #   col=1
    #   # }
    #   if (missing(xlim)) {
    #     xlim = c(-4, log2(max(nsd)))
    #   }
    #   plot(log2(nsd), -log10(pval), main=h, pch=16, xlim=xlim, ...)
    #   # if (label_col_names %in% colnames(m2s)) {
    #   #   legend("bottomright", col=unique(col), legend=unique(col), pch=16)
    #   # }
    #   abline(v=log2(c(2,3)))
    #   abline(h=-log10(0.05))
    # },
    # do_m2s_analysis = function(probe_names, exp_grp_key, ctrl_name, nb_perm=100, MONITORED=FALSE, ...) {
    #   "Performs permutation test to detect right shifted expression groups for a given `probe_name` and a given `exp_grp_key`."
    #   # Before starting...
    #   # exp_grp = .self$get_exp_grp()
    #   tmp_data = .self$get_data()
    #   tmp_data = tmp_data[probe_names,rownames(exp_grp)[!is.na(exp_grp[[exp_grp_key]])]]
    #   # histos
    #   histo_names = na.omit(unique(exp_grp[[exp_grp_key]]))
    #   histo_names = histo_names[-which(histo_names == ctrl_name)]
    #   # perm_sample_names
    #   perm_sample_names = sapply(1:nb_perm, function(i) {
    #     set.seed(i)
    #     sample(colnames(tmp_data))
    #   })
    #   # do permutations
    #   if (MONITORED) {
    #     apply_function_name = "monitored_apply"
    #   } else {
    #     apply_function_name = "apply"
    #   }
    #   m2s_perm_data = get(apply_function_name)(t(0:nb_perm), 2, function(perm){
    #     cur_data = tmp_data
    #     if (perm > 0) {
    #       colnames(cur_data) = perm_sample_names[,perm]
    #     }
    #     # ctrl
    #     ctrl = t(apply(cur_data[probe_names, na.omit(rownames(exp_grp)[exp_grp[[exp_grp_key]] == ctrl_name])], 1, function(line) {
    #       mean = mean(line)
    #       sd = sd(line)
    #       return(c(mean, sd))
    #     }))
    #     colnames(ctrl) = c("ctrl_m", "ctrl_sd")
    #     # mean of histos
    #     means = sapply(histo_names, function(h) {
    #       sample_names = na.omit(rownames(exp_grp)[exp_grp[[exp_grp_key]] == h])
    #       ret = apply(cur_data[probe_names, sample_names], 1, mean)
    #       return(ret)
    #     })
    #     colnames(means) = paste(histo_names, "m", sep="_")
    #     m2s = cbind(ctrl, means)
    #     # nb of sd
    #     nsd = sapply(histo_names, function(h) {
    #       (m2s[,paste(h, "m", sep="_")] - m2s[,"ctrl_m"]) / m2s[,"ctrl_sd"]
    #     })
    #     colnames(nsd) = paste(histo_names, "nsd", sep="_")
    #     m2s = cbind(m2s, nsd)
    #     # fold change
    #     if (perm < 0) {
    #       fc = sapply(histo_names, function(h) {
    #         s = sign(m2s[,paste(h, "m", sep="_")] - m2s[,"ctrl_m"])
    #         ret = s * apply(cbind(m2s[,paste(h, "m", sep="_")] / m2s[,"ctrl_m"], m2s[,"ctrl_m"] / m2s[,paste(h, "m", sep="_")]), 1, max)
    #       })
    #       colnames(fc) = paste(histo_names, "fc", sep="_")
    #       m2s = cbind(m2s, fc)
    #     }
    #     # m > m+2sd ?
    #     bool = sapply(histo_names, function(h) {
    #       m2s[,paste(h, "nsd", sep="_")] > 2
    #     })
    #     colnames(bool) = paste(histo_names, "m2s", sep="_")
    #     m2s = cbind(m2s, bool)
    #     # return
    #     m2s = m2s[,sort(colnames(m2s))]
    #     return (data.frame(m2s))
    #   }, ...)
    #   # extract m2s data
    #   m2s = m2s_perm_data[[1]]
    #   # p_nsd
    #   perm_nsd = sapply(histo_names, function(h) {
    #     nsd_h = m2s_perm_data[[1]][, paste(h, "nsd", sep="_")]
    #     nsd_h_perm = lapply(1:nb_perm, function(i) {
    #       m2s_perm_data[[i+1]][, paste(h, "nsd", sep="_")]
    #     })
    #     nsd_h_perm = do.call(cbind, nsd_h_perm)
    #     bool_nsd_h_perm = apply(nsd_h_perm, 2, function(col) {
    #       col >= nsd_h
    #     })
    #     sum_bool_nsd_h_perm = apply(bool_nsd_h_perm, 1, sum) / nb_perm
    #     return(sum_bool_nsd_h_perm)
    #   })
    #   colnames(perm_nsd) = paste(histo_names, "p_nsd", sep="_")
    #   m2s = cbind(m2s, perm_nsd)
    #   # max h ?
    #   max_m = apply(t(m2s[, paste(histo_names, "m", sep="_")]), 2, function(l) {
    #     max_l = max(l)
    #     ret = histo_names[which(l == max_l)[1]]
    #     return(ret)
    #   })
    #   m2s = cbind(m2s, max_m)
    #   # sort
    #   m2s = m2s[,sort(colnames(m2s))]
    #   return(m2s)
    # },
    plot_boxplot = function(probe_name, exp_grp_key,  gene_name, ylim, las=2, col="grey", border="black", bp_function_name="boxplot", ...) {
      "Draw the box plot for a given `probe_name` and a given `exp_grp_key`."
    	idx_sample = rownames(.self$get_exp_grp())
    	# Dealing with ratio data, reduce data to interesting probes and samples
    	filtred_bp_data = .self$get_data()[probe_name, idx_sample]
    	# Box plots
      if (missing(ylim)) {
    	  ylim = c(min(filtred_bp_data), max(filtred_bp_data))
      }
      if (missing(gene_name)) {
        main = probe_name
      } else {
        main = paste(gene_name, probe_name, sep="@")
      }
	    # boxplot_filename <- paste(study_dirname, "/", gene, "_", probe_name, "_", exp_grp_key, ".pdf", sep="")
	    # pdf(file=boxplot_filename, height=10, width=10)# open jpeg device with specified dimensions
      bp_function = get(bp_function_name)
	    bp_function(filtred_bp_data~.self$get_exp_grp()[,exp_grp_key], # open boxplot
          # what=c(1,1,1,0),
	      ylim = ylim,                     # fixed vertical scale
	      col = col, border = border,  # colors of inside and border of box
	      las = las,                         # written vertically
	      xlab = exp_grp_key,
	      main = main, # title
        ...
	    )
	    # dev.off()
    }  
  )
)




