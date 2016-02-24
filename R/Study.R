#' Study Factory
#'
#' Intanciates an object that extends the epimedtools internal class Study_abstract.
#'
#' This function intanciates an object that extends the epimedtools internal abstract class Study_abstract. You can learn more about Study_abstract class by typing ?epimedtools::Study_abstract
#' 
#' @param cache_filename A character string that describes the study cache file name.
#' @param Study_RC_name A character string that specify the RC to use to instanciate the object.
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
create_study = function(cache_filename, Study_RC_name="Study_geo") {
  if (missing(cache_filename)) {
    study = get(Study_RC_name)()
  } else {
    study = get(Study_RC_name)(cache_filename)
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
#' @export
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @import methods
Study_abstract = setRefClass(
  "Study_abstract",
  fields = list(
    "platform" = "ANY",
    "platform_filename" = "character",
    "cache_filename" = "character"
  ),
  methods = list(
    initialize = function(cfn) {
      "Constructor."
      if (!missing(cfn)) {
        .self$cache_filename = cfn
        if (!file.exists(.self$cache_filename)) {
          cache_dir = paste(rev(rev(
            strsplit(.self$cache_filename, "/")[[1]]
          )[-1]), collapse = "/")
          dir.create(cache_dir, showWarnings = FALSE, recursive =
                       TRUE)
          .self$save()
        } else {
          print("reifiyng study...")
          s2 = readRDS(.self$cache_filename)
          # print(.self)
          for (f in names(s2$getRefClass()$fields())) {
            .self[[f]] = s2[[f]]
          }
          .self$cache_filename = cfn
          print("done.")
        }
      }
    },
    save = function(cache_filename) {
      "Writes the study on the disk."
      if (!missing(cache_filename)) {
        .self$cache_filename = cache_filename
      }
      if (length(.self$cache_filename) == 1) {
        print("caching study...")
        saveRDS(.self, .self$cache_filename)
        print("done.")
      }
    },
    get_platform = function(CACHE = FALSE, MEMOISE = TRUE, dest_dir = "data/platforms") {
      "Computes if not and returns the platform field."
      if (is.null(dim(.self$platform))) {
        if (length(.self$platform_filename) == 1) {
          if (file.exists(.self$platform_filename)) {
            print(
              paste(
                "Launching platform informations from ", .self$platform_filename, "...", sep =
                  ""
              )
            )
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
          dir.create(dest_dir, showWarnings = FALSE, recursive =
                       FALSE)
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
        if (CACHE) {
          .self$save()
        }
      }
      return(.self$platform)
    },
    get_platform_name = function(...) {
      "Abstract method that gives access to the platform name."
      stop("Call to the abstract method get_platform_name of the class Study_abstract.")
    },
    get_exp_grp = function(...) {
      "Abstract method that gives access to the experimental grouping."
      stop("Call to the abstract method get_exp_grp of the class Study_abstract.")
    },
    get_data = function(...) {
      "Abstract method that gives access to the data."
      stop("Call to the abstract method get_data of the class Study_abstract.")
    },
    # Analysis method dealing with Study_abstract class intances...
    plot_qc = function(method = "boxplot", ...) {
      "Plot the quality control of the study."
      get(method)(t(na.omit(.self$get_data())) ~ colnames(.self$get_data()), las =
                    2, ...)
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
    gs_to_probe = function(gene, ALL_PF_COL=FALSE, pf_col_name="Gene Symbol", ...) {
      "Return probes name for a given gene."
      pf = .self$get_platform(...)
      if (ALL_PF_COL) {
        probes = unique(unlist(sapply(colnames(pf), function(cn) {
          as.character(pf[grep(gene, pf[[cn]]), ]$ID)
        })))
      } else {
        probes = as.character(pf[grep(gene, pf[[pf_col_name]]), ]$ID)    
      }
    },
    get_probe_gene_tab = function(gene, DEEP_SEARCH=FALSE, pf_col_name="Gene Symbol", ...) {
      "Build a gene / probe table from a gene list."
      pf = .self$get_platform(...)
      probe_gene_tab = lapply(genes, function(gene) {
        probe = gs_to_probe(gene, pf_col_name=pf_col_name, ...)
        if (length(probe) == 0 & DEEP_SEARCH) {
            probe = gs_to_probe(gene, ALL_PF_COL=TRUE, ...)
        }
        if (length(probe) == 0) {
          return(NULL)
        }      
        df = data.frame(probe=probe)
        df$gene = gene
        return(df)
      })
      probe_gene_tab = do.call(rbind,probe_gene_tab)
    },
    do_sw = function(sample_names, probe_names) {
      "Performs the Shapiro-Wilk test of normality over for each probe names.Perform shaanova test for a given `probe_name` and a given `factor_name`."
      data = .self$get_data()
      if (missing(probe_names)) {
        probe_names = rownames(data)
      }
      bar = data[probe_names, sample_names]
      apply()

      foo = msapply(probe_names, function(probe_name) {
        s = shapiro.test(data[probe_name, sample_names])
        s_pval = -log10(s$p.value)
        return(s_pval)
      })
    	idx_sample = rownames(.self$get_exp_grp())
    	# Dealing with ratio data, reduce data to interesting probes and samples
    	filtred_bp_data = .self$get_data()[probe_name, idx_sample]
      m = aov(filtred_bp_data~.self$get_exp_grp()[,factor_name])
      # tests
      s = shapiro.test(residuals(m))
      shap = -log10(s$p.value)
      f_kc = m$coefficients[2]
      pval = -log10(summary(m)[[1]][["Pr(>F)"]][1])
      ret = list(shap=shap, pval=pval) 
      return(ret)
    },
    plot_m2s_analysis = function(m2s, histo) {
      h = histo
      nsd = m2s[[paste(h, "nsd", sep="_")]]
      idx = nsd > 0
      nsd = nsd[idx]
      pval = m2s[[paste(h, "p_nsd", sep="_")]]
      pval[pval==0] = 1/(1/min(pval[pval!=0]) + 1)
      pval = pval[idx]
      plot(log2(nsd), -log10(pval), main=h, pch=16, col=adjustcolor(1, alpha.f=0.4))
      abline(v=log2(c(2,3)))
      abline(h=-log10(0.05))      
    }, 
    do_m2s_analysis = function(probe_names, exp_grp_key, ctrl_name, nb_perm=1000, MONITORED=FALSE, ...) {
      "Performs permutation test to detect right shifted exprtession groups for a given `probe_name` and a given `factor_name`."
      # Before starting...
      exp_grp = .self$get_exp_grp()
      data = .self$get_data()
      data = data[probe_names,rownames(exp_grp)[!is.na(exp_grp[[exp_grp_key]])]]
      # histos
      histo_names = na.omit(unique(exp_grp[[exp_grp_key]]))
      histo_names = histo_names[-which(histo_names == ctrl_name)]
      # perm_sample_names
      perm_sample_names = sapply(1:nb_perm, function(i) {
        set.seed(i)
        sample(colnames(data))
      })
      # do permutations
      if (MONITORED) {
        apply_function_name = "monitored_apply"
      } else {
        apply_function_name = "apply"        
      }
      m2s_perm_data = get(apply_function_name)(t(0:nb_perm), 2, function(perm){
        cur_data = data
        if (perm > 0) {
          colnames(cur_data) = perm_sample_names[,perm] 
        }
        # ctrl
        ctrl = t(apply(cur_data[probe_names, na.omit(rownames(exp_grp)[exp_grp[[exp_grp_key]] == ctrl_name])], 1, function(line) {
          mean = mean(line)
          sd = sd(line)
          return(c(mean, sd))
        }))
        colnames(ctrl) = c("ctrl_m", "ctrl_sd")
        # mean of histos
        means = sapply(histo_names, function(h) {
          sample_names = na.omit(rownames(exp_grp)[exp_grp[[exp_grp_key]] == h])
          ret = apply(cur_data[probe_names, sample_names], 1, mean)
          return(ret)
        })
        colnames(means) = paste(histo_names, "m", sep="_")
        m2s = cbind(ctrl, means)
        # nb of sd
        nsd = sapply(histo_names, function(h) {
          (m2s[,paste(h, "m", sep="_")] - m2s[,"ctrl_m"]) / m2s[,"ctrl_sd"]
        })
        colnames(nsd) = paste(histo_names, "nsd", sep="_")
        m2s = cbind(m2s, nsd)
        # fold change
        if (perm < 0) {
          fc = sapply(histo_names, function(h) {
            s = sign(m2s[,paste(h, "m", sep="_")] - m2s[,"ctrl_m"])
            ret = s * apply(cbind(m2s[,paste(h, "m", sep="_")] / m2s[,"ctrl_m"], m2s[,"ctrl_m"] / m2s[,paste(h, "m", sep="_")]), 1, max)
          })
          colnames(fc) = paste(histo_names, "fc", sep="_")
          m2s = cbind(m2s, fc)    
        }
        # m > m+2sd ?
        bool = sapply(histo_names, function(h) {
          m2s[,paste(h, "nsd", sep="_")] > 2
        })
        colnames(bool) = paste(histo_names, "m2s", sep="_")
        m2s = cbind(m2s, bool)    
        # return
        m2s = m2s[,sort(colnames(m2s))]
        return (data.frame(m2s))
      }, ...)
      # extract m2s data
      m2s = m2s_perm_data[[1]]
      # p_nsd
      perm_nsd = sapply(histo_names, function(h) {
        nsd_h = m2s_perm_data[[1]][, paste(h, "nsd", sep="_")]
        nsd_h_perm = lapply(1:nb_perm, function(i) {
          m2s_perm_data[[i+1]][, paste(h, "nsd", sep="_")]
        })
        nsd_h_perm = do.call(cbind, nsd_h_perm)
        bool_nsd_h_perm = apply(nsd_h_perm, 2, function(col) {
          col >= nsd_h
        })
        sum_bool_nsd_h_perm = apply(bool_nsd_h_perm, 1, sum) / nb_perm
        return(sum_bool_nsd_h_perm)  
      })
      colnames(perm_nsd) = paste(histo_names, "p_nsd", sep="_")
      m2s = cbind(m2s, perm_nsd)
      # max h ?
      max_m = apply(t(m2s[, paste(histo_names, "m", sep="_")]), 2, function(l) {
        max_l = max(l)
        ret = histo_names[which(l == max_l)[1]]
        return(ret)
      })
      m2s = cbind(m2s, max_m)    
      # sort
      m2s = m2s[,sort(colnames(m2s))]
      return(m2s)
    },
    do_anova = function(probe_name, factor_name) {
      "Performs anova test for a given `probe_name` and a given `factor_name`."
    	idx_sample = rownames(.self$get_exp_grp())
    	# Dealing with ratio data, reduce data to interesting probes and samples
    	filtred_bp_data = .self$get_data()[probe_name, idx_sample]
      m = aov(filtred_bp_data~.self$get_exp_grp()[,factor_name])
      # tests
      s = shapiro.test(residuals(m))
      shap = -log10(s$p.value)
      f_kc = m$coefficients[2]
      pval = -log10(summary(m)[[1]][["Pr(>F)"]][1])
      ret = list(shap=shap, pval=pval) 
      return(ret)
    },
    plot_boxplot = function(probe_name, factor_name, ylim, las=2, col="grey", border="black", bp_function_name="boxplot", ...) {
      "Draw the box plot for a given `probe_name` and a given `factor_name`."
    	idx_sample = rownames(.self$get_exp_grp())
    	# Dealing with ratio data, reduce data to interesting probes and samples
    	filtred_bp_data = .self$get_data()[probe_name, idx_sample]
    	# Box plots
      if (missing(ylim)) {
    	  ylim = c(min(filtred_bp_data), max(filtred_bp_data))
      }
	    gene_name = "gene_name"
	    # boxplot_filename <- paste(study_dirname, "/", gene, "_", probe_name, "_", factor_name, ".pdf", sep="")
	    # pdf(file=boxplot_filename, height=10, width=10)# open jpeg device with specified dimensions
      bp_function = get(bp_function_name)
	    bp_function(filtred_bp_data~.self$get_exp_grp()[,factor_name], # open boxplot
          # what=c(1,1,1,0),
	      ylim = ylim,                     # fixed vertical scale
	      col = col, border = border,  # colors of inside and border of box
	      las = las,                         # written vertically
	      xlab = factor_name,
	      main = paste(gene_name, probe_name, sep="@"), # title
        ...
	    )
	    # dev.off()
    }  
  )
)




