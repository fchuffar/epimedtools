#' Update vector of genes using epimed API
#'
#' This function ppdates vector of genes using epimed API
#' @param genes_to_update a vector of genes to be updated
#' @param jobid the id of the job if job exists
#' @param taxid the name of the TCGA project to retrieve (9606 for homo spaiens, 10090 for mus musculus)
#' @param WAIT a bolean indicating if call to service is synchrone or assynchrone
#' @importFrom httr POST
#' @importFrom jsonlite fromJSON
#' @export
epimed_api_update = function(genes_to_update, jobid, taxid=9606, WAIT=TRUE) {
  if (missing(jobid)) {
    url = "http://epimed.univ-grenoble-alpes.fr/database/query/jobid"
    jobid = jsonlite::fromJSON(url)    
    print(paste0("create job ", jobid))
    url = "http://epimed.univ-grenoble-alpes.fr/database/query/genes/update"
    body = list(jobid=jobid, symbols=paste(genes_to_update, collapse=", "), taxid=taxid)
    response = httr::POST(url, body = body, encode = "form")
  }

  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobstatus?jobid=", jobid)
  job = jsonlite::fromJSON(url)
  while (job$status != "success" & WAIT) {
    print(paste0("job:", job$jobid, " ", job$status, " ", job$current, "/", job$total))
    system("sleep 2")
    job = jsonlite::fromJSON(url)
  }
  print(paste0("job:", job$jobid, " ", job$status, " ", job$current, "/", job$total))

  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobs?jobid=", jobid)
  foo = read.csv2(url, header=TRUE, sep=";", stringsAsFactors=FALSE)
  return(foo)
}



#' Retrieve exp_grp of a given TCGA project from epimeddb
#'
#' This function retrieves exp_grp of a given TCGA project from epimeddb
#' @param tcga_project the name of the TCGA project to retrieve
#' @export
get_tcga_exp_grp = function(tcga_project) {
  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/parameters/",tcga_project)
  df1 = read.csv2(url, header=TRUE, sep=";", stringsAsFactors=FALSE, dec=".", na.strings="")
  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/expgroup/",tcga_project)
  df2 = read.csv2(url, header=TRUE, sep=";", stringsAsFactors=FALSE, dec=".", na.strings="")
  df = merge(df2,df1,by=1)
  # identify rownames
  df$rn_id = substr(df[,1],1,15)
  # dup rn
  dup_rn_id = df$rn_id[duplicated(df$rn_id)]
  dup_rn_id
  foo = df[df$rn_id %in% dup_rn_id,]
  idx = !is.na(foo[1,] == foo[2,]) & foo[1,] != foo[2,]
  n = colnames(idx)[idx]
  foo[,n]
  # remove FFPE
  df$is_ffpe = as.logical(df$is_ffpe)
  sum(df$is_ffpe)
  df = df[!df$is_ffpe,]
  # dup rn
  dup_rn_id = df$rn_id[duplicated(df$rn_id)]
  dup_rn_id
  # # dedup rn
  # df = df[!duplicated(df$rn_id),]
  # df$rn_id = NULL
  rownames(df) = substr(df[,1],1,15)

  # exp_grp
  exp_grp = df
  exp_grp$dead = as.logical(exp_grp$dead)
  exp_grp$os = survival::Surv(exp_grp$os_month, exp_grp$dead)
  return(exp_grp)
}

#' gdc-client wrapperChecking data study properties
#'
#' This function wraps gdc-client.
#' @param manifest_filename the path to the manifest file
#' @param dest_dir the directory to put data
#' @export
gdc_client_wrapper = function(manifest_filename, dest_dir, files) {
  manifest = read.table(manifest_filename, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (!missing(files)) {
    head(manifest)
    # sum(!files %in% manifest$filename)
    idx = manifest$filename %in% files
    manifest = manifest[idx,]
  }
  dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)
  # Go
  AGAIN = TRUE
  while (AGAIN) {
    foo = epimedtools::monitored_apply(mod=10, manifest, 1, function(l) {
      # l = manifest[1,]
      full_file_name = paste0(dest_dir, "/", l[["id"]], "/", l[["filename"]])
      if (file.exists(full_file_name)) {
        # full_file_name = "../data/count_data/7da031d1-04f3-407b-af7b-61bf4edd38e7/9f0897a4-db0c-43f8-9d75-8debe9b6c847.htseq.counts.gz"
        # if (tools::md5sum(full_file_name) != l[["md5"]]) {
        #   print(full_file_name)
        #   return(TRUE)
        # } else
          {
            return(FALSE)
        }
      }  else {
        print(full_file_name)
        return(TRUE)
      }
    })
    command="gdc-client"
    if (sum(foo) > 0) {
      write.table(manifest[foo,], "tmp_manifest.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
      args=paste0("download -m tmp_manifest.txt -d ", dest_dir)
      paste(command, args)
      system2(command, args)
    } else {
      AGAIN = FALSE
    }
  }
}

#' Checking data study properties
#'
#' This function checks study properties.
#' @param study The study to be checked.
#' @export
check_study =function (study) {
  if (!is.matrix(study$data)) {
    stop("$data is not a matrix.")
  }
  if (!is.data.frame(study$platform)) {
    stop("$platform is not a data.frame.")
  }
  if (!is.data.frame(study$exp_grp)) {
    stop("$exp_grp is not a data.frame.")
  }
  if (sum(!rownames(study$platform) %in% rownames(study$data))) {
    stop("$platform rownames are not all included in $data rownames.")
  }
  if (sum(!rownames(study$exp_grp) %in% colnames(study$data))) {
    stop("$exp_grp rownames are not all included in $data colnames.")
  }
  if (sum(rownames(study$data) != rownames(study$platform))) {
    warning("$platform rownames are not strictly equal to $data rownames.")
  }
  if (sum(colnames(study$data) != rownames(study$exp_grp))) {
    warning("$exp_grp rownames are not strictly equal to $data colnames.")
  }
  print("study checked.")
}


#' Fuming data
#'
#' This function fums data.
#' @param data           data.frame extracted from json
#' @param colname        column name to flat.
#' @export
tcga_fum_data = function(data, colname) {
  # data_json_file = "data.project-TCGA-BRCA.2018-03-02.json"
  # data = jsonlite::fromJSON(txt=data_json_file, flatten=TRUE)
  data_orig = data
  # filter row with empty col named colname
  len = sapply(data_orig[[colname]], length)
  data = data_orig[len==max(len),]
  # flat col named colname
  foo = data_orig[len==max(len),colname]
  sapply(foo, dim)
  diagnoses = do.call(rbind, data_orig[len==max(len),colname])
  colnames(diagnoses) = paste0(colname, "_", colnames(diagnoses))
  if (dim(diagnoses)[1] == nrow(data)) {
    data = cbind(data, diagnoses)
    data[[colname]] = NULL
  } else {
    stop(paste0("can't fum ", colname, " in data."))
  }
  # insert filtered rows
  for (i in which(len!=max(len))){
    keys = colnames(data)[colnames(data) %in% colnames(data_orig)]
    data[nrow(data)+1,] = NA
    data[nrow(data), keys] = data_orig[i,keys]
  }
  return(data)
}

#' Exploring data
#'
#' This function explores data
#' @param df           data.frame extracted from json
#' @param lev        level of explorer, increased when recursiv call
#' @export
tcga_explore_data = function(df, lev=0) {
  unlist(sapply(names(df), function(n) {
    # print(n)
    l = length(df[[1,n]])
    if (l > 1) {
      print(paste(paste(rep(" ", lev*2), collapse=""), ">>>", n))
      tcga_explore_data(df[[1,n]], lev = lev+1)
      return(n)
    } else {
      print(paste(paste(rep(" ", lev*2), collapse=""), n))
      return(NULL)
    }
  }))[1]
}

#' Merging metadata cart and clinicals data
#'
#' This function merge metadata cart and clinicals data
#' @param metadata_json_file        json filename for metadata cart
#' @param clinical_json_file        json filename for for clinical data
#' @param col_to_remove             names of columns to be removed
#' @export
#' @importFrom utils read.table
#' @importFrom utils head
tcga_merge_metadata_cart_and_clinicals = function(metadata_json_file, clinical_json_file, col_to_remove) {
  # clinical
  # clinical_json_file = "clinical.project-TCGA-BRCA.2018-03-02.json"
  clinical = jsonlite::fromJSON(txt=clinical_json_file, flatten=TRUE)
  print("Exploring clinical")
  foo = tcga_explore_data(clinical[1,])
  while (!is.null(foo)) {
    # metadata_cart = tcga_fum_data(metadata_cart, "annotations")
    clinical = tcga_fum_data(clinical, foo)
    foo = tcga_explore_data(clinical[1,])
  }
  head(clinical)
  rownames(clinical) = clinical$case_id
  dim(clinical)

  # metadata_cart
  metadata_cart = jsonlite::fromJSON(txt=metadata_json_file, flatten=TRUE)
  print("Exploring metadata_cart")
  foo = tcga_explore_data(metadata_cart[1,])
  if (!missing(col_to_remove)) {
    print(paste0("removing ", col_to_remove))
    metadata_cart = metadata_cart[,-which(colnames(metadata_cart) %in% col_to_remove)]
    foo = tcga_explore_data(metadata_cart[1,])
  }
  while (!is.null(foo)) {
    # metadata_cart = tcga_fum_data(metadata_cart, "annotations")
    metadata_cart = tcga_fum_data(metadata_cart, foo)
    foo = tcga_explore_data(metadata_cart[1,])
  }
  head(metadata_cart)
  rownames(metadata_cart) = metadata_cart$file_name
  dim(metadata_cart)

  # add clinical to metadata_cart
  foo = clinical[metadata_cart$associated_entities_case_id,]
  sum(colnames(foo) %in% colnames(metadata_cart))
  metadata_cart = cbind(metadata_cart, foo)
  dim(metadata_cart)

  # # which king of data ?
  # metadata_cart$diagnoses_site_of_resection_or_biopsy
  # metadata_cart$diagnoses_days_to_death
  # metadata_cart$diagnoses_tumor_stage
  # metadata_cart$diagnoses_days_to_death
  # metadata_cart$diagnoses_days_to_last_follow_up # fut
  # metadata_cart$diagnoses_vital_status #dead

  # deal with metadata_cart$associated_entities_entity_submitter_id
  metadata_cart$project     = substr(metadata_cart$associated_entities_entity_submitter_id, 1, 4)
  metadata_cart$tss         = substr(metadata_cart$associated_entities_entity_submitter_id, 6, 7)
  metadata_cart$participant = substr(metadata_cart$associated_entities_entity_submitter_id, 9, 12)
  metadata_cart$sample      = substr(metadata_cart$associated_entities_entity_submitter_id, 14, 15)
  metadata_cart$vial        = substr(metadata_cart$associated_entities_entity_submitter_id, 16, 16)
  metadata_cart$portion     = substr(metadata_cart$associated_entities_entity_submitter_id, 18, 19)
  metadata_cart$analyte     = substr(metadata_cart$associated_entities_entity_submitter_id, 20, 20)
  metadata_cart$plate       = substr(metadata_cart$associated_entities_entity_submitter_id, 22, 25)
  metadata_cart$center      = substr(metadata_cart$associated_entities_entity_submitter_id, 27, 28)
  head(metadata_cart)
  unique(metadata_cart$sample)

  # tissue_status
  metadata_cart$tissue_status = NA
  metadata_cart[metadata_cart$sample %in% c("11"),]$tissue_status = "normal"
  metadata_cart[metadata_cart$sample %in% c("01", "06"),]$tissue_status = "tumoral"

  # annotations
  if (!is.null(metadata_cart$annotations)) {
    idx =  sapply(metadata_cart$annotations, is.null)
    sum(!idx)
    metadata_cart = metadata_cart[idx,]
  }

  # survival
  # os_month
  metadata_cart$os_month = sapply(metadata_cart$diagnoses_days_to_last_follow_up, function(fut){
    if (is.na(fut)) {
      return(NA)
    } else {
      fut / 30
    }
  })
  # month_to_death
  metadata_cart$month_to_death = sapply(metadata_cart$diagnoses_days_to_death, function(days_to_death){
    if (is.na(days_to_death)) {
      return(NA)
    } else {
      days_to_death / 30
    }
  })
  # dead
  metadata_cart$dead = sapply(metadata_cart$diagnoses_vital_status, function(da){
    if (is.na(da)) {
      return(NA)
    }
    if (da == "dead") {
      return(TRUE)
    } else if (da == "alive") {
      return(FALSE)
    } else {
      return(NA)
    }
  })
  # replace os_month by month_to_death when month_to_death is not NA
  metadata_cart[!is.na(metadata_cart$month_to_death),]$os_month = metadata_cart[!is.na(metadata_cart$month_to_death),]$month_to_death
  # survival in day
  metadata_cart$os_day = sapply(metadata_cart$diagnoses_days_to_last_follow_up, function(fut){
    if (is.na(fut)) {
      return(NA)
    } else {
      fut
    }
  })
  # month_to_death
  metadata_cart$day_to_death = sapply(metadata_cart$diagnoses_days_to_death, function(days_to_death){
    if (is.na(days_to_death)) {
      return(NA)
    } else {
      days_to_death
    }
  })
  # replace os_day by day_to_death when day_to_death is not NA
  metadata_cart[!is.na(metadata_cart$day_to_death),]$os_day = metadata_cart[!is.na(metadata_cart$day_to_death),]$day_to_death
  # check
  foo = metadata_cart[,c("diagnoses_days_to_last_follow_up", "diagnoses_days_to_death", "os_month", "month_to_death", "diagnoses_vital_status", "dead")]
  # foo = metadata_cart[,c("diagnoses_days_to_last_follow_up", "diagnoses_vital_status")]
  rownames(foo) = colnames(foo) = NULL
  foo[is.na(foo[,1]),]
  foo[is.na(foo[,2]),]
  # survival to NA for NTL
  metadata_cart[metadata_cart$case_id %in% metadata_cart$case_id[duplicated(metadata_cart$case_id)][1],]
  metadata_cart[metadata_cart$case_id %in% metadata_cart$case_id[duplicated(metadata_cart$case_id)][2],]
  metadata_cart[metadata_cart$sample == "11",]$os_month = NA
  metadata_cart[metadata_cart$sample == "11",]$dead = NA
  # # sensoring to 120 month
  # idx = !is.na(metadata_cart$os_month) & metadata_cart$os_month > 120
  # metadata_cart[idx,]$os_month = 120
  # metadata_cart[idx,]$dead = FALSE
  # build survival
  metadata_cart$os = survival::Surv(metadata_cart$os_month, metadata_cart$dead)
  metadata_cart$osd = survival::Surv(metadata_cart$os_day, metadata_cart$dead)
  return(metadata_cart)
}




#' A Function That Perform Intersection Between Sets.
#'
#'
#' This function use recurssively intersect function
#' @param x set: vectors (of the same mode) containing a sequence of items (conceptually) with no duplicated values.
#' @param y set: vectors (of the same mode) containing a sequence of items (conceptually) with no duplicated values.
#' @param ... set: vectors (of the same mode) containing a sequence of items (conceptually) with no duplicated values.
#' @export
intersect_rec = function(x, y, ...){
  if (missing(...)) intersect(x, y)
  else intersect(x, intersect_rec(y, ...))
}

#' A Function That Computes the Reverse Complement of a DNA Sequence
#'
#'
#' This function returns the Reverse complement of a DNA sequence.
#' @param x a string specifying the DNA sequence
#' @importFrom seqinr comp
#' @export
revcomp=function(x){
  toupper(paste(rev(seqinr::comp(strsplit(x, "")[[1]])), collapse=""))
}

#' A Function That Computes the Number of Line a a Text File
#'
#'
#' This function returns the number of line a a text fileplots heatmaps.
#' @param x a string specifying the text file
#' @export
get_nbline = function(x){
  as.integer(system2("wc",
    args = c("-l",
    x,
    " | awk '{print $1}'"),
    stdout = TRUE))
}


#' A Function That Plots Heatmaps.
#'
#'
#' This function plots heatmaps.
#' @param study An epimedtools Study RC object
#' @param pf_col A platform column name use to label gene axis
#' @param exp_grp_col_label An exp_grp column name use to label sample
#' @param exp_grp_idx A vector indexing samples
#' @param platform_idx A vector indexing genes
#' @param method_hclust A string specifying hclust method to use
#' @param nb_clust An interger specifying the number of cluster
#' @param var A string defining the distance used to cluster data
#' @param PLOT_HM_raw A boolean specifying if heatmap needs to be ploted
#' @param USE_CLUST A boolean specifying if dendogram needs to be ploted
#' @param CLUST_COL A boolean specifying if column dendogram needs to be ploted
#' @param CLUST_ROW A boolean specifying if row dendogram needs to be ploted
#' @param BREAK_TOO_EXPENSIVE A boolean specifying if too expensive operation must be breaked.
#' @param ... Parameters passed to gplots::heatmap.2 function
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @importFrom survival survfit
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
#' @importFrom gplots heatmap.2
#' @export
plot_hm2 = function(study, pf_col, exp_grp_col_label, exp_grp_idx, platform_idx, method_hclust="complete", nb_clust=2, var="raw", PLOT_HM_raw=TRUE, USE_CLUST=FALSE, CLUST_ROW=FALSE, CLUST_COL=FALSE, BREAK_TOO_EXPENSIVE=TRUE, colors=c("royalblue", "springgreen", "yellow", "red"), ...) {
  if (missing(pf_col)) {
    pf_col = colnames(study$platform)[1]
  }
  if (missing(exp_grp_col_label)) {
    exp_grp_col_label = colnames(study$exp_grp)[1]
  }
  if (missing(exp_grp_idx)) {
    exp_grp_idx = rownames(study$exp_grp)
  }
  if (missing(platform_idx)) {
    platform_idx = rownames(study$platform[!is.na(study$platform[[pf_col]]), ])
    platform_idx = platform_idx[order(study$platform[platform_idx, pf_col])]
  }
  if (var=="samples") {
    d = study$data[platform_idx, exp_grp_idx]
    if (sum(sapply(apply(d, 2, unique), length) == 1) > 0) {
      d = jitter(d)
    }
    if (dim(d)[2] > 1000 & BREAK_TOO_EXPENSIVE) {
      warning(paste0("Can't correlate, matrix is too big: ", dim(d)[1], "x", dim(d)[2], "."))
      return(NULL)
    }
    d = cor(d, method="sp")
    hc = hclust(dist(1 - d), method=method_hclust)
    if (USE_CLUST) {
      Colv = Rowv = as.dendrogram(hc)
      dendrogram="both"
    } else {
      Colv = Rowv = FALSE
      dendrogram="none"
    }
    rownames(d) = study$exp_grp[exp_grp_idx,exp_grp_col_label]
  } else if (var=="genome") {
    # dat = study$data[platform_idx, exp_grp_idx]
    # idx = apply(dat, 1, function(l){
    #   length(unique(l)) > 1
    # })
    # d = dat[idx,]
    d = study$data[platform_idx, exp_grp_idx]
    if (sum(sapply(apply(d, 2, unique), length) == 1) > 0) {
      d = jitter(d)
    }
    d = t(d)
    if (dim(d)[2] > 1000 & BREAK_TOO_EXPENSIVE) {
      warning(paste0("Can't correlate, matrix is too big: ", dim(d)[1], "x", dim(d)[2], "."))
      return(NULL)
    }
    d = cor(d, method="sp")
    # sum(is.na(d))
    # apply(!is.na(d), 1, sum) == 0
    # length(d)
    hc = hclust(dist(1 - d), method=method_hclust)
    if (USE_CLUST) {
      Colv = Rowv = as.dendrogram(hc)
      dendrogram="both"
    } else {
      Colv = Rowv = FALSE
      dendrogram="none"
      colnames(d) = NULL
      colnames(d)[!duplicated(study$platform[platform_idx,pf_col])] = study$platform[platform_idx,pf_col][!duplicated(study$platform[platform_idx,pf_col])]
    }
  } else {
    d = t(study$data[platform_idx, exp_grp_idx])
    colnames(d) = NULL
    colnames(d)[!duplicated(study$platform[platform_idx,pf_col])] = study$platform[platform_idx,pf_col][!duplicated(study$platform[platform_idx,pf_col])]
    hc = hclust(dist(d), method=method_hclust)
    if (USE_CLUST) {
      Rowv = as.dendrogram(hc)
      Colv = FALSE
      dendrogram="row"
    } else {
      Colv = Rowv = FALSE
      dendrogram="none"
    }
    if (CLUST_COL) {
      tmp_d = study$data[platform_idx, exp_grp_idx]
      if (sum(sapply(apply(tmp_d, 2, unique), length) == 1) > 0) {
        tmp_d = jitter(tmp_d)
      }
      tmp_d = t(tmp_d)
      tmp_d = cor(tmp_d, method="sp")
      hc_row = hclust(dist(1 - tmp_d), method=method_hclust)
      Colv = as.dendrogram(hc_row)
      dendrogram="col"
    }
    rownames(d) = study$exp_grp[rownames(d), exp_grp_col_label]
  }
  grps = cutree(hc, k = nb_clust)

  if (PLOT_HM_raw) {
    patientcolors = rainbow(length(unique(grps)))[grps]
    q = sort(unique(quantile(d, probs=seq(0,1,length.out=20))))
    q = q[ abs(diff(q)) > (max(q) - min(q)) / length(q) / 500]
    cols = colorRampPalette(colors)(length(q)-1)
    if ((dim(d)[1] > 1000 | dim(d)[2] > 1000) & BREAK_TOO_EXPENSIVE) {
      warning(paste0("Can't draw heatmap, matrix is too big: ", dim(d)[1], "x", dim(d)[2], "."))
      return(NULL)
    }
    hm = gplots::heatmap.2(d, dendrogram=dendrogram, Rowv=Rowv, Colv=Colv, scale='none', trace="none", density.info="density", margin=c(5,10), col=cols, breaks=q, RowSideColors=patientcolors, symkey=FALSE, main=paste(dim(d), collapse="x"), ...)
  } else {
    plot(as.dendrogram(hc))
  }
  return(list(grps=grps,hm=hm,hc=hc))
}


#' A Function That Overviews an Experiment Grouping
#'
#' This function overviews an experiment grouping.
#' @param e An expriment grouping
#' @param FULL A boolean to specifying if homogenous or fully heterogenous colums need to be reported
#' @return Nothing
#' @export
overview = function(e, FULL=FALSE) {
  invisible(sapply(colnames(e), function(cn) {
    col = e[,cn]
    if (length(unique(col)) > 1 & length(unique(col)) < length(col) & !FULL) {
      print(paste("============= ", cn, " (", length(unique(col)), ") =============", sep=""))
      if (length(unique(col)) < 100) {
        print(unique(col))
      }
    }
  }))
}

# #' A Function That Requests Epimed Database
# #'
# #' This function requests epimed database.
# #' @param query A character string indicating to SQL query to execute
# #' @param dbname A character string integer the nane of the database
# #' @param host A character string the host of the database
# #' @param port An integer indicating the port for the connexion to the database
# #' @param user A character string indicating the user tuse for the connexion to the database
# #' @param password A character string indicating the password to use for the connexion to the database
# #' @return a dataframe corresponding to data requested
# #' @importFrom DBI dbDriver
# #' @importFrom RPostgreSQL dbConnect
# #' @importFrom RPostgreSQL dbGetQuery
# #' @importFrom RPostgreSQL dbDisconnect
# #' @importFrom RPostgreSQL dbUnloadDriver
# #' @export
# req = function(query, dbname="epimed_prod", host="epimed-db.imag.fr", port=5432, user="epimedtools", password="epimedtools") {
#  # load the PostgreSQL driver
#  drv = dbDriver("PostgreSQL")
#  # create a connection to the postgres database
#  con = dbConnect(drv, dbname=dbname, host=host, port=port, user=user, password=password)
#  # load data
#  info = dbGetQuery(con, query)
#  # close the connection
#  dbDisconnect(con)
#  dbUnloadDriver(drv)
#  return(info)
# }


#' A Function That Draws a Volcano Plot
#'
#' This function draws volcano plot, logratio vs. p-value.
#' @param res_key_lr the column name of the result dataframe to considere as logratio
#' @param res_key_pval the column name of the result dataframe to considere as logratio
#' @param anova_mw_res A dataframe with ANOVA and Mann-Whitney results
#' @param fc_thres A vector of integer corresponding to the foldchange thresholds
#' @param pval_thres An integer corresponding to the p-value thresholds
#' @param legend_place A character string to specify where to put the legend
#' @param LEGEND A boolean to specifying if legend need to be plotted
#' @param ... Parameters passed to plot function
#' @return NULL
#' @importFrom gtools foldchange2logratio
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @export
plot_volcano = function(res_key_lr, res_key_pval, anova_mw_res, fc_thres=c(-1.5, -1.2, 1.2, 1.5), pval_thres=0.05, legend_place="bottomright", LEGEND=TRUE, ...){
  lr_thres = foldchange2logratio(fc_thres)
  plot(anova_mw_res[[res_key_lr]], -log10(anova_mw_res[[res_key_pval]]), xlab=res_key_lr , ylab=paste("-log10(", res_key_pval,")", sep=""), xlim=range(c(lr_thres, anova_mw_res[[res_key_lr]])), ...)
  abline(h=-log10(pval_thres), v=lr_thres, col=2, lty=3)
  if (LEGEND) {
    legend(legend_place, legend=paste("fc thres (", paste(fc_thres, collapse=","), ")", sep=""), col=2, lty=3, cex=0.5)
  }
}

#' A Function That Plots Heatmaps of Differentialy Expressed Genes Accros Conditions
#'
#' This function plots heatmaps of differentialy expressed genes across conditions.
#' @param exp_grp_key the column name of the exprimental grouping to considere
#' @param case_fctr The label of the group that we considere as differentialy expressed
#' @param anova_mw_res A dataframe with ANOVA and Mann-Whitney results
#' @param study An Reference classes object that contains data and metadata
#' @param ctrl_fctr The label of the group that we considere as the reference
#' @param main A character string to explicit the title of the plot
#' @param mw_pval the the suffix the column name corresponding to the Mann-Whitney p-value used to
#' @param gene_pf_colname A character string specifying the name of the platform column to use for the gene name
#' @param PLOT_GS A boolean set to TRUE if genes names to be ploted
#' @param PLOT_PVAL A boolean set to TRUE if ANOVA and Mann-Whitney need to be ploted
#' @param fc_thres A vector of integer coprresponding to the foldchange thresholds
#' @param key A character string that will suffix the column names of the resulting data frame
#' @return probes used to plot the heatmap
#' @importFrom gtools foldchange2logratio
#' @importFrom graphics hist
#' @importFrom graphics image
#' @importFrom graphics points
#' @export
plot_hm = function(exp_grp_key, case_fctr, anova_mw_res, study, ctrl_fctr, main,  mw_pval="mw_pval_adj", gene_pf_colname="gene_name", key, PLOT_GS=FALSE, fc_thres = c(-1.5, -1.2, 1.2, 1.5), PLOT_PVAL=TRUE) {
  lr_thres = gtools::foldchange2logratio(fc_thres)
  if (missing(main)) {main=""}
  if (missing(ctrl_fctr)) {
    ctrl_fctr = "others"
    key_lr_an = paste("logratio_", exp_grp_key, "_", case_fctr, sep="")
    key_fc_an = paste("foldchange_", exp_grp_key, "_", case_fctr, sep="")
    key_pval_an = paste("adj_pval_", exp_grp_key, sep="")
    key_lr_mw = key_lr_an
    key_fc_mw = key_fc_an
    key_pval_mw = key_pval_an
    main = paste(main, " ", case_fctr, sep="")
  } else {
    mw_key = paste(ctrl_fctr, case_fctr, sep=".")
    key_lr_an = paste("logratio_", exp_grp_key, "_", mw_key, sep="")
    key_lr_mw = paste(mw_key, ".logratio", sep="")
    key_fc_an = paste("foldchange_", exp_grp_key, "_", mw_key, sep="")
    key_fc_mw = paste(mw_key, ".foldchange", sep="")
    key_pval_an = paste("adj_pval_", exp_grp_key, sep="")
    key_pval_mw = paste(mw_key, ".", mw_pval, sep="")
    anova_mw_res[[key_lr_an]] = anova_mw_res[[paste("logratio_", exp_grp_key, "_", case_fctr, sep="")]] - anova_mw_res[[paste("logratio_", exp_grp_key, "_", ctrl_fctr, sep="")]]
    anova_mw_res[[key_fc_an]] = gtools::logratio2foldchange(anova_mw_res[[key_lr_an]])
    main = paste(main, " ", mw_key, sep="")
  }

  if (!missing(key)) {
    key_pval_an = paste(key_pval_an, key, sep="_")
    key_pval_mw = paste(key_pval_mw, key, sep="_")
  }

  ## which genes?
  # foo = anova_mw_res[anova_mw_res[[key_pval_an]] < 0.05,]
  # foo = foo[order(foo[[key_pval_an]]), ]
  foo = anova_mw_res[anova_mw_res[[key_pval_mw]] < 0.05,]
  foo = foo[order(foo[[key_pval_mw]]), ]
  idx_1 = rownames(foo)[which(foo[[key_fc_mw]] > 1.5  & foo[[key_fc_an]] > 1.5)  ]
  idx_2 = rownames(foo)[which(foo[[key_fc_mw]] < -1.5 & foo[[key_fc_an]] < -1.5) ]
  # idx_probes = c(idx_1)
  idx_probes = c(idx_2,idx_1)

  # wich samples?
  fact_vals = na.omit(unique(study$exp_grp[[exp_grp_key]]))
  idx_samples = rownames(study$exp_grp)[study$exp_grp[[exp_grp_key]] %in% fact_vals]
  # idx_samples = rownames(study$exp_grp)

  ## data logratio?
  # data_logratio = apply(study$data[idx_probes,], 1, function(l) {
  #   # log2(2^l / mean(2^l))
  #   l - mean(l)
  # })[idx_samples,]

  tmp_data = study$data[idx_probes,idx_samples]
  if (nrow(t(tmp_data)) == 1) {
    tmp_data=t(tmp_data)
    rownames(tmp_data) = idx_probes
  }
  data_logratio = apply(t(tmp_data), 2, function(l) {
    # 2^l / mean(2^l)
    # log2(2^l / mean(2^l))
    l - mean(l)
  })
  if (nrow(t(data_logratio)) == 1) {
    data_logratio=t(data_logratio)
    rownames(data_logratio) = idx_probes
  }

  data_logratio_mean = apply(t(data_logratio[,idx_1]), 2, mean)
  mean_idx = names(data_logratio_mean)[rev(order(data_logratio_mean))]
  # mean_idx = names(data_logratio_mean)
  ordered_samples = sapply(fact_vals, function(f){rownames(study$exp_grp[mean_idx, ])[study$exp_grp[mean_idx, exp_grp_key] == f]})
  idx_samples = unlist(ordered_samples)
  data_logratio = data_logratio[idx_samples,]


  # Go!
  main = paste(main, " (", ncol(data_logratio), ")", sep="")
  hm_carpet = data_logratio
  q = max(abs(quantile(data_logratio, probs=c(0.01, 0.99))))
  # q = max(abs(quantile(data_logratio, probs=c(0.1, 0.9))))

  hm_breaks = c(-max(abs(data_logratio)), seq(-q, q, length.out=49), max(abs(data_logratio)))
  nr = ncol(hm_carpet)
  nc = nrow(hm_carpet)
  col=c("green", "black", "red")
  col_ramp = colorRampPalette(col)(50)

  # plotting...
  # layout(matrix(1:6, 2, byrow=TRUE), respect=TRUE)
  layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)

  if (PLOT_PVAL) {
    plot(-log10(anova_mw_res[[key_pval_an]]), -log10(anova_mw_res[[key_pval_mw]]), pch=".", col=(rownames(anova_mw_res) %in% idx_probes+1), xlab=paste("-log10(",key_pval_an,")", sep=""), ylab=paste("-log10(",key_pval_mw,")", sep=""))
    points(-log10(anova_mw_res[idx_probes, key_pval_an]), -log10(anova_mw_res[idx_probes, key_pval_mw]), pch=16, col=adjustcolor(2, alpha.f=0.3))
    abline(h=-log10(0.05), v=-log10(0.05), col="grey", lty=2)
  }
  # plot(anova_mw_res[[key_lr_an]], anova_mw_res[[key_lr_mw]], col=(rownames(anova_mw_res) %in% c(idx_1, idx_2)) + 1, xlab=key_lr_an, ylab=key_lr_mw, pch=".")
  # abline(h=lr_thres, v=lr_thres, col="grey", lty=2)

  # range
  hist(data_logratio, breaks=hm_breaks, col=col_ramp)
  abline(v=lr_thres, col=2, lty=2)
  # heatmap
  image(1:nc, 1:nr, hm_carpet, xlim = 0.5 + c(0, nc), ylim = 0.5 +
      c(0, nr), axes = FALSE, xlab = "", ylab="", main=main, col=col_ramp,
      breaks = hm_breaks, useRaster=TRUE
  )
  tab_fact = table(study$exp_grp[idx_samples,exp_grp_key])[fact_vals]
  abline(v=cumsum(tab_fact) + 0.5, h=length(idx_2) + 0.5, col=adjustcolor("white", alpha.f=0.7))
  axis(1, tick=FALSE, at=cumsum(tab_fact) - tab_fact/2 + 0.5, labels=fact_vals)
  tab_idx = c(length(idx_2), length(idx_1))
  axis(4, tick=FALSE, at=cumsum(tab_idx) - tab_idx/2 + 0.5, labels=c(paste(ctrl_fctr, ">>", case_fctr, " (", length(idx_2), ")", sep=""), paste(case_fctr, ">>", ctrl_fctr, " (", length(idx_1), ")", sep="")))
  if (PLOT_GS) {
    axis(2, tick=TRUE, at=1:length(idx_probes), labels=study$platform[idx_probes,gene_pf_colname], las=2, cex.axis=0.6)
  }
  return(idx_probes)
}


#' A Function That Plots Survival Panel
#'
#' This function plots ...
#' @param probe_name probe to consider
#' @param sample_names index of samples to consider
#' @param exp_grp_key the column name of the exprimental grouping to considere
#' @param study An Reference classes object that contains data and metadata
#' @param nb_q An integer specifying the number of quantile to initialize discretized expression data
#' @param anova_mw_res A dataframe with ANOVA and Mann-Whitney results
#' @param gene_pf_colname A character string specifying the name of the platform column to use for the gene name.
#' @param ss_key A character string specifying the experimental grouping column to use for survival
#' @param ylim A vector of 2 numeric specifying the y axis range
#' @param cex.axis interger specifying the size of the axis decorations
#' @return NULL
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline
#' @importFrom graphics layout
#' @importFrom stats density
#' @importFrom utils combn
#' @importFrom beanplot beanplot
#' @export
plot_survival_panel = function(probe_name, sample_names, exp_grp_key, study, nb_q=5, anova_mw_res, gene_pf_colname="lower_gs", ss_key="os",cex.axis=0.7, ylim) {
  if (missing(sample_names)) {
    sample_names = rownames(study$exp_grp[!is.na(study$exp_grp[[ss_key]]),])
  }
  data = study$data[probe_name,sample_names]
  exp_grp = study$exp_grp[sample_names,]
  if (!missing(gene_pf_colname)) {
    gene_name = study$platform[probe_name,][[gene_pf_colname]]
  } else {
    gene_name = ""
  }
  # graphicals
  if (missing(anova_mw_res)) {
    design=exp_grp
    model=as.formula(paste("val~", exp_grp_key))
    anova_mw_res = epimedtools::perform_anova_gen(design=design, model=model, data=data)
    rownames(anova_mw_res) = probe_name
  }
  if (missing(ylim)) {
    ylim=range(data[rownames(exp_grp)[!is.na(exp_grp[[exp_grp_key]])]])
  }



  v = study$data[probe_name, sample_names]
  dv = density(v)
  fact = study$exp_grp[sample_names, exp_grp_key]
  pval_cox = coxres(study$exp_grp[sample_names,ss_key], v)[1]
  # quantiles
  vd_all = discr(v, nb_q)
  b_all = attr(vd_all, "breaks")
  # optimize quantiles
  tbc = 2:(nb_q)
  ibs = unlist(lapply( 1:length(tbc), function(n) {
   as.list(as.data.frame(combn(tbc,n)))
  }), recursive=FALSE)
  pv_logranks = sapply(ibs, function(ib){
   vd_it = discr(v, breaks=b_all[c(1,ib,length(b_all))])
   # scurve(ss, vd, main=paste(gene_name, probe_name, sep="@"))
   coxres(study$exp_grp[sample_names,ss_key], vd_it)[1]
  })
  ib = ibs[[which(pv_logranks == min(pv_logranks))]]
  vd_opt = discr(v, breaks=b_all[c(1,ib,length(b_all))])
  b_opt = attr(vd_opt, "breaks")
  # graphicals
  layout(matrix(c(1,1,1,1,2,3,3,3,4,4,5,5), 3, byrow=TRUE))
  # print(anova_mw_res)
  # print(exp_grp_key)
  # print(paste("adj_pval", exp_grp_key, sep="_"))
  # print(anova_mw_res[probe_name, paste("adj_pval", exp_grp_key, sep="_")])
  if (!missing(anova_mw_res)) {
    anova_pval_leg = paste("p_anova=", signif(anova_mw_res[probe_name, paste("adj_pval", exp_grp_key, sep="_")],3), sep="")
  } else {
    anova_pval_leg = ""
  }
  bp = beanplot(study$data[probe_name,rownames(study$exp_grp)]~study$exp_grp[[exp_grp_key]], las=2, log="", ylim=range(study$data),
    main=paste(gene_name, "@", probe_name, " ", anova_pval_leg, sep=""), ylab="log2(exprs)")
  mw_pval_colnames = colnames(anova_mw_res)[grep("mw_pval_adj", colnames(anova_mw_res))]
  case_names = do.call(rbind, strsplit(mw_pval_colnames, ".", fixed=TRUE))[,2]

  mtc = match(case_names, bp$names)
  rng = rep(range(study$data)[2], length(mtc))
  dup = duplicated(paste(mtc,rng, sep="_"))
  while(sum(dup) > 0) {
    rng = rng - dup/2
    dup = duplicated(paste(mtc,rng, sep="_"))
  }
  duplicated(paste(mtc,rng, sep="_"))
  if (length(mtc) > 0) {
    text(mtc, rng, paste(mw_pval_colnames, "=", signif(anova_mw_res[probe_name,mw_pval_colnames],3), sep=""))
  }
  # sapply(case_names, function(case_name) {
  #   cn = colnames(anova_mw_res)[grep(paste(case_name, ".mw_pval_adj", sep=""), colnames(anova_mw_res))]
  #   text(which(case_name == bp$names),range(study$data)[2], paste(cn, "=", signif(anova_mw_res[probe_name,cn],3), sep=""))
  # })
  plot(dv$y, dv$x, ylim=range(v), xlim=rev(range(dv$y)), type='l', ylab="log2(exprs)")
  abline(h=b_all, lty=1, lwd=1, col=adjustcolor(1, alpha.f=0.1))
  abline(h=b_opt, lty=2, lwd=1)
  beanplot(v~fact, las=2, log="", ylim=range(v), main=paste(gene_name, "@", probe_name, " p_lrk=", signif(as.numeric(pval_cox),3)
    # , " p_anova=", signif(as.numeric(pval_cox),3)
    , sep=""))
  abline(h=b_all, lty=1, lwd=1, col="grey")
  abline(h=b_opt, lty=2, lwd=1)
  scurve(study$exp_grp[sample_names,ss_key], vd_all, main=paste(gene_name, probe_name, sep="@"))
  scurve(study$exp_grp[sample_names,ss_key], vd_opt, main=paste(gene_name, probe_name, sep="@"))
  return(NULL)
}










#' A Function That Plots Expression as Beanplots
#'
#' This function plots ...
#' @param probe_name probe to consider
#' @param sample_names index of samples to consider
#' @param exp_grp_key the column name of the exprimental grouping to considere
#' @param study An Reference classes object that contains data and metadata
#' @param anova_mw_res A dataframe with ANOVA and Mann-Whitney results
#' @param gene_pf_colname A character string specifying the name of the platform column to use for the gene name.
#' @param cex.axis A numeric specifying the size of axis text
#' @param ylim A vector of 2 numeric specifying the y axis range
#' @return NULL
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline
#' @importFrom graphics layout
#' @importFrom stats density
#' @importFrom stats as.formula
#' @importFrom utils combn
#' @importFrom beanplot beanplot
#' @export
plot_bean_expr = function(probe_name, sample_names, exp_grp_key, study, anova_mw_res, gene_pf_colname="lower_gs",cex.axis=0.7, ylim) {
  if (missing(sample_names)) {
    sample_names = rownames(study$exp_grp)
  }
  data = study$data[probe_name,sample_names]
  exp_grp = study$exp_grp[sample_names,]
  gene_name = study$platform[probe_name,][[gene_pf_colname]]
  # graphicals
  if (missing(anova_mw_res)) {
    design=exp_grp
    model=as.formula(paste("val~", exp_grp_key))
    anova_mw_res = epimedtools::perform_anova_gen(design=design, model=model, data=data)
    rownames(anova_mw_res) = probe_name
  }
  if (missing(ylim)) {
    ylim=range(data[rownames(exp_grp)[!is.na(exp_grp[[exp_grp_key]])]])
  }
  anova_pval_leg = paste("p_anova=", signif(anova_mw_res[probe_name, paste("adj_pval", exp_grp_key, sep="_")],3), sep="")
  bp = beanplot(data~exp_grp[[exp_grp_key]], las=2, log="", ylim=ylim,
    main=paste(gene_name, "@", probe_name, " ", exp_grp_key, " ", anova_pval_leg, sep=""), ylab="log2(exprs)",cex.axis=cex.axis)
  mw_pval_colnames = colnames(anova_mw_res)[grep("mw_pval_adj", colnames(anova_mw_res))]
  case_names = do.call(rbind, strsplit(mw_pval_colnames, ".", fixed=TRUE))[,2]

  mtc = match(case_names, bp$names)
  rng = rep(range(data)[2], length(mtc))
  dup = duplicated(paste(mtc,rng, sep="_"))
  while(sum(dup) > 0) {
    rng = rng - dup/2
    dup = duplicated(paste(mtc,rng, sep="_"))
  }
  duplicated(paste(mtc,rng, sep="_"))
  if (length(mtc) > 0) {
    text(mtc, rng, paste(mw_pval_colnames, "=", signif(anova_mw_res[probe_name,mw_pval_colnames],3), sep=""))
  }
  return(NULL)
}








#' A Function That Plots Survival Panel
#'
#' This function plots ...
#' @param probe_name probe to consider
#' @param sample_names index of samples to consider
#' @param study An Reference classes object that contains data and metadata
#' @param nb_q An integer specifying the number of quantile to initialize discretized expression data
#' @param gene_pf_colname A character string specifying the name of the platform column to use for the gene name
#' @param ss_key A character string specifying the experimental grouping column to use for survival
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @param main A character string to explicit the title of the plot
#' @param ... Parameters passed to plot_survival_panel_simple2 function.
#' @return NULL
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline
#' @importFrom graphics layout
#' @importFrom stats density
#' @importFrom utils combn
#' @export
plot_survival_panel_simple = function(probe_name, sample_names, study, nb_q=5, gene_pf_colname="lower_gs", ss_key="os", colors=c("deepskyblue", "black", "red"), main, ...) {
  if (missing(sample_names)) {
    sample_names = rownames(study$exp_grp[!is.na(study$exp_grp[[ss_key]]),])
  }
  if (gene_pf_colname %in% colnames(study$platform)) {
    gene_name = study$platform[probe_name,][[gene_pf_colname]]
  } else {
    gene_name = ""
  }
  if (missing(main)) {
    main = paste(gene_name, "@", probe_name, " (", ss_key, ")", sep="")
  } else {
    main = paste(main, " ", gene_name, "@", probe_name, " (", ss_key, ")", sep="")
  }
  v = study$data[probe_name, sample_names]
  ss = study$exp_grp[sample_names,ss_key]
  ret = plot_survival_panel_simple2(ss, v, nb_q=nb_q, gene_pf_colname=gene_pf_colname, colors=colors, main=main, ...)
  return(ret)
}



#' A Function That Plots Survival Panel
#'
#' This function plots ...
#' @param ss A survival structure such as produced by function Surv of package survival.
#' @param v A (discretized) vector indexed as ss.
#' @param nb_q An integer specifying the number of quantile to initialize discretized expression data
#' @param gene_pf_colname A character string specifying the name of the platform column to use for the gene name
#' @param colors Colors to interpolate; must be a valid argument to col2rgb()
#' @param main A character string to explicit the title of the plot
#' @param thresh A vector of integer used fopr thresholds
#' @param ... Parameters passed to scurve function
#' @param PLOT A boolean set to TRUE if data needs to be plotted.
#' @return NULL
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline
#' @importFrom graphics layout
#' @importFrom stats density
#' @importFrom utils combn
#' @export
plot_survival_panel_simple2 = function(ss, v, nb_q=5, gene_pf_colname="lower_gs", colors=c("deepskyblue", "black", "red"), main, thresh, PLOT=TRUE, ...) {
  if (missing(main)) {
    main=""
  }
  idx = !is.na(ss)
  ss = ss[idx]
  v = v[idx]
  dv = density(v)
  hz_cox = coxres(ss, v)
  pval_cox = hz_cox[1]
  pval_glo = pval_cox
  if (missing(thresh)) {
    # quantiles
    vd_all = discr(v, nb_q)
  } else {
    vd_all = discr(v, nb_q, thresh)
    nb_q = length(unique(vd_all))
  }
  b_all = attr(vd_all, "breaks")
  # optimize quantiles
  tbc = 2:(nb_q)
  if (nb_q>2) {
    ibs = unlist(lapply( 1:length(tbc), function(n) {
     as.list(as.data.frame(combn(tbc,n)))
    }), recursive=FALSE)
    pv_logranks = sapply(ibs, function(ib){
      vd_it = discr(v, breaks=b_all[c(1,ib,length(b_all))])
      if (length(unique(vd_it)) < 2) {
        return(1)
      } else {
        return(coxres(ss, vd_it)[1])
      }
    })
    ib = ibs[[which(pv_logranks == min(pv_logranks))[1]]]
    vd_opt = discr(v, breaks=b_all[c(1,ib,length(b_all))])
    hz_opt = coxres(ss,vd_opt)
    pval_opt = min(pv_logranks)
    vd_opt = discr(v, breaks=b_all[c(1,ib,length(b_all))])
  } else {
    vd_opt = discr(v, breaks=b_all)
    hz_opt = coxres(ss,vd_opt)
    pval_opt = hz_opt[1]
  }
  b_opt = attr(vd_opt, "breaks")
  # graphicals
  if (PLOT) {
    layout(matrix(c(1,1,1,1,2,2,3,3,2,2,3,3), 3, byrow=TRUE), respect=TRUE)
    main2 = paste(main, " p_lrk=", signif(as.numeric(pval_cox),3), sep="")
    plot( dv$x, dv$y, xlim=range(v), ylim=range(dv$y), type='l', xlab="log2(exprs)", ylab="", main=main2)
    abline(v=b_all, lty=1, lwd=1, col=adjustcolor(1, alpha.f=0.1))
    abline(v=b_opt, lty=2, lwd=1)
    # print(vd_all)
    # vd_all <<- vd_all
    scurve(ss, vd_all, main=main, colors=colors, ...)
    scurve(ss, vd_opt, main=main, colors=colors, ...)
  }
  return(list(boundaries=b_opt, pval_glo=pval_glo, pval_opt=pval_opt, hz_opt=hz_opt, card=table(vd_opt)))
}



# exp_grp_key = "histo"
# exp_grp_key = "histo.y"
# nb_q = 5
# probe_name = probe_names[1]
# plot_survival_panel(probe_name, sample_names, exp_grp_key, study)



#' A Function That Computes Survival Table.
#'
#' This function computes the false discovery rate from a vector of independent p-values.
#' It extracts the cutoff corresponding to the specified FDR. See Benjamini & Hochberg 1995 paper.
#' @param probe_names index of probes to consider
#' @param sample_names index of samples to consider
#' @param exp_grp_key the column name of the exprimental grouping to considere
#' @param study An Reference classes object that contains data and metadata
#' @param suffix A character string used to suffix the cache file nam
#' @param USE_CACHE A boolean set to TRUE if cache could be used
#' @param PLOT_SCURVE A boolean set to TRUE if Survival curve of the corresponding `exp_grp_key` needs to be ploted.
#' @return A dataframe taht contains survival infirmations
#' @export
compute_survival_table = function(probe_names, sample_names, exp_grp_key, study, suffix="", USE_CACHE=FALSE, PLOT_SCURVE=FALSE) {
  # ss
  ss = study$exp_grp[sample_names,]$ss
  if (PLOT_SCURVE) {
    cox_results = scurve(ss, study$exp_grp[sample_names, exp_grp_key], main=paste("Survival on ", exp_grp_key, sep=""))
    print(cox_results)
  }
  # all probes
  survival_res_filename = paste("survival_", suffix, "_res.rds", sep="")
  if (!file.exists(survival_res_filename) | !USE_CACHE) {
    survival_res = monitored_apply(t(probe_names), 2, function(probe_name) {
      v = study$data[probe_name, sample_names]
      return(coxres(ss, v))
    })
    survival_res = data.frame(t(survival_res), stringsAsFactors=FALSE)
    rownames(survival_res) = probe_names
    if (USE_CACHE) {
      saveRDS(survival_res, survival_res_filename)
    }
    # survival_res$lower_gs = study$platform[rownames(survival_res), ]$lower_gs
  }
  if (USE_CACHE) {
    survival_res = readRDS(survival_res_filename)
  }
  return(survival_res)
}

# suffix = "lung"
# exp_grp_key = "histo.y"
# sample_names = rownames(study$exp_grp[!is.na(study$exp_grp$fut),])
# probe_names = probe_gene_tab$probe
# survival_lung_res = compute_survival_table(probe_names, sample_names, study, exp_grp_key, suffix, USE_CACHE=FALSE)
# rownames(survival_lung_res[survival_lung_res$pv_logrank < 0.001,])





#' A Function that Plots Pairs Correpaltion
#'
#' This function enhances pairs fucntion by adding decoration in classic pairs graphic.
#'
#' @param ... Parameters passed to pairs function.
#' @importFrom graphics pairs
#' @importFrom graphics par
#' @importFrom graphics strwidth
#' @importFrom graphics text
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @importFrom stats symnum
#' @export

et_pairs = function(...) {
  panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = abs(cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex = 0.8/strwidth(txt)

      test = cor.test(x,y)
      # borrowed from printCoefmat
      Signif = symnum(test$p.value, corr = FALSE, na = FALSE,
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))

      text(0.5, 0.5, txt, cex = cex * r)
      text(.8, .8, Signif, cex=cex, col=2)
  }
  pairs(..., upper.panel=panel.cor)

}


#' A Function that Plots Quality Control for an Expression Data Matrix
#'
#' This function plots a boxplot for each sample in an expression data matrix.
#'
#' @param data A matrix of exrpression values for probe/gene (lines) and confditions (colums).
#' @param USE_LOG2_FOR_EXPR A boolean set to TRUE if we want to plot the log of the expression.
#' @param ... Parameters passed to plot function.
#' @importFrom graphics axis
#' @importFrom graphics boxplot
#' @export
qc_expr = function(data, USE_LOG2_FOR_EXPR=TRUE, ...) {
  if (USE_LOG2_FOR_EXPR) {
    trans= log2
    ylab = "log2(expr)"
  } else {
    trans = identity
    ylab = "expr_ratio"
  }
  plot(0,0, col=0, xlim=c(1,ncol(data)), ylim=range(trans(data)), xaxt="n", xlab="", ylab=ylab, ...)
  axis(1, at=1:ncol(data), labels=colnames(data), tick=TRUE, las=2)
  for (i in 1:ncol(data)) {
    boxplot(trans(data[,i]), at=i, add=TRUE, axes=FALSE)
  }
}





#' A Function Performs Differential Expression Analysis Based On ANOVA
#'
#' This function peforms an ANOVA for each probe/gene (lines) across confditions (colums) of a data matrix.
#' It returns a dataframe that includes many metrics as a results (betas, p-values, ...).
#' Genes are indexed by `g` in `G`,
#' conditions are indexed by `h` in `H` (histology, sex...) and
#' samples are indexed by `i` in `I`.
#' For each gene `g`:
#' \itemize{
#'   \item we compute a linear model of the level of expression value according to histological groups: `log2(exprs_{g, i}) = beta_{g, h} + epsilon_{g, i}`.
#'   \item we perform ANOVA test on this model and harvest corresponding p-value
#'   \item we extract of each group: `beta_{g, h} = frac{sum_{g in G} mean(log2(exprs_{g,i}))}{|G|} - mean(log2(exprs_{g,i}))`
#'   \item we compute corresponding logratio scores: `logratio_{g, h} = (1 + frac{1}/{|H| - 1}) beta_{g, h}`
#'   \item we compute corresponding foldchange scores using `gtools::foldchange2logratio` function
#' }
#' These values are retruned as a dataframe which fields are:
#' \itemize{
#'   \item probename:            the name of the corresponding probe
#'   \item ad_pval:              the result of the Anderson-Darling normality test (hypothesis for ANOVA)
#'   \item intercept:            the mean of means of each group
#'   \item beta_hhh_FFF:         the previously describe beta value for the factor `FFF` of the experimental grouping field `hhh`
#'   \item logratio_hhh_FFF :    the previously describe logratio value for the factor `FFF` of the experimental grouping field `hhh`
#'   \item foldchange_hhh_FFF:   the previously describe foldchange value for the factor `FFF` of the experimental grouping field `hhh`
#'   \item pval_hhh:             the p-value associated to the ANOVA test
#'   \item adj_pval_hhh:         the adjusted p-value associated to the ANOVA test using Benjamini-Hochberg procedure to the threshold of 0.05.
#' }
#' @param design A dataframe that describes the design of the experiment.
#' @param model A model describe the test to apply.
#' @param data A matrix of exrpression values for probe/gene (lines) and confditions (colums).
#' @param key A character string that will suffix the column names of the resulting data frame
#' @param MONITORED_APPLY A boolean set to TRUE if we want to monitor the loop on lines (usefull for debugging).
#' @param AD_TEST A boolean set to TRUE if we want to perform normality test.
#' @importFrom nortest ad.test
#' @importFrom stats residuals
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats p.adjust
#' @importFrom stats anova
#' @importFrom gtools logratio2foldchange
#' @export
perform_anova_gen = function(design, model, data, key=NULL, MONITORED_APPLY=FALSE, AD_TEST=TRUE) {
  options(contrasts=c("contr.sum", "contr.poly"))
  if (nrow(t(data)) == 1) {
    data=t(data)
    data = t(data[,rownames(design)])
  } else {
    data = data[,rownames(design)]
  }
  anova_res = monitored_apply(data, 1, function(l) {
    # print(l)
    d = data.frame(val=l)
    exp_grp_key = as.character(model)[3]
    d[[exp_grp_key]] = as.factor(design[[exp_grp_key]])
    m = lm(model, data=d)
    # Anderson-Darling (normality)
    if (AD_TEST) {
      ad = nortest::ad.test(residuals(m))
      ad_pval = ad$p.value
    } else {
      ad_pval = 1
    }
    aov_coeff = c(intercept=m$coefficients[[1]])
    for (xlevels_name in names(m$xlevels)) {
      for (lev in m$xlevels[[xlevels_name]]) {
        idx = which(lev == m$xlevels[[xlevels_name]])
        if (lev == rev(m$xlevels[[xlevels_name]])[1]) {
          aov_coeff[[paste("beta", xlevels_name, lev, sep="_")]] = -sum(m$coefficients[paste(xlevels_name, 1:(idx - 1), sep="")])
        } else {
          aov_coeff[[paste("beta", xlevels_name, lev, sep="_")]] = m$coefficients[paste(xlevels_name, idx, sep="")]
        }
      }
    }

    anov = anova(m)
    aov_pval = anov[names(m$xlevels),5]
    names(aov_pval) = paste("pval", names(m$xlevels), sep="_")
    # Mann Withney
    wm = wilcox.test(model, data = d)
    mw_pval = wm$p.value
    ret = c(ad_pval=ad_pval, aov_coeff, aov_pval, mw_pval)
    if (!is.null(key)) {
      names(ret) = paste(names(ret), key, sep="_")
    }
    return(ret)
  })
  anova_res = data.frame(t(anova_res))
  # adj_pval
  for (colname in colnames(anova_res)[which(strtrim(colnames(anova_res), 5) == "pval_")]) {
    anova_res[[paste("adj", colname, sep="_")]] = p.adjust(anova_res[[colname]], method="BH")
  }
  anova_res$mw_padj = p.adjust(anova_res$mw_pval) 
  # logratio and foldchange
  beta_colnames = colnames(anova_res)[which(strtrim(colnames(anova_res), 5) == "beta_")]
  nb_lev = length(beta_colnames)
  for (colname in beta_colnames) {
    suffix = substr(colname, 6, 10000)
    anova_res[[paste("logratio", suffix, sep="_")]] = anova_res[[colname]] * (1 + 1 / (nb_lev - 1))
    anova_res[[paste("foldchange", suffix, sep="_")]] = logratio2foldchange(anova_res[[paste("logratio", suffix, sep="_")]])
  }
  return(anova_res)
}


#' A Function Performing ANOVAs
#'
#' This function peforms an ANOVA for each probe/gene (lines) across confditions (colums) of a data matrix. Its return a dataframe that includes many metrics as a results (beta, pval...).
#'
#' @param design A dataframe that describes the design of the experiment.
#' @param model A model describe the test to apply.
#' @param data A matrix of exrpression values for probe/gene (lines) and confditions (colums).
#' @param key A character string that will suffix the column names of the resulting data frame
#' @param correction A numeric vector (typically in [-1, 1]) to correct the beta signs according to the wanted reference.
#' @param MONITORED_APPLY A boolean set to TRUE if we want to monitor the loop on lines (usefull for debugging).
#' @param AD_TEST A boolean set to TRUE if we want to perform normality test.
#' @importFrom nortest ad.test
#' @importFrom stats residuals
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats p.adjust
#' @importFrom stats anova
#' @importFrom gtools logratio2foldchange
#' @export
perform_anova = function(design, model, data, key, correction = 1, MONITORED_APPLY=FALSE, AD_TEST=TRUE) {
  options(contrasts=c("contr.sum", "contr.poly"))
  anova_res = monitored_apply(data[, design$sample], 1, function(l) {
    # print(l)
    design$val = log2(l)
    m = lm(model, data=design)
    # Shapiro-Wilks (normality)
    # sw = shapiro.test(residuals(m))
    # sw_pval = sw$p.value
    # Anderson-Darling (normality)
    if (AD_TEST) {
      ad = ad.test(residuals(m))
      ad_pval = ad$p.value
    } else {
      ad_pval = 1
    }
    # Brown-Forsythe (homoskedasticity)
    # if (length(grep("*", "a*a", fixed=TRUE) > 0)) {
    #
    # }
    # bf = lawstat::levene.test(residuals(m), do.call(paste,design[,as.character(model[c(3:min(length(model), 4))])]))
    # bf_pval = bf$p.value
    # boxplot(val~sp, data=design)
    anov = anova(m)
    # print(anov)
    m1 = aov(model, data=design)
    summ_aov = summary(m1)
    # model.tables(m1)
    aov_coeff = m1$coefficients[-1] * correction
    len = length(aov_coeff)
    aov_pval = summ_aov[[1]][1:len,5]
    if (len > 1) {
      names(aov_coeff) = paste("beta_", key, 1:len, sep="")
      names(aov_pval) = paste("pval_", key, 1:len, sep="")
    } else {
      names(aov_coeff) = paste("beta", key, sep="_")
      names(aov_pval) = paste("pval", key, sep="_")
    }
    # ret = c(ad_pval=ad_pval, bf_pval=bf_pval, aov_coeff, aov_pval)
    # names(ret)[1:2] = paste(names(ret)[1:2], key, sep="_")
    ret = c(ad_pval=ad_pval, aov_coeff, aov_pval)
    names(ret)[1] = paste(names(ret)[1], key, sep="_")
    # ret
    # print(ret)
    return(ret)
  })
  anova_res = data.frame(t(anova_res))
  len = length(grep("beta", names(anova_res)))
  if (len > 1) {
    for (i in 1:len) {
      anova_res[[paste("adj_pval_", key, i, sep="")]] = p.adjust(anova_res[[paste("pval_", key, i, sep="")]], method="BH")
      anova_res[[paste("logratio_", key, i, sep="")]] = 2 * anova_res[[paste("beta_", key, i, sep="")]]
      anova_res[[paste("foldchange_", key, i, sep="")]] = logratio2foldchange(anova_res[[paste("logratio_", key, i, sep="")]])
    }
  } else {
    anova_res[[paste("adj_pval_", key, sep="")]] = p.adjust(anova_res[[paste("pval_", key, sep="")]], method="BH")
    anova_res[[paste("logratio_", key, sep="")]] = 2 * anova_res[[paste("beta_", key, sep="")]]
    anova_res[[paste("foldchange_", key, sep="")]] = logratio2foldchange(anova_res[[paste("logratio_", key, sep="")]])
  }
  return(anova_res)
}


#' A Function That Plots Survival Curves.
#'
#' This function takes a survival structure such as produced by function Surv of
#' package survival. Takes a (discretized) vector v, indexed as ss.
#' Performs the Cox proportional hazard rate
#' test of ss against v. It plots the survival curves for each
#' level of v. It also prints the p-value of the log-likelihood test.
#'
#' @param ss A survival structure such as produced by function Surv of package survival.
#' @param v A (discretized) vector indexed as ss.
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @param main A character string to explicit the title of the plot
#' @param legend A vector of character to explicit the legend of the plot
#' @param nb_sign An integer  indicating the number of significant digits to be used
#' @param legend_place A character string to specify where to put the legend
#' @param PLOT_LEGEND A boolean to specifying if legend need to be plotted
#' @param censoring date at which data are censored
#' @param cex.leg A numeric specifying the size of the legend
#' @param ... Parameters passed to plot function.
#' @importFrom survival survfit
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
#' @export
scurve = function(ss, v, colors=c("deepskyblue", "black", "red"), main="Survival", legend, nb_sign=3, legend_place="topright",PLOT_LEGEND=TRUE , cex.leg=1, censoring, ...) {

  idx = !is.na(ss)
  ss = ss[idx]
  if (missing(v)) {
    v = rep(1, sum(idx))
  } else {
    v = v[idx]
  }

  if (!missing(censoring)) {
    date = ss[,1]
    state = ss[,2]
    # censoring to 120 month
    idx = !is.na(date) & date > censoring
    date[idx] = censoring
    state[idx] = 0
    ss =  survival::Surv(date, state)
  }

  if (length(unique(v)) == 1) {
    cox_results = 1
  } else {
    cox_results = coxres(ss,as.character(v))
  }
  pv = cox_results[1]
  if (pv<1e-100) {
    pvt = "<1e-100"
  } else {
    pvt = format(pv,
    digits=3,scientific=TRUE)
  }
  sf=survfit(ss~v)
  levels = length(na.omit(unique(v)))
  col = colorRampPalette(colors)(levels)
  main= paste(main, " p_lrk=", signif(as.numeric(pvt), nb_sign), sep="")
  plot(sf, col=col, main=main, ...)
  tab = table(v)
  tab=tab[tab!=0]
  if (PLOT_LEGEND) {
    if (missing(legend)) {
      if ("breaks" %in% names(attributes(v))) {
        b = signif(attr(v, "breaks"),3)
        legend = paste("[", b[1:(length(b)-1)], ",", b[2:length(b)], c(rep("[", length(b)-2), "]"), sep="")
      } else {
        legend = names(tab)
      }
    }
    legend = paste(legend, " (", tab, ")", sep="")
    legend(legend_place, legend=legend, col=col, pch=3, lty=1, cex=cex.leg)
  }
  return(cox_results)
}

#' A Function That Fits the Cox Regression Model
#'
#' This function takes a survival structure such as produced by function Surv of
#' package survival. Takes a (discretized) vector v, indexed as ss.
#' Fits the Cox regression model of ss against v, and tests the
#' proportional hazards assumption.
#' Returns as a named vector of length 5:
#'      pv_logrank: the significance p-value (logrank)
#'      pvhz: the p-value for the validity of the model
#'      hrlb: the lower bound of the 95% confidence interval for hr
#'      hr: the hazard ratio
#'      hrub: the upper bound of the 95% confidence interval for hr
#' @param censoring date at which data are censored
#' @param ss A survival structure such as produced by function Surv of package survival.
#' @param v A (discretized) vector indexed as ss.
#' @return A named vector of length 5.
#' @importFrom survival coxph
#' @importFrom survival cox.zph
#' @export
coxres = function(ss,v) {
  f = suppressWarnings(coxph(ss~v))
  sf =  summary(f)
  pv_logrank = sf$sctest[3]
  tf = cox.zph(f)
  tf = tf$table
  pvhz = tf[dim(tf)[1],3]
  hr = sf$conf.int[1,1]
  hrlb = sf$conf.int[1,3]
  hrub = sf$conf.int[1,4]
  res = c(pv_logrank,pvhz,hrlb,hr,hrub)
  names(res) = c("pv_logrank","pvhz","hrlb","hr","hrub")
  return(res)
}

#' A Function Discretize a Vector of Numeric
#'
#' Takes a vector v of numeric. Replaces the values of v by integers
#' in (1:nd), according to values of breaks (increasing vector of numeric).
#' If breaks is missing, nd regular quantiles are taken.
#
#' @param v A vector of numeric.
#' @param nd The number of classes.
#' @param breaks An increasing vector of numeric.
#' @return A vector of factors.
#' @importFrom stats quantile
#' @export
discr = function(v, nd=5, breaks){
  if (missing(breaks)) {
    b = quantile(v,(0:nd)/nd)
  } else {
    b = c(min(c(v,breaks)), breaks, max(c(v, breaks)))
  }
  b = sort(unique(b))
  vd = cut(v, b, include.lowest=TRUE, right=FALSE, labels=(1:(length(b)-1)))
  # rnd = floor(log10(abs(b[-length(b)]-b[-1])))
  # print(rnd)
  labels = paste("[", signif(b[-length(b)], 3), ", ", signif(b[-1],3), "[", sep="")
  labels[length(labels)] = paste(substr(labels[length(labels)], 1, nchar(labels[length(labels)],)-1), "]", sep="")
  attributes(vd)$levels = labels
  attributes(vd)$class = c("factor", "ordered")
  return(structure(vd, breaks=b))
}

#' A Function That Returns the Longest Common Prefix.
#'
#' This function returns the longest common prefix.
#'
#' @param s A vector of character.
#' @param split character passed to strsplit function as split.
#' @param fixed boolean passed to gsub function as fixed.
#' @importFrom stats na.omit
#' @export
longest_common_prefix = function(s, split="", fixed=FALSE) {
  s = na.omit(s)
  l = min(sapply(strsplit(unique(s), split), length))
  f = lapply(strsplit(s, split), function(x) {x[1:l]})
  m = do.call(rbind, f)
  vl = apply(m, 2, function(col) {
    length(unique(col))
  })
  d_idx = which(vl!=1)
  if (length(d_idx) == 0) {
    d_idx = length(vl)
  }
  fd_idx = d_idx[1] - 1
  if (fd_idx > 0) {
    lcp = paste(c(f[[1]][1:fd_idx], ""), collapse=split)
  } else {
    lcp = ""
  }
  return(lcp)
}

#' A Function That Summarizes Expriment Grouping.
#'
#' This function summarizes expriment grouping.
#'
#' @param exp_grp a dataframe that describes the experimental grouping.
#' @export
summarize_exp_grp = function(exp_grp) {
  foo = sapply(colnames(exp_grp), function(cn) {
    nb_fact = length(unique(exp_grp[[cn]]))
    if (nb_fact > 1 & nb_fact < length(exp_grp[[cn]])) {
      print(paste("____________________________", cn, "____________________________"))
      print(table(exp_grp[[cn]]))
    }
  })
}
#' A Function That Simplifies Expriment Grouping.
#'
#' This function simplyfies expriment grouping.
#'
#' @param exp_grp a dataframe that describes the experimental grouping.
#' @return a simplified experimental grouping.
#' @export
simplify_exp_grp = function(exp_grp) {
  colnames(exp_grp) = simplify_column_names(colnames(exp_grp))
  col_div = apply(exp_grp, 2, function(col) {
    length(unique(col))
  })
  exp_grp = exp_grp[,col_div != 1]
  for (n in colnames(exp_grp)) {
    exp_grp[[n]] = simplify_column_names(exp_grp[[n]])

    exp_grp[[n]] = simplify_factor_names(exp_grp[[n]])
  }
  return(exp_grp)
}

#' A Function That Computes de False Discovery Rate.
#'
#' This function computes the false discovery rate from a vector of independent p-values.
#' It extracts the cutoff corresponding to the specified FDR. See Benjamini & Hochberg 1995 paper.
#' @param x A vector of independent p-values.
#' @param FDR The false discovery rate.
#' @return The corresponding cutoff.
#' @importFrom stats na.omit
#' @export
FDR = function (x, FDR)
{
    x <- sort(na.omit(x))
    N = length(x)
    i = 1
    while (N * x[i]/i < FDR & i <= N) i = i + 1
    if (i == 1)
        return(NA)
    else return(x[i - 1])
}

#' A Function That Builds a Fake Study.
#'
#' This function builds a fake study.
#'
#' @param nb_samples An integer that describes the number of samples.
#' @param nb_probes An integer that describes the number of probes.
#' @return a fake study.
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom survival Surv
#' @export
get_fake_study = function(nb_samples = 12, nb_probes = 10) {
  # generate data
  # data = matrix(c(rnorm(nb_probes * floor(nb_samples/2)), rnorm(nb_probes * ceiling(nb_samples/2), 3,1)), nb_probes)
  data = matrix(rnorm(nb_probes * nb_samples), nb_probes)
  data = apply(data, 2, function(c,c_shift) {c + c_shift}, runif(nb_probes, 1,10))
  data = signif(data,3)
  colnames(data) = c(paste(rep("ctrl", ceiling(nb_samples/2)), 1:ceiling(nb_samples/2), sep=""), paste(rep("case", floor(nb_samples/2)), 1:floor(nb_samples/2), sep=""))
  rownames(data) = paste(rep("prb", nb_probes), 1:nb_probes, sep="")

  # generate exp_grp
  exp_grp = data.frame(
    sex=ifelse(runif(nb_samples)>0.5, "Male", "Female"),
    age=round(runif(nb_samples,20,30)),
    tabac=c(rep("Smoker", ceiling(ceiling(nb_samples/2)/2)), rep("Non_Smoker", floor(ceiling(nb_samples/2)/2)), rep("Smoker", ceiling(floor(nb_samples/2)/2)), rep("Non_Smoker", floor(floor(nb_samples/2)/2))),
    treatment = c(rep("0ug", ceiling(nb_samples/2)), rep("15ug", floor(nb_samples/2))),
    histo=rep("lung", nb_samples),
    fut = round(runif(nb_samples,20,120)),
    da = round(runif(nb_samples))
  )
  sampled_rownames = sample(rownames(data))
  rownames(exp_grp) = colnames(data)

  # generate platform
  platform = data.frame(
    gene_name = sapply(1:nb_probes, function(i) {
      paste(toupper(c(sample(letters[1:26],round(runif(1,2,4))),sample(0:9,1), sample(letters[1:26],round(runif(1,0,1))))), collapse="")
    }),
    GOID = rep("...",nb_probes),
    stringsAsFactors=FALSE
  )
  rownames(platform) = rownames(data)

  # up and dwn with smokers
  idx_probes_up = sampled_rownames[1:ceiling(nb_probes/4)]
  idx_probes_dwn = rev(sampled_rownames)[1:floor(nb_probes/4)]
  idx_smoker_ctrl = rownames(exp_grp)[exp_grp$tabac == "Smoker" & exp_grp$treatment == "0ug"]
  idx_smoker_case = rownames(exp_grp)[exp_grp$tabac == "Smoker" & exp_grp$treatment == "15ug"]
  idx_nsmoker_case = rownames(exp_grp)[exp_grp$tabac == "Non_Smoker" & exp_grp$treatment == "15ug"]
  smoking_shift = runif(length(c(idx_probes_up, idx_probes_dwn)), 1.2,3)
  smoking_shift = smoking_shift * c(rep(1, length(idx_probes_up)), rep(-1, length(idx_probes_dwn)))

  # update dta
  data[c(idx_probes_up, idx_probes_dwn), idx_smoker_ctrl] = data[c(idx_probes_up, idx_probes_dwn), idx_smoker_ctrl] + smoking_shift
  data[c(idx_probes_up, idx_probes_dwn), idx_smoker_case] = data[c(idx_probes_up, idx_probes_dwn), idx_smoker_case] + smoking_shift/2
  data[c(idx_probes_up, idx_probes_dwn), idx_nsmoker_case] = data[c(idx_probes_up, idx_probes_dwn), idx_nsmoker_case] - smoking_shift/8
  data = data - min(data) + 1

  # update exp_grp
  exp_grp[idx_smoker_ctrl,]$da = ifelse(runif(length(idx_smoker_ctrl))>0.75, 0, 1)
  exp_grp[idx_smoker_case,]$da = ifelse(runif(length(idx_smoker_case))>0.66, 0, 1)
  exp_grp[idx_nsmoker_case,]$da = ifelse(runif(length(idx_nsmoker_case))>0.55, 0, 1)
  exp_grp$ss = Surv(exp_grp$fut, exp_grp$da)

  # update platform
  platform[idx_probes_up,]$gene_name = paste(platform[idx_probes_up,]$gene_name, "u", sep="_")
  platform[idx_probes_dwn,]$gene_name = paste(platform[idx_probes_dwn,]$gene_name, "d", sep="_")

  # build study
  study = create_study()
  study$data = data
  study$exp_grp = exp_grp
  study$platform = platform
  return(study)
}
#' A Function That Computes `mean` + 2 * `sd` on a Numeric Vector.
#'
#' This function computes `mean` + 2 * `sd` on a numeric vector.
#'
#' @param ctrl A numeric vector
#' @return `mean` + 2 * `sd` of the inpuit vector
#' @importFrom stats sd
#' @export
m2sd = function(ctrl) {
  mean(ctrl) + 2 * sd(ctrl)
}
#' A Function That Simplifies Sample Names.
#'
#' This function simplyfies sample names.
#'
#' @param sample_names A character vector that describes the sample names.
#' @return A character vector that describes simplified sample names.
#' @export
simplify_sample_names = function(sample_names) {
  sample_names = sub(".CEL.gz", "", sample_names, ignore.case = TRUE)
  sample_names = sub(".CEL", "", sample_names, ignore.case = TRUE)
  tmp_sample_names = as.vector(sapply(sample_names, function(gsm) {
    as.list(strsplit(gsm, "_")[[1]][1])
  }))
  if (sum(duplicated(tmp_sample_names)) == 0) {
    sample_names = tmp_sample_names
  }
  # sample_names = do.call(cbind, t(sample_names))
  return(sample_names)
}
#' A Function That Simplifies Factor Names.
#'
#' This function simplyfies factor names.
#'
#' @param factor_names A character vector that describes the factor names.
#' @param split character passed to strsplit function as split.
#' @param fixed boolean passed to gsub function as fixed.
#' @return A character vector that describes simplified factor names.
#' @export
simplify_factor_names = function(factor_names, split="", fixed=FALSE) {
  factor_names = as.character(factor_names)
  one_way = function(factor_names, split, fixed) {
    # print(factor_names)
    l = min(sapply(strsplit(unique(factor_names), split), length))
    f = lapply(strsplit(factor_names, split), function(x) {x[1:l]})
    m = do.call(rbind, f)
    vl = apply(m, 2, function(col) {
      length(unique(col))
    })
    d_idx = which(vl!=1)
    if (length(d_idx) == 0) {
      d_idx = length(vl)
    }
    fd_idx = d_idx[1] - 1
    if (fd_idx > 0) {
      to_be_remove = paste(c(f[[1]][1:fd_idx], ""), collapse=split)
      if (!fixed) {
        to_be_remove = paste("^", to_be_remove, sep="")
      }
      factor_names = gsub(to_be_remove, "", factor_names, fixed=fixed)
      # factor_names = gsub(to_be_remove, "", factor_names, fixed=fixed)
    }
    # factor_names = substr(factor_names, d_idx, 100000L)
    # print(factor_names)
    return(factor_names)
  }
  str_rev = function(factor_names) {
    sapply(lapply(strsplit(factor_names, ""), rev), function(s) { paste(s, collapse='')})
  }
  factor_names = one_way(factor_names, split, fixed)
  factor_names = str_rev(factor_names)
  factor_names = one_way(factor_names, split, fixed)
  factor_names = str_rev(factor_names)
  # print(factor_names)

  return(factor_names)
}
#' A Function That Simplifies Column Names.
#'
#' This function simplyfies column names.
#'
#' @param column_names A character vector that describes the column names.
#' @return A character vector that describes simplified column names.
#' @export
simplify_column_names = function(column_names) {
  n = tolower(column_names)
  n = gsub("[.|:| ]", "_", n)
  n = gsub("_+", "_", n)
  n = gsub("^_", "", n)
  n = gsub("_$", "", n)
  return(n)
}
#' Retrieve Directory's .CEL.gz Files.
#'
#' This function takes a directory as parameter and return a string vector of the .CEL.gz files
#' that it contains.
#'
#' @param dir A character string to consider.
#' @param full.names A logical value.
#' @param ERROR_IF_EMPTY A logical value that allow to not throw an error if dir do not contain any .CEL.gz file.
#' @return It returns a string vector of the .CEL.gz files.
#' @export
get_cel_filenames = function(dir, full.names=TRUE, ERROR_IF_EMPTY=TRUE) {
  raw_dir = paste(dir, "/raw", sep="")
  cel_filenames = list.files(raw_dir, pattern="*.CEL.gz", full.names=full.names, ignore.case = TRUE)
  if (length(cel_filenames) == 0) {
    cel_filenames = list.files(dir, pattern="*.CEL.gz", full.names=full.names, ignore.case = TRUE)
  }
  if (ERROR_IF_EMPTY & length(cel_filenames) == 0 ) { stop(paste("ERROR! No .CEL.gz file in ", paste(dir, collapse=" "), ".")) }
  return(cel_filenames)
}

#' Fuse two experimental grouping
#'
#' This function tales two experimental grouping as parameter and merge them according to their row names.
#' @param exp_grp1 The first experimental grouping.
#' @param exp_grp2 The second experimental grouping.
#' @param by Specifications of the columns used for merging. See 'Details' in ?merge.
#' @return The fused experimental grouping
#' @export
fuse_exp_grp = function(exp_grp1, exp_grp2, by="row.names") {
  simplify_exp_grp_rownames = function(exp_grp2, exp_grp1) {
    rownames(exp_grp1) = sub(".CEL.gz", "", rownames(exp_grp1), ignore.case = TRUE)
    rownames(exp_grp2) = sub(".CEL.gz", "", rownames(exp_grp2), ignore.case = TRUE)
    new_names = sapply(rownames(exp_grp2), function(n) {
      rownames_grep = unlist(sapply(rownames(exp_grp1), grep, n))
      if (length(rownames_grep) == 1) {
        return(names(rownames_grep))
      } else if (length(rownames_grep) > 1) {
        stop(paste("More than one rowname match ", n, ".", sep=""))
      } else {
        return(n)
      }
    })
    return(new_names)
  }
  rownames(exp_grp1) = simplify_exp_grp_rownames(exp_grp1, exp_grp2)
  rownames(exp_grp2) = simplify_exp_grp_rownames(exp_grp2, exp_grp1)

  fused_exp_grp = merge(exp_grp1, exp_grp2, by=by, all=TRUE)
  rownames(fused_exp_grp) = fused_exp_grp[,"Row.names"]
  fused_exp_grp[,"Row.names"] = NULL
  return(fused_exp_grp)
}



#' Retrieving _RAW.tar archive from GEO.
#'
#' This function retrieves _RAW.tar archgive from GEO to the file system.

#' @param gse A GSE id
#' @param datashare_dir A string describing the targeted directory.
#' @export
download_gse_raw_tar = function (gse, datashare_dir="~/projects/datashare") {
  dest_dir = paste0(datashare_dir, "/", gse)
  raw_dest_dir = paste0(dest_dir, "/raw")
  if (!file.exists(raw_dest_dir)) {
    command = "mkdir"
    args = paste("-p", raw_dest_dir)
    print(paste(command, args))
    system2(command, args)

    # gse = "GSE6791"
    gse_tar_filename = paste0(dest_dir, "/", gse, "_RAW.tar")
    gse_pref = paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
    url = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", gse_pref, "/", gse, "/suppl/", gse, "_RAW.tar")
    command = "wget"
    args = paste(url, "-P", dest_dir)
    print(paste(command, args))
    if (!file.exists(gse_tar_filename)) {
      system2(command, args)
    }

    tar_filename = paste0(dest_dir, "/", gse, "_RAW.tar")
    command = "tar"
    args = paste("xfv", tar_filename, "-C", raw_dest_dir)
    print(paste(command, args))
    system2(command, args)

    tar_filename = paste0(dest_dir, "/", gse, "_RAW.tar")
    command = "rm"
    args = tar_filename
    print(paste(command, args))
    system2(command, args)
  } else {
    print(paste(gse, " download ever done."))
  }
  return(dest_dir)
}



#' Retrieving .CEL.gz Files from GEO.
#'
#' This function retrieves .CEL.gz files from GEO to the file system. It takes as parameters a character vector of GMSs to be retrived and a targeted directory. It returns a character vector of .CEL.gz files.
#' @param gsms A character vector of GMSs to be retrived.
#' @param basedir A string describing the targeted directory.
#' @return A character vector of .CEL.gz files.
#' @export
#' @importFrom GEOquery getGEOSuppFiles
get_gsm = function(gsms, basedir) {
  "Retrieve .CEL.gz from NCBI GEO web site using GSM identifiers."
  dir.create(path=basedir, showWarnings=FALSE, recursive=TRUE)
  cel_filenames = get_cel_filenames(basedir, ERROR_IF_EMPTY=FALSE)
  sapply(gsms, function(gsm)  {
    gsm_grep = grep(gsm, cel_filenames)
    if (length(gsm_grep) == 0) {
      getGEOSuppFiles(gsm, baseDir=basedir, makeDirectory=FALSE)
    } else if (length(gsm_grep) != 1) {
      stop(paste("More than one .CEL.gz file match ", gsm, " in ", basedir, sep=""))
    } else {
      print(paste("Using locally cached version: ", cel_filenames[gsm_grep], " for ", gsm, ".", sep=""))
    }
  })
  return(get_cel_filenames(basedir))
}

#' Broken-stick Distribution for 'p' Pieces.
#'
#' Compute the expected values of the broken-stick distribution for `p` pieces.
#
#  @author Pierre Legendre
#' @param p An interger that specify the number of pieces
#' @return broken-stick distribution for `p` pieces.
#' @examples
#' broken_stick_out_20 = broken_stick(20)
#' @export
broken_stick <- function(p)
{
result = matrix(0,p,2)
colnames(result) = c("j","E(j)")
for(j in 1:p) {
   E = 0
   for(x in j:p) E = E+(1/x)
   result[j,1] = j
   result[j,2] = E/p
   }
return(result)
}



#' A memoised version of read.csv2
#'
#' This function offer a memoised version of read.csv2.
#'
#' @param ... parameters passed to read.csv2 function.
#' @return csv file content.
#' @importFrom GEOquery getGEO
#' @importFrom memoise memoise
mread.csv2 = memoise(function(...) {
  read.csv2(...)
})

#' A memoised version of getGEO
#'
#' This function offer a memoised version of getGEO. It is used to get GEO data
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to getGEO function.
#' @return GEO content.
#' @examples
#' # foo = mgetGEO("GSE42707", getGPL=FALSE)
#' @importFrom GEOquery getGEO
#' @importFrom memoise memoise
mgetGEO = memoise(function(...) {
  getGEO(...)
})

#' An other memoised version of getGEO
#'
#' This function offer a memoised version of getGEO. It returns a Table
#' version of a GEO query response. It is used to get platform metadatas
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to getGEO function.
#' @return GEO content.
#' @examples
#' # foo = mgetGEO("GPL13534", getGPL=FALSE)
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @importFrom memoise memoise
mtgetGEO = memoise(function(...) {
  Table(getGEO(...))
})


#' A Monitored Version of `apply`
#'
#' This function offer a monitored version of `apply`.
#'
#' @param X A matrix usualy gives as `apply` first parameter.
#' @param MARGIN An integer describing the marginal to use, usualy gives as  `apply` second parameter.
#' @param FUN A function usualy gives as `sapply` third parameter.
#' @param mod An integer that define the frequency of the monitoring.
#' @param ... Parameters passed to `func` function.
#' @return A vector of the application of the `func` function to each element of the `vec` vector.
# ' @examples
# ' # foo = monitored_apply(matrix(rnorm(50), 10), 1, function(v) {print(length(v)); Sys.sleep(1); return(c(mean(v), sd(v)))}, mod = 1)
#' @export
monitored_apply = function(X, MARGIN=1, FUN, mod=100, ...) {
  mat = X
  marg = MARGIN
  func = FUN
  epimedtools_global_cnt <<- 0
  nb_it = dim(mat)[marg]
  d1 = as.numeric(Sys.time())
  foo = apply(mat, marg, function(vec, ...) {
    epimedtools_global_cnt <<- epimedtools_global_cnt+1
    id = epimedtools_global_cnt
    ret = func(vec, ...)
    # monitoring...
    if (id %% mod == 0) {
      d2 = Sys.time()
      d2_num = as.numeric(Sys.time())
      elapse = d2_num - d1
      remain = (nb_it - id) * elapse / id
      end_at = d2 + remain
      print(paste("~", round(remain), " seconds remaining, finishing at ~", end_at, " (", id, "/", nb_it, ")." , sep=""))
    }
    return(ret)
  }, ...)
}
#' A Monitored Version of `sapply`
#'
#' This function offer a monitored version of `sapply`.
#'
#' @param vec A vector usualy gives as `sapply` first parameter.
#' @param func A function usualy gives as `sapply` second parameter.
#' @param ... Parameters passed to `monitored_apply` function.
#' @return A vector of the application of the `func` function to each element of the `vec` vector.
# ' @examples
# ' # foo = monitored_sapply(1:10, function(i) {print(i); Sys.sleep(1); return(c(i*i, i+i))}, mod = 1)
#' @export
monitored_sapply = function(vec, func, ...) {
  return(monitored_apply(t(vec), 2, func, ...))
}

#' A memoised version of readRDS
#'
#' This function offer a memoised version of readRDS. It is used to get GEO data
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to readRDS function.
#' @return RDS content.
#' @examples
#' # foo = mreadRDS("GSE42707", getGPL=FALSE)
#' @importFrom memoise memoise
mreadRDS = memoise(function(...) {
  readRDS(...)
})