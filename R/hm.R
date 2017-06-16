plot_rna_sig = function(d_rna, pf, col=1, ylim, main="", ADD=FALSE, NORM=FALSE, PLOT_RAW=FALSE, PLOT_CI=FALSE, ALL_CURVE=FALSE, kern_dim=1, ...) {
  if (NORM) {
    main = paste(main, "density", sep=" ")
    d_rna  = t(t(d_rna)/apply(d_rna, 2, sum))
  } else {
    main = paste(main, "raw", sep=" ")
  }
  raw_sig = apply(d_rna, 1, function(l) {
    tt = t.test(l)
    if (tt$estimate==0) {
      ret = c(0,0,0,l)
    } else {
      ret = c(tt$estimate, tt$conf.int, l)
    }
    return(ret)
  })
  if (missing(ylim)) {
    ylim=range(raw_sig[1:3,], na.rm=TRUE)
  }
  k =kernel("daniell", kern_dim)
  smooth_sig = raw_sig
  for(i in 1:nrow(raw_sig)) {
    smooth_sig_tmp = kernapply(raw_sig[i,], k)
    nb_na = (length(raw_sig[i,]) - length(smooth_sig_tmp)) / 2
    smooth_sig[i,] = c(rep(NA, floor(nb_na)), smooth_sig_tmp, rep(NA, ceiling(nb_na)))
  }
  if (!ADD) {
    plot(0,0, col=0, xlim=c(1,ncol(raw_sig)), ylim=ylim, xaxt="n", main=main, ylab="signal (au)", ...)
    foo = pf[rownames(d_rna),]$chr[!is.na(pf[rownames(d_rna),]$chr)]
    idx = 1:length(foo)
    ndup = !duplicated(foo)
    abline(v=idx[ndup], col="grey", lty=3)
    axis(1, at=idx[ndup], labels=foo[ndup], las=2)
  }
  if (PLOT_RAW) {
    lines(raw_sig[1,], type="l", col="grey")
  }
  lines(smooth_sig[1,], type="l", col=adjustcolor(col, alpha.f=0.7), lty=1, lwd=3)
  # lines(smooth_sig[2,], type="l", col=col, lty=2)
  # lines(smooth_sig[3,], type="l", col=col, lty=2)
  if (ALL_CURVE) {
    apply(smooth_sig[-(1:3),], 1, lines, col=adjustcolor(col, alpha.f=0.7))
  }
  if (PLOT_CI) {
    y_sig_poly = c(smooth_sig[2,], rev(smooth_sig[3,]))
    x_sig_poly = c(1:ncol(smooth_sig), ncol(smooth_sig):1)
    x_sig_poly = x_sig_poly[!is.na(y_sig_poly)]
    y_sig_poly = y_sig_poly[!is.na(y_sig_poly)]
    polygon(x_sig_poly, y_sig_poly, lty=2, col=adjustcolor(col, alpha.f=0.3), border=col)
  }
}

#' A Function That runs Heatmap App
#'
#' This function runs heatmap app.
#' @param study An epimedtools Study RC object
#' @param height string specifying the height oh the app
#' @param var A vector defining the distances used to cluster data
#' @param USE_CLUST A boolean specifying if dendogram needs to be plotted
#' @param kern_dim the size of the kernel used to smooth RNAseq signal
#' @param ALL_CURVE A boolean specifying if individual RNAseq signal needs to be plotted
#' @param nb_clust An interger specifying the number of cluster
#' @param pf_col_label A platform column name use to label gene axis
#' @param exp_grp_col An exp_grp column name use to label sample
#' @param exp_grp_col_label exprimental
#' @param second exprimental
#' @param ... argument passed to plot_rna_sig function
#' @importFrom grDevices rainbow
#' @importFrom graphics lines
#' @importFrom graphics polygon
#' @importFrom stats as.dendrogram
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats kernapply
#' @importFrom stats kernel
#' @importFrom stats t.test
#' @importFrom shiny shinyApp
#' @importFrom shiny renderPlot
#' @importFrom shiny fluidPage
#' @importFrom shiny uiOutput
#' @importFrom shiny fluidRow
#' @importFrom shiny column
#' @importFrom shiny selectInput
#' @importFrom shiny checkboxGroupInput
#' @importFrom shiny renderPlot
#' @importFrom shiny updateCheckboxGroupInput
#' @importFrom shiny renderDataTable
#' @importFrom shiny actionLink
#' @importFrom shiny checkboxGroupInput
#' @importFrom shiny checkboxInput
#' @importFrom shiny column
#' @importFrom shiny dataTableOutput
#' @importFrom shiny div
#' @importFrom shiny fluidPage
#' @importFrom shiny fluidRow
#' @importFrom shiny h3
#' @importFrom shiny h4
#' @importFrom shiny insertUI
#' @importFrom shiny numericInput
#' @importFrom shiny observe
#' @importFrom shiny plotOutput
#' @importFrom shiny radioButtons
#' @importFrom shiny reactiveValues
#' @importFrom shiny removeUI
#' @importFrom shiny renderDataTable
#' @importFrom shiny renderPlot
#' @importFrom shiny renderUI
#' @importFrom shiny selectInput
#' @importFrom shiny uiOutput
#' @importFrom shiny updateCheckboxGroupInput
#' @export
get_hm_app = function(study, height="5000", var=c("samples", "raw", "genome"), USE_CLUST=c(TRUE, FALSE), kern_dim=1, ALL_CURVE=FALSE, nb_clust=3, pf_col_label="exon", exp_grp_col="histo", exp_grp_col_label, second, ...) {
  if (missing(second)) {
    second = exp_grp_col
  }
  if (missing(exp_grp_col_label)) {
    exp_grp_col_label = exp_grp_col
  }
  

  hm_app = shiny::shinyApp(

    ui = fluidPage (title="FASTKD1 namalwan",
      # selectInput("study", label="Cohort:", choices=study$cache_filename, selected=study$cache_filename),

      uiOutput("index_controls"), 
      
      selectInput("exp_grp_col_label",
        label="Experimental grouping column to use for labeling:",
        choices=colnames(study$exp_grp), selected=exp_grp_col_label
      ),
      selectInput("pf_col_label",
        label="Platform column to use for labeling:",
        choices=colnames(study$platform), selected=pf_col_label
      ),

      plotOutput("gene_genome"),

      selectInput("method_hclust",
        label="clustering method",
        choices=c("single", "complete", "mcquitty", "ward.D", "ward.D2", "average", "median", "centroid"),
        selected="complete"),
      numericInput("nb_clust", "nb. clusters:", 1, min = 1, max = 50, value = nb_clust),
      checkboxInput("PLOT_HM_raw", "plot raw heatmap", TRUE),
      radioButtons("raw_or_cor", "distance to use:", var),
      radioButtons("USE_CLUST", "use clustering", USE_CLUST),
      plotOutput("heatmap_raw"),

      # plotOutput("survival_clust_raw"),
      #
      # checkboxInput("MEAN_SIG_raw", "mean signal",    FALSE),
      # checkboxInput("PLOT_CI_raw", "conf. int.",      FALSE),
      # checkboxInput("ALL_CURVE", "indiv. curves", ALL_CURVE),
      # checkboxGroupInput("selected_clusters", "selected clusters:", 1:10),
      # numericInput("kern_dim", "smoothing:", 1, min = 1, max = 50, value = kern_dim),
      # plotOutput("rna_seq_clust_raw"),

      dataTableOutput("clusters")

      # , downloadButton("report", "Generate report")

    ),


    server = function(input, output, session) {
      react_vals = reactiveValues()
      filter_keys = c(exp_grp="exp_grp_col_filter", platform="platform_col_filter")

      output$index_controls = renderUI({
        lapply(names(filter_keys), function(key){
          snipet_df_filter(key)$ui
        })
      })      

      snipet_df_filter = function(key) {
        idx = apply(study[[key]], 2, function(c){
          length(unique(c)) != 1 & length(unique(c)) != nrow(study[[key]])
        })
        idx[1] = TRUE
        ui = fluidRow(
          column(div(
            h3(paste0("Filtering ", key)), 
            checkboxGroupInput(paste0(key, "_col_filter"),
              label="Which filter?",
              choices=colnames(study[[key]])[idx]
            )),
          width=6),
          column(div(id=paste0(key, "_col_filter_items")), 
          width=6)
        )
        obs = observe({})
        return(list(ui=ui, obs=obs))        
      }

      update_index_col_filter_items = observe({
        lapply(names(filter_keys), function(key){        
          if (!is.null(react_vals[[paste0(key, "_col_filter")]])) {
            to_insert = input[[paste0(key, "_col_filter")]][!input[[paste0(key, "_col_filter")]] %in% react_vals[[paste0(key, "_col_filter")]]]
            to_remove = react_vals[[paste0(key, "_col_filter")]][!react_vals[[paste0(key, "_col_filter")]] %in% input[[paste0(key, "_col_filter")]]]
          } else {
            to_insert = input[[paste0(key, "_col_filter")]]
            to_remove = NULL
          }
          # print(paste0("to_insert:", to_insert))
          # print(paste0("to_remove:", to_remove))
          for (col_name in to_insert) {
            snipet = snipet_colum_filter(col_name, key)
            obs = snipet$obs
            insertUI(
              selector = paste0("#", key, "_col_filter_items"),
              where = "afterEnd",
              ui = snipet$ui
            )
          }
          for (col_name in to_remove) {
            removeUI(selector = paste0("#", key, "_col_filter_div_", col_name))
          }
          react_vals[[paste0(key, "_col_filter")]] = input[[paste0(key, "_col_filter")]]            
        })
      })

      update_indexes = observe({
        for (key in names(filter_keys)) { 
          if (is.null(react_vals[[paste0(key, "_idx")]])) {
            react_vals[[paste0(key, "_idx")]] = rep(TRUE, nrow(study[[key]]))
          }
          if (!is.null(react_vals[[paste0(key, "_col_filter")]])) {
            # print(react_vals[[paste0(key, "_col_filter")]])
            if (length(react_vals[[paste0(key, "_col_filter")]]) == 0) {
              ret = rep(TRUE, nrow(study[[key]]))
            } else {
              ret = sapply(react_vals[[paste0(key, "_col_filter")]], function(col_name){
                # print(paste(key, react_vals[[paste0(key, "_col_filter")]], input[[paste0(key, "_col_filter_", col_name)]]))
                idx = !is.na(study[[key]][[col_name]]) & study[[key]][[col_name]] %in% input[[paste0(key, "_col_filter_", col_name)]]
                return(idx)
              })
             ret = apply(ret,1,all)
            }
            if (!all(react_vals[[paste0(key, "_idx")]] == ret)) {
              react_vals[[paste0(key, "_idx")]] = ret                            
            }
          }
        }
      })

      snipet_colum_filter = function (col_name, key) {
        ui = div(id=paste0(key, "_col_filter_div_", col_name),
          h4(col_name), 
          checkboxGroupInput(
            inputId = paste0(key, "_col_filter_", col_name), , 
            label   = NULL,
            choices = sort(unique(study[[key]][[col_name]])), 
            selected = sort(unique(study[[key]][[col_name]]))
          ),
          actionLink(
            inputId = paste0("selectall_", key, "_col_filter_", col_name), 
            label   = "Select All"
          )
        ) 
        obs = observe({
            if (!is.null(input[[paste0("selectall_", key, "_col_filter_", col_name)]])) {
              if (input[[paste0("selectall_", key, "_col_filter_", col_name)]] == 0) {
              } else if (input[[paste0("selectall_", key, "_col_filter_", col_name)]]%%2 == 0) {
                updateCheckboxGroupInput(session, paste0(key, "_col_filter_", col_name),
                  choices = sort(unique(study[[key]][[col_name]]))
                )
              } else {
                updateCheckboxGroupInput(session, paste0(key, "_col_filter_", col_name),
                  choices = sort(unique(study[[key]][[col_name]])),
                  selected = sort(unique(study[[key]][[col_name]]))
                )
              }                
            }
        })                 
        return(list(ui=ui, obs=obs))
      }


        ##############
       # GRAPHICALS #
      ##############
      
      output$gene_genome = renderPlot({
        # default
        platform_idx = rownames(study$platform)
        exp_grp_idx = rownames(study$exp_grp)
        pf_col_label = rev(colnames(study$platform))[1]
        # MVC
        if (!is.null(react_vals$exp_grp_idx) &
            !is.null(react_vals$platform_idx) &
            !is.null(input$pf_col_label)) {
          if (sum(react_vals$exp_grp_idx )!=0 &
              sum(react_vals$platform_idx)!=0) {
            pf_col_label = input$pf_col_label
            platform_idx = react_vals$platform_idx
            exp_grp_idx = react_vals$exp_grp_idx

            # plot
            d = study$data[platform_idx, exp_grp_idx]
            xs = 1:nrow(study$platform[platform_idx, ])
            ys = apply(d, 1, mean)
            layout(1)
            plot(0,0,col=0, xlim=range(xs), ylab="signal (au)", xlab="", main=paste0("signal over genome ", nrow(d), "x", ncol(d)), ylim=range(ys), xaxt="n")
            labs = study$platform[rownames(d),pf_col_label]
            idx_lab = !duplicated(labs)
            axis(1, at=xs[idx_lab], labels=labs[idx_lab], las=2)
            abline(v=xs[idx_lab], lty=3, col="grey")
            lines(xs, ys, lwd=2, col=1, type="s")
          }
        }
      })

      output$heatmap_raw = renderPlot({
        # default
        exp_grp_col = rev(colnames(study$exp_grp))[1]
        second = exp_grp_col
        exp_grp_col_label = exp_grp_col
        pf_col_label = rev(colnames(study$platform))[1]
        exp_grp_idx = rownames(study$exp_grp)
        platform_idx = rownames(study$platform)
        method_hclust = "complete"
        nb_clust = 2
        var = "genome"
        PLOT_HM_raw = TRUE
        USE_CLUST = FALSE
        if (!is.null(react_vals$exp_grp_idx) &
            !is.null(react_vals$platform_idx)) {
          if (sum(react_vals$exp_grp_idx )!=0 &
              sum(react_vals$platform_idx)!=0) {
            # MVC
            exp_grp_col = input$exp_grp_col
            exp_grp_col_label = input$exp_grp_col_label
            second = input$second
            pf_col_label = input$pf_col_label
            exp_grp_idx = react_vals$exp_grp_idx
            platform_idx = react_vals$platform_idx
            method_hclust = input$method_hclust
            USE_CLUST = input$USE_CLUST
            nb_clust = input$nb_clust
            var = input$raw_or_cor
            PLOT_HM_raw = input$PLOT_HM_raw

            # Heatmap
            hm = plot_hm2(study, pf_col=pf_col_label, exp_grp_col_label=exp_grp_col_label, exp_grp_idx=exp_grp_idx, platform_idx=platform_idx, method_hclust=method_hclust, nb_clust=nb_clust, var=var, PLOT_HM_raw=PLOT_HM_raw, USE_CLUST=USE_CLUST)

            # MVC
            react_vals$grps = hm$grps
            updateCheckboxGroupInput(session, "selected_clusters",
              choices = unique(hm$grps),
              selected = unique(hm$grps)
            )
          }
        }
      })

      output$survival_clust_raw = renderPlot({
        if (!is.null(react_vals$grps) & 
            ("efs" %in% names(study$exp_grp) |
             "os" %in% names(study$exp_grp))
        ) {
          grps = react_vals$grps
          layout(matrix(1:2, 1), respect=TRUE)
          cols = rainbow(length(unique(grps)))[1:length(unique(grps))]
          if (!is.null(study$exp_grp[names(grps ),"efs"])) {
            scurve(ss=study$exp_grp[names(grps ),"efs"], v=as.factor(grps), main="efs", colors=cols)
          }
          if (!is.null(study$exp_grp[names(grps ),"os"])) {
            scurve(ss=study$exp_grp[names(grps ),"os"], v=as.factor(grps), main="os", colors=cols)
          }
        }
      })

      output$rna_seq_clust_raw = renderPlot({
        # default
        exp_grp_idx = rownames(study$exp_grp)
        platform_idx = rownames(study$platform)
        MEAN_SIG_raw = FALSE
        PLOT_CI_raw = TRUE
        USE_CLUST = FALSE
        ALL_CURVE = FALSE
        raw_or_cor = "raw"
        kern_dim = 1
        selected_clusters = 1:2

        # MVC
        if (!is.null(react_vals$exp_grp_idx) &
            !is.null(react_vals$platform_idx)) {
          if (length(react_vals$exp_grp_idx)!=0 &
              length(react_vals$platform_idx  )!=0) {
            selected_clusters = input$selected_clusters
            grps = react_vals$grps
            kern_dim = input$kern_dim
            exp_grp_idx = react_vals$exp_grp_idx
            platform_idx = react_vals$platform_idx
            MEAN_SIG_raw = input$MEAN_SIG_raw
            PLOT_CI_raw = input$PLOT_CI_raw
            ALL_CURVE = input$ALL_CURVE
            selected_clusters = input$selected_clusters
            raw_or_cor = input$raw_or_cor
          }
        }

        # plot
        if (!is.null(react_vals$grps)) {
          if (raw_or_cor != "genome") {
            col_rainbow = rainbow(10)
            layout(matrix(1:2, 1), respect=TRUE)
            for (NORM in c(FALSE, TRUE)) {
              col_mean_sig = ifelse(MEAN_SIG_raw, "black", "white")
              plot_rna_sig(study$data[platform_idx, exp_grp_idx], study$platform[platform_idx, ], NORM=NORM, col=col_mean_sig,  PLOT_CI=PLOT_CI_raw, ALL_CURVE=ALL_CURVE, kern_dim=kern_dim)
              is = as.numeric(selected_clusters)
              if (ALL_CURVE) {
                for (i in is) {
                  idx_tmp = names(grps)[grps==i]
                  plot_rna_sig(study$data[platform_idx, idx_tmp], , study$platform[platform_idx, ], col=rainbow(length(unique(grps)))[i], ADD=TRUE, NORM=NORM, PLOT_CI=FALSE, ALL_CURVE=TRUE, kern_dim=kern_dim)
                }
              }
              if (PLOT_CI_raw) {
                for (i in is) {
                  idx_tmp = names(grps)[grps==i]
                  plot_rna_sig(study$data[platform_idx, idx_tmp], , study$platform[platform_idx, ], col=rainbow(length(unique(grps)))[i], ADD=TRUE, NORM=NORM, PLOT_CI=TRUE, ALL_CURVE=FALSE, kern_dim=kern_dim)
                }
              }
              for (i in is) {
                idx_tmp = names(grps)[grps==i]
                plot_rna_sig(study$data[platform_idx, idx_tmp], , study$platform[platform_idx, ], col=rainbow(length(unique(grps)))[i], ADD=TRUE, NORM=NORM, PLOT_CI=FALSE, ALL_CURVE=FALSE, kern_dim=kern_dim)
              }
            }
          }          
        } else {
          return(NULL)
        }
      })

      output$clusters = renderDataTable({
        # default
        # MVC
        grps = react_vals$grps
        exp_grp_col_label = input$exp_grp_col_label
        # plot
        exp_grp = study$exp_grp
        exp_grp$grps = NA
        exp_grp[names(grps), "grps"] = grps
        exp_grp = exp_grp[order(exp_grp$grps),]
        table(exp_grp[,"grps"], exp_grp[,exp_grp_col_label])
      })

    },
    options = list(height=height),
  )

  return(hm_app)
}