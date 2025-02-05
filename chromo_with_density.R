
####### REPO update procedure: human genome from ensembl
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # connection to ENSEMBL
# all_features <- getBM(
#   attributes = c("entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name", "start_position", "end_position"),
#   mart = ensembl
# )
# write.table(all_features, "hsapiens_gene_ensembl.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# test <- read.delim("hsapiens_gene_ensembl.tsv")



###############################
#####    chromo class    ######
###############################

setClass(
  Class = "Chromo",
  slots = list(
    data = "data.frame",
    columns = "ANY",
    classification = "ANY",
    genome = "ANY",
    composition = "ANY",
    density = "ANY",
    interactions = "ANY"
  )
)




######################
# include cytobands df for non-human species? (its not widely spread, so no!) if(chrObj@genome$inp_annot != "default_human"){dont include cytobands!!!}


################################
#####    chromoInitiate    #####
################################

# 4 OPTIONS: use default ensembl (human), download new ensembl (any), provide ensembl (any), provide gtf (any)

chromoInitiate <- function(
    inp_annot = "default_human", # c("default_human", "new_biomart", "provide_biomart", "annotated", "gtf")
    species, # only use if using "new_biomart"
    ensembl_df, # only use if using "provide_biomart"
    gtf, # future
    DEdf,
    pval_cutoff = 0.05,
    log2fc_cutoff = 1,
    gene_col,
    fc_col,
    p_col,
    transcript_type = c("Mt_tRNA","Mt_rRNA","protein_coding","lncRNA","snRNA","rRNA","pseudogene","misc_RNA","processed_pseudogene",
                     "transcribed_unprocessed_pseudogene","rRNA_pseudogene","unprocessed_pseudogene","TEC","miRNA","transcribed_processed_pseudogene",
                     "snoRNA","unitary_pseudogene","transcribed_unitary_pseudogene","sRNA","IG_V_gene","IG_C_pseudogene","TR_J_gene",
                     "TR_V_gene","TR_J_pseudogene","TR_D_gene","TR_C_gene","IG_V_pseudogene","scaRNA","ribozyme","artifact","IG_C_gene",
                     "IG_J_gene","IG_J_pseudogene","IG_D_gene","translated_processed_pseudogene","TR_V_pseudogene","IG_pseudogene","vault_RNA","scRNA")
){

    if(inp_annot == "default_human"){

      all_features <- read.delim(system.file("extdata", "hsapiens_gene_ensembl.tsv", package = "chromo")) # repo updated frequently
      chromoObject <- new("Chromo", genome = list(inp_annot = inp_annot))

    }else if(inp_annot == "new_biomart"){

      ensembl <- useMart("ensembl", dataset = species) # connection to ENSEMBL
      all_features <- getBM(
        attributes = c("entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name", "start_position", "end_position"),
        filters = "external_gene_name",
        values = DEdf[[gene_col]],
        mart = ensembl
      )
      chromoObject <- new("Chromo", genome = list(inp_annot = inp_annot, species = species))

    }else if(inp_annot == "provide_biomart"){

      all_features <- ensembl_df
      chromoObject <- new("Chromo", genome = list(inp_annot = inp_annot))

    }else if(inp_annot == "gtf"){

      #gtf processing and cleaning first
      # all_features <- ...
      chromoObject <- new("Chromo", genome = list(inp_annot = inp_annot))
    }

    # if ensembl ever changes the name of the columns just change here to adapt the code!!!
    all_features <- all_features %>%
      filter(
        complete.cases(.),
        gene_biotype %in% c(transcript_type),
        chromosome_name %in% c(as.character(seq(1, 22)), "X", "Y"), #e se não for humano!!!! falar com helder se faz só humano ou não!!! Colocar um if no pipe se sim!
        !duplicated(external_gene_name) | !duplicated(external_gene_name, fromLast = TRUE) # choosing one of duplicated gene names
      ) %>%
      mutate(
        chromosome_name = factor(chromosome_name, levels = c(as.character(seq(1, 22)), "X", "Y")),
        gene_length = end_position - start_position,
        avg_position = (end_position + start_position)/2
      )

    DEdf <- DEdf %>%
      left_join(all_features, by = setNames("external_gene_name", gene_col)) %>%
      filter(
        complete.cases(.),
      ) %>%
      mutate(
        DEG = case_when(
          !!sym(fc_col) > log2fc_cutoff & !!sym(p_col) < pval_cutoff ~ "UP",
          !!sym(fc_col) < -log2fc_cutoff & !!sym(p_col) < pval_cutoff ~ "DOWN",
          TRUE ~ "NO"
        ),
        altered = case_when(
          !!sym(fc_col) > 0 & !!sym(p_col) < pval_cutoff ~ "UP",
          !!sym(fc_col) < 0 & !!sym(p_col) < pval_cutoff ~ "DOWN",
          TRUE ~ "NO"
        )
      ) %>%
      filter(!duplicated(!!sym(gene_col))) %>%
      arrange(chromosome_name, start_position)

    chromoObject@columns <- list( # if ensembl ever changes the name of the columns just change here to adapt the code!!! All downstream function should work!!!
      gene_col = gene_col, fc_col = fc_col, p_col = p_col,
      chromosome = "chromosome_name", start_position = "start_position", end_position = "end_position",
      avg_position = "avg_position", gene_length = "gene_length", DEG = "DEG", altered = "altered"
    )
    chromoObject@classification <- list(pval_cutoff = pval_cutoff, log2fc_cutoff = log2fc_cutoff)
    chromoObject@data <- DEdf

  return(chromoObject)
}

###################################
#####    chromoReclassify    ######
###################################

chromoReclassify <- function(
    chromoObject,
    pval_cutoff = 0.05,
    log2fc_cutoff = 1
){
  chromoObject@data <- chromoObject@data %>%
    mutate(
      DEG = case_when(
        !!sym(chromoObject@columns$fc_col) > log2fc_cutoff & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "UP",
        !!sym(chromoObject@columns$fc_col) < -log2fc_cutoff & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "DOWN",
        TRUE ~ "NO"
      ),
      altered = case_when(
        !!sym(chromoObject@columns$fc_col) > 0 & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "UP",
        !!sym(chromoObject@columns$fc_col) < 0 & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "DOWN",
        TRUE ~ "NO"
      )
    )

  chromoObject@classification <- list(pval_cutoff = pval_cutoff, log2fc_cutoff = log2fc_cutoff)

  return(chromoObject)
}

###################################
###      chromo Composition     ###
###################################

chromoComposition <- function(
    chromoObject,
    alteration = "DEG", # DEG or altered
    separate_by = "chromosome_name", #separate group
    only_expr_features = F, #calculates percentage/enrichment based on only expressed features (for Seurat only! When separate_by = "celltype") #### colocar isso numa lista?
    score_method = "pct", # "pct" or "hyp" or "hyp_padj"
    padj_method = "BH" # existing params of p.adjust() function # only if using "hyp_padj"
){

  aux <- chromoObject@data %>%
    {if (only_expr_features) filter(., pct.1 + pct.2 != 0) else .} # needs improvement or removal!!!

  if(score_method == "pct"){
    compo_df <- aux %>%
      group_by(!!sym(separate_by), !!sym(alteration)) %>%
      summarise(total = n()) %>%
      group_by(!!sym(separate_by)) %>%
      mutate(
        compo = (total / sum(total)) * 100
      ) %>%
      ungroup() %>%
      filter(
        !!sym(alteration) != "NO"
      )

  }else if(score_method %in% c("hyp","hyp_padj")){
    compo_df <- aux %>%
      group_by(!!sym(separate_by), !!sym(alteration)) %>%
      summarise(total = n())

    totals <- list()
    totals$UP <- aux %>% filter(!!sym(alteration) == "UP") %>% nrow()
    totals$DOWN <- aux %>% filter(!!sym(alteration) == "DOWN") %>% nrow()
    totals$NO <- aux %>% filter(!!sym(alteration) == "NO") %>% nrow()

    # pvalue calculation
    compo_df$compo <- apply(compo_df, 1, function(x){
      phyper(
        as.numeric(x["total"]) - 1,
        totals[[x[alteration]]],
        nrow(aux) - totals[[x[alteration]]],
        sum(compo_df$total[compo_df[[separate_by]] == x[separate_by]]),
        lower.tail = FALSE
      )
    })

    compo_df <- compo_df %>%
      filter(!!sym(alteration) != "NO")

    if(score_method == "hyp_padj"){
      compo_df$compo <- p.adjust(compo_df$compo, method = padj_method)
    }
  }

  chromoObject@composition <- list(
    compo_df = compo_df,
    alteration = alteration,
    separate_by = separate_by,
    only_expr_features = only_expr_features,
    score_method = score_method
  )

  if (score_method == "hyp_padj") {
    chromoObject@composition$padj_method <- padj_method
  }

  return(chromoObject)
}





#######################################
###      chromo Composition Plot    ###
#######################################

chromoCompositionPlot <- function(
    chromoObject,
    highlight_features = "top_by_separator", # "top_by_separator", "top_overall" or a list of names
    show_if_not_deg = T, # only if providing a list of names
    n_top_features = 1, # only if using "top_by_separator" or "top_overall"
    fc_line = T,
    title_xaxis = "Chromosome",
    title_yaxis = "Log2 fold change",

    color_dot_down = "#3771c8aa",
    color_dot_up = "#ff2200aa",
    color_dot_no = "#dddddd33",
    color_bar_down = "#3771c866",
    color_bar_up = "#ff220066",
    color_gene_name_down = "#2D5EAA",
    color_gene_name_up = "#aa0000",
    color_gene_name_no = "#555555",
    color_score_down = "#2D5EAA",
    color_score_up = "#aa0000",
    color_line_up = "#aa0000",
    color_line_down = "#2D5EAA",
    color_xaxis_text = "black",
    color_xaxis_label = "black",
    color_yaxis_text = "black",
    color_yaxis_label = "black",

    size_dot_alt = 1.2,
    size_dot_no = 0.8,
    size_bar = 0.8,
    size_gene_name = 3.2,
    size_score = ifelse(chromoObject@composition$score_method %in% c("hyp", "hyp_padj"), 7, 4),
    size_line = 0.4,
    size_xaxis_text = 14,
    size_xaxis_label = 12,
    size_yaxis_text = 12,
    size_yaxis_label = 12,

    style_gene_name = "plain",
    style_score = "bold",
    style_line = 2, # 1 - continuous, 2- dashed, 3 - dotted, etc.
    style_xaxis_text = "bold",
    style_xaxis_label = "plain",
    style_yaxis_text = "bold",
    style_yaxis_label = "plain",
    style_spacing = "" # Not working
){

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  alteration <- chromoObject@composition$alteration
  separate_by <- chromoObject@composition$separate_by
  compo_df <- chromoObject@composition$compo_df

  if(!is.null(chromoObject@composition)){

    aux <- chromoObject@data %>%
      {if (chromoObject@composition$only_expr_features) filter(., pct.1 + pct.2 != 0) else .} # needs improvement or removal!!!

    compo_df <- compo_df %>%
      mutate(
        compo = case_when(
          chromoObject@composition$score_method %in% c("hyp","hyp_padj") ~ -log10(compo),
          TRUE ~ compo
        )
      )

    max_compo <- max(compo_df$compo, na.rm = T) # debugged

    compo_df <- compo_df %>%
      mutate(
        proportion = (compo/max_compo * (max(abs(aux[[fc_col]])))),
        proportion = case_when(
          !!sym(alteration) == "UP" ~ proportion,
          !!sym(alteration) == "DOWN" ~ -proportion,
          TRUE ~ proportion
        ),
        compo = case_when(
          chromoObject@composition$score_method %in% c("hyp","hyp_padj") ~ ifelse(compo > -log10(0.001),"***",ifelse(compo > -log10(0.01),"**",ifelse(compo > -log10(0.05),"*",""))),
          TRUE ~ paste0(round(compo, 1), "%")
        ),
        y_axis = case_when(
          !!sym(alteration) == "UP" ~ 1.2 * max(aux[[fc_col]]),
          !!sym(alteration) == "DOWN" ~ 1.2 * min(aux[[fc_col]]),
          TRUE ~ NA
        ),
        color = case_when(
          !!sym(alteration) == "UP" ~ color_score_up,
          !!sym(alteration) == "DOWN" ~ color_score_down,
          TRUE ~ NA
        )
      )

    # highlight features
    if(highlight_features %in% c("top_by_separator", "top_overall")){
      max_genes <- aux %>%
        filter(!!sym(alteration) == "UP") %>%
        { if (highlight_features == "top_by_separator") group_by(., !!sym(separate_by)) else . } %>%
        slice_max(order_by = !!sym(fc_col), n = n_top_features, with_ties = FALSE) %>%
        { if (highlight_features == "top_by_separator") ungroup(.) else . } %>%
        dplyr::select(!!sym(gene_col), !!sym(separate_by))

      min_genes <- aux %>%
        filter(!!sym(alteration) == "DOWN") %>%
        { if (highlight_features == "top_by_separator") group_by(., !!sym(separate_by)) else . } %>%
        slice_min(order_by = !!sym(fc_col), n = n_top_features, with_ties = FALSE) %>%
        { if (highlight_features == "top_by_separator") ungroup(.) else . } %>%
        dplyr::select(!!sym(gene_col), !!sym(separate_by))

      hlf <- rbind(max_genes, min_genes) %>%
        mutate(gene_sep = paste0(!!sym(gene_col), "_", !!sym(separate_by))) %>%
        pull(gene_sep)
    }else{
      hlf <- highlight_features
    }

    local_aux <- aux %>%
      mutate(
        highlight = case_when(
          highlight_features %in% c("top_by_separator", "top_overall") & paste0(!!sym(gene_col), "_", !!sym(separate_by)) %in% hlf ~ !!sym(gene_col),
          !highlight_features %in% c("top_by_separator", "top_overall") & show_if_not_deg & !!sym(gene_col) %in% hlf ~ !!sym(gene_col),
          !highlight_features %in% c("top_by_separator", "top_overall") & !show_if_not_deg & !!sym(gene_col) %in% hlf & !!sym(alteration) %in% c("UP", "DOWN") ~ !!sym(gene_col),
          TRUE ~ NA
        ),
        color = case_when(
          !is.na(highlight) & !!sym(alteration) == "UP" ~ color_gene_name_up,
          !is.na(highlight) & !!sym(alteration) == "DOWN" ~ color_gene_name_down,
          !is.na(highlight) & !!sym(alteration) == "NO" ~ color_gene_name_no,
          TRUE ~ NA
        )
      )

    set.seed(42) # gene name position

    deg_plot <- ggplot()+
      geom_bar(
        data = compo_df,
        aes(x = !!sym(separate_by), y = proportion, fill = !!sym(alteration)),
        stat = "identity",
        width = size_bar
      ) +
      scale_fill_manual(values = c("UP" = color_bar_up, "DOWN" = color_bar_down)) # Custom colors for the bars

    if(fc_line){
      deg_plot <- deg_plot +
        geom_hline(yintercept = chromoObject@classification$log2fc_cutoff, color = color_line_up, size = size_line, linetype = style_line) +
        geom_hline(yintercept = -chromoObject@classification$log2fc_cutoff, color = color_line_down, size = size_line, linetype = style_line)
    }

    deg_plot <- deg_plot +
      geom_jitter( # not DEGs
        data = local_aux %>% filter(!!sym(alteration) == "NO"),
        aes(x = !!sym(separate_by), y = !!sym(fc_col), color = !!sym(alteration)),
        position = position_jitterdodge(
          jitter.width = 0.4,
          jitter.height = 0,
          dodge.width = 0,
          seed = 42 # reproducible results (fixed random number generator)
        ),
        size = size_dot_no
      )+
      geom_jitter( # DEGs
        data = local_aux %>% filter(!!sym(alteration) %in% c("UP", "DOWN")),
        aes(x = !!sym(separate_by), y = !!sym(fc_col), color = !!sym(alteration)),
        position = position_jitterdodge(
          jitter.width = 0.4,
          jitter.height = 0,
          dodge.width = 0,
          seed = 42 # reproducible results (fixed random number generator)
        ),
        size = size_dot_alt
      )+
      labs(
        x = title_xaxis,
        y = title_yaxis
      )+
      scale_color_manual(values = c("DOWN" = color_dot_down, "NO" = color_dot_no, "UP" = color_dot_up)) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = color_xaxis_text, size = size_xaxis_text, face = style_xaxis_text),
        axis.title.x = element_text(color = color_xaxis_label, size = size_xaxis_label, face = style_xaxis_label),
        axis.text.y = element_text(color = color_yaxis_text, size = size_yaxis_text, face = style_yaxis_text),
        axis.title.y = element_text(color = color_yaxis_label, size = size_yaxis_label, face = style_yaxis_label),
        legend.position = "none" # Remove legends
      )+
      geom_text_repel(
        data = local_aux,
        aes(x = !!sym(separate_by), y = !!sym(fc_col), color = !!sym(alteration), label = highlight),
        max.overlaps = Inf,
        color = local_aux$color,
        size = size_gene_name,
        fontface = style_gene_name
      ) +
      annotate( #score
        geom = "text",
        x = compo_df[[separate_by]],
        y = compo_df$y_axis,
        label = compo_df$compo,
        color = compo_df$color,
        size = size_score,
        fontface = style_score
      )

    return(deg_plot)

  }else{
    stop("Run chromoComposition first!")
  }
}


################################
#####    chromoDensity    ######
################################

chromoDensity <- function(
    chromoObject,
    bandwidth = "nrd0",
    cluster_threshold = 20, # 20%
    separate_up_down = F,
    calculate_bands = F, ### fix!

    weight_by = "none", #c("none", "length", "foldchange")
    scale_density = T,
    score_from_padj = T
){
  cytobands <- read.delim(system.file("extdata", "hsapiens_cytogenicbands.tsv", package = "chromo")) # bases obtained from UCSC's table browser

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  chromosome <- chromoObject@columns$chromosome
  start_position <- chromoObject@columns$start_position
  end_position <- chromoObject@columns$end_position
  avg_position <- chromoObject@columns$avg_position
  gene_length <- chromoObject@columns$gene_length
  DEG <- chromoObject@columns$DEG

  calculate_density <- function(subset, chro) {
    dens <- density(
      subset[[avg_position]],
      n = nrow(df), #reference cytobands?
      kernel = "gaussian",
      from = 1,
      to = max(cytobands[cytobands$chr == chro,"baseEnd"]),
      bw = bandwidth,
      weights = if (weight_by == "none") NULL else if (weight_by == "length") subset[[gene_length]] else subset[[fc_col]]
    )
    return(data.frame(x = dens$x, y = dens$y))
  }

  # All features density
  ALL_density <- chromoObject@data %>%
    group_by(!!sym(chromosome)) %>%
    do({
      dens_df <- calculate_density(., chro = unique(.data[[chromosome]]))
      dens_df
    }) %>%
    mutate(
      y = if_else(y < max(y)*(cluster_threshold/100), 0, y) ####### remover isso aqui?
    ) %>%
    ungroup() %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome)))

  # adding 0s on borders
  for (i in levels(ALL_density[[chromosome]])) {
    ALL_density <- ALL_density %>%
      add_row(!!sym(chromosome) := i, x = 0, y = 0) %>%
      add_row(!!sym(chromosome) := i, x = max(cytobands[cytobands$chr == i,"baseEnd"]), y = 0)
  }

  # ordering
  ALL_density <- ALL_density %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = levels(chromoObject@data[[chromosome]]))) %>%
    arrange(!!sym(chromosome), x)

  # Clustering
  ALL_clusters <- ALL_density[ALL_density$y != 0 | (ALL_density$y == 0 & (c(TRUE, ALL_density$y[-length(ALL_density$y)] != 0) | c(ALL_density$y[-1] != 0, TRUE))), , drop = FALSE]
  ALL_clusters <- as.data.frame(ALL_clusters)
  if(ALL_clusters$y[1] == 0 & ALL_clusters$y[2] == 0){ALL_clusters <- ALL_clusters[-1,]}
  if(ALL_clusters$y[nrow(ALL_clusters)] == 0 & ALL_clusters$y[nrow(ALL_clusters)-1] == 0){ALL_clusters <- ALL_clusters[-nrow(ALL_clusters),]}

  # list
  ALL_density_list <- list()
  boundaries <- which(ALL_clusters$y == 0)
  for (i in seq(1, length(boundaries), by=2)) {
    ALL_density_list[[(i+1)/2]] <- ALL_clusters[boundaries[i]:boundaries[i + 1], ]
  }

  # making clusters
  ALL_clusters <- map2_dfr(ALL_density_list, seq_along(ALL_density_list), ~ {
    df <- .x
    cluster_num <- .y
    data.frame(
      chromosome = unique(df[[chromosome]]),
      cluster_num = cluster_num,
      end_position = max(df$x),
      start_position = min(df$x)
    )
  }) %>%
    mutate(size = end_position - start_position)

  if(!separate_up_down){

    # DEGs density
    DEG_density <- chromoObject@data %>%
      filter(!!sym(DEG) != "NO") %>%
      group_by(!!sym(chromosome)) %>%
      do({
        dens_df <- calculate_density(., chro = unique(.data[[chromosome]]))
        dens_df
      }) %>%
      mutate(
        y = if_else(y < max(y)*(cluster_threshold/100), 0, y) # by chromosome and NOT or overall!!!
      ) %>%
      ungroup() %>%
      mutate(!!sym(chromosome) := factor(!!sym(chromosome)))

    # adding 0s on borders
    for (i in levels(DEG_density[[chromosome]])) {
      DEG_density <- DEG_density %>%
        add_row(!!sym(chromosome) := i, x = 0, y = 0) %>%
        add_row(!!sym(chromosome) := i, x = max(cytobands[cytobands$chr == i,"baseEnd"]), y = 0)
    }

    # ordering
    DEG_density <- DEG_density %>%
      mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = levels(chromoObject@data[[chromosome]]))) %>%
      arrange(!!sym(chromosome), x)

    # Clustering
    DEG_clusters <- DEG_density[DEG_density$y != 0 | (DEG_density$y == 0 & (c(TRUE, DEG_density$y[-length(DEG_density$y)] != 0) | c(DEG_density$y[-1] != 0, TRUE))), , drop = FALSE]
    DEG_clusters <- as.data.frame(DEG_clusters)
    if(DEG_clusters$y[1] == 0 & DEG_clusters$y[2] == 0){DEG_clusters <- DEG_clusters[-1,]}
    if(DEG_clusters$y[nrow(DEG_clusters)] == 0 & DEG_clusters$y[nrow(DEG_clusters)-1] == 0){DEG_clusters <- DEG_clusters[-nrow(DEG_clusters),]}

    # list
    DEG_density_list <- list()
    boundaries <- which(DEG_clusters$y == 0)
    for (i in seq(1, length(boundaries), by=2)) {
      DEG_density_list[[(i+1)/2]] <- DEG_clusters[boundaries[i]:boundaries[i + 1], ]
    }

    # remaking clusters df
    DEG_clusters <- map2_dfr(DEG_density_list, seq_along(DEG_density_list), ~ {
      df <- .x
      cluster_num <- .y
      data.frame(
        chromosome = unique(df[[chromosome]]),
        cluster_num = cluster_num,
        end_position = max(df$x),
        start_position = min(df$x)
      )
    }) %>%
      mutate(size = end_position - start_position)

    DEG_clusters$all_features <- NA
    DEG_clusters$DEGs <- NA
    for (i in 1:nrow(DEG_clusters)) { # change to apply!
      DEG_clusters$all_features[i] <- chromoObject@data %>%
        filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")
      DEG_clusters$DEGs[i] <- chromoObject@data %>%
        filter(!!sym(DEG) != "NO") %>%
        filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")
    }

    DEG_clusters <- DEG_clusters %>%
      mutate(
        n_features = str_count(all_features, ";") + 1,
        n_DEG = str_count(DEGs, ";") + 1,
        n_not_DEG = n_features - n_DEG,
      )

    # Fischer exact test
    total_DEG <- chromoObject@data %>% filter(!!sym(DEG) != "NO") %>% nrow()
    total_no_DEG <- chromoObject@data %>% filter(!!sym(DEG) == "NO") %>% nrow()
    DEG_clusters$pval <- apply(DEG_clusters, 1, function(row){
      matrix_data <- matrix(c(row["n_DEG"], row["n_not_DEG"], total_DEG - as.numeric(row["n_DEG"]), total_no_DEG - as.numeric(row["n_not_DEG"])), nrow = 2)
      matrix_data <- apply(matrix_data, 2, as.numeric)
      fish_res <- fisher.test(matrix_data)
      return(fish_res[["p.value"]])
    })

    # padj and score
    DEG_clusters$padj <- p.adjust(DEG_clusters$pval, method = "BH")
    DEG_clusters$score <- -log10(DEG_clusters$padj)

    # enrichment by band
    BANDS_clusters <- cytobands
    BANDS_clusters$all_features <- NA
    BANDS_clusters$DEGs <- NA
    for (i in 1:nrow(BANDS_clusters)) { # change to apply!
      BANDS_clusters$all_features[i] <- chromoObject@data %>%
        filter(!!sym(chromosome) == BANDS_clusters$chr[i],
               !!sym(avg_position) < BANDS_clusters$baseEnd[i],
               !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")

      BANDS_clusters$n_features[i] <- chromoObject@data %>%
        filter(!!sym(chromosome) == BANDS_clusters$chr[i],
               !!sym(avg_position) < BANDS_clusters$baseEnd[i],
               !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
        nrow()

      BANDS_clusters$DEGs[i] <- chromoObject@data %>%
        filter(!!sym(DEG) != "NO") %>%
        filter(!!sym(chromosome) == BANDS_clusters$chr[i],
               !!sym(avg_position) < BANDS_clusters$baseEnd[i],
               !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")

      BANDS_clusters$n_DEG[i] <- chromoObject@data %>%
        filter(!!sym(DEG) != "NO") %>%
        filter(!!sym(chromosome) == BANDS_clusters$chr[i],
               !!sym(avg_position) < BANDS_clusters$baseEnd[i],
               !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
        nrow()
    }

    BANDS_clusters <- BANDS_clusters %>%
      mutate(
        n_not_DEG = n_features - n_DEG
      ) %>%
      as.data.frame()

    # Fischer exact test
    total_DEG <- chromoObject@data %>% filter(!!sym(DEG) != "NO") %>% nrow()
    total_no_DEG <- chromoObject@data %>% filter(!!sym(DEG) == "NO") %>% nrow()
    BANDS_clusters$pval <- apply(BANDS_clusters, 1, function(row){
      matrix_data <- matrix(c(row["n_DEG"], row["n_not_DEG"], total_DEG - as.numeric(row["n_DEG"]), total_no_DEG - as.numeric(row["n_not_DEG"])), nrow = 2)
      matrix_data <- apply(matrix_data, 2, as.numeric)
      fish_res <- fisher.test(matrix_data)
      return(fish_res[["p.value"]])
    })

    # padj and score
    BANDS_clusters$padj <- p.adjust(BANDS_clusters$pval, method = "BH")
    BANDS_clusters$score <- -log10(BANDS_clusters$padj)

    # scaling
    if(scale_density){
      for(i in 1:nlevels(ALL_density[[chromosome]])){
        # aux vars
        max_y_all_density <- max(ALL_density[ALL_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
        if(nrow(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i],]) == 0){
          max_y_deg_density <- 0
        }else{
          max_y_deg_density <- max(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
        }
        max_density <- max(max_y_deg_density, max_y_all_density)

        # all clusters
        clusters_names_all <- ALL_clusters[ALL_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]
        max_chr_all <- bind_rows(ALL_density_list[clusters_names_all])
        max_chr_all <- max(max_chr_all$y)
        for (j in clusters_names_all) {
          ALL_density_list[[j]]$y = (ALL_density_list[[j]]$y) * (max_density) / max_chr_all
        }

        # DEGs
        clusters_names_degs <- DEG_clusters[DEG_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]
        if(length(clusters_names_degs) > 0){
          max_chr_deg <- bind_rows(DEG_density_list[clusters_names_degs])
          max_chr_deg <- max(max_chr_deg$y)
          for (j in clusters_names_degs) {
            DEG_density_list[[j]]$y = (DEG_density_list[[j]]$y) * (max_density) / max_chr_deg
          }
        }
      }
    }

    chromoObject@density$notSeparated = list(
      DEG_clusters = DEG_clusters,
      BANDS_clusters = BANDS_clusters,
      ALL_clusters = ALL_clusters, ####### remover isso aqui?
      DEG_density = DEG_density,
      ALL_density = ALL_density,
      parameters = list(bandwidth = bandwidth, threshold = cluster_threshold, weight_by = weight_by, scaled = scale_density)
    )

  }else if(separate_up_down){

    DEG_clusters_list <- list()
    BANDS_clusters_list <- list()
    DEG_density_aux <- list()
    ALL_density_aux <- list()

    for(k in c("UP", "DOWN")){
      # DEGs density
      DEG_density <- chromoObject@data %>%
        filter(!!sym(DEG) == k) %>%
        group_by(!!sym(chromosome)) %>%
        do({
          dens_df <- calculate_density(., chro = unique(.data[[chromosome]]))
          dens_df
        }) %>%
        mutate(
          y = if_else(y < max(y)*(cluster_threshold/100), 0, y) # by chromosome and NOT or overall!!!
        ) %>%
        ungroup() %>%
        mutate(!!sym(chromosome) := factor(!!sym(chromosome))) # não tirar!

      # adding 0s on borders
      for (i in levels(DEG_density[[chromosome]])) {
        DEG_density <- DEG_density %>%
          add_row(!!sym(chromosome) := i, x = 0, y = 0) %>%
          add_row(!!sym(chromosome) := i, x = max(cytobands[cytobands$chr == i,"baseEnd"]), y = 0)
      }

      # ordering
      DEG_density <- DEG_density %>%
        mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = levels(chromoObject@data[[chromosome]]))) %>%
        arrange(!!sym(chromosome), x)

      # Clustering
      DEG_clusters <- DEG_density[DEG_density$y != 0 | (DEG_density$y == 0 & (c(TRUE, DEG_density$y[-length(DEG_density$y)] != 0) | c(DEG_density$y[-1] != 0, TRUE))), , drop = FALSE]
      DEG_clusters <- as.data.frame(DEG_clusters)
      if(DEG_clusters$y[1] == 0 & DEG_clusters$y[2] == 0){DEG_clusters <- DEG_clusters[-1,]}
      if(DEG_clusters$y[nrow(DEG_clusters)] == 0 & DEG_clusters$y[nrow(DEG_clusters)-1] == 0){DEG_clusters <- DEG_clusters[-nrow(DEG_clusters),]}
      boundaries <- which(DEG_clusters$y == 0)
      DEG_density_list <- list()
      for (i in seq(1, length(boundaries), by=2)) {
        DEG_density_list[[(i+1)/2]] <- DEG_clusters[boundaries[i]:boundaries[i + 1], ]
      }

      DEG_clusters <- map2_dfr(DEG_density_list, seq_along(DEG_density_list), ~ {
        df <- .x
        cluster_num <- .y
        data.frame(
          chromosome = unique(df[[chromosome]]),
          cluster_num = cluster_num,
          end_position = max(df$x),
          start_position = min(df$x)
        )
      }) %>%
        mutate(size = end_position - start_position)

      DEG_clusters$all_features <- NA
      DEG_clusters$DEGs <- NA
      for (i in 1:nrow(DEG_clusters)) { # change to apply!
        DEG_clusters$all_features[i] <- chromoObject@data %>%
          filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
          pull(!!sym(gene_col)) %>%
          paste(collapse = ";")
        DEG_clusters$DEGs[i] <- chromoObject@data %>%
          filter(!!sym(DEG) == k) %>%
          filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
          pull(!!sym(gene_col)) %>%
          paste(collapse = ";")
      }

      DEG_clusters <- DEG_clusters %>%
        mutate(
          n_features = str_count(all_features, ";") + 1,
          n_DEG = str_count(DEGs, ";") + 1,
          n_not_DEG = n_features - n_DEG,
        )

      # Fischer exact test
      total_DEG <- chromoObject@data %>% filter(!!sym(DEG) == k) %>% nrow()
      total_no_DEG <- chromoObject@data %>% filter(!!sym(DEG) != k) %>% nrow()
      DEG_clusters$pval <- apply(DEG_clusters, 1, function(row){
        matrix_data <- matrix(c(row["n_DEG"], row["n_not_DEG"], total_DEG - as.numeric(row["n_DEG"]), total_no_DEG - as.numeric(row["n_not_DEG"])), nrow = 2)
        matrix_data <- apply(matrix_data, 2, as.numeric)
        fish_res <- fisher.test(matrix_data)
        return(fish_res[["p.value"]])
      })

      # padj and score
      DEG_clusters$padj <- p.adjust(DEG_clusters$pval, method = "BH")
      DEG_clusters$score <- -log10(DEG_clusters$padj)

      # saving df
      DEG_clusters_list[[k]] <- DEG_clusters

      # enrichment by band
      BANDS_clusters <- cytobands
      BANDS_clusters$all_features <- NA
      BANDS_clusters$DEGs <- NA
      for (i in 1:nrow(BANDS_clusters)) { # change to apply!
        BANDS_clusters$all_features[i] <- chromoObject@data %>%
          filter(!!sym(chromosome) == BANDS_clusters$chr[i],
                 !!sym(avg_position) < BANDS_clusters$baseEnd[i],
                 !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
          pull(!!sym(gene_col)) %>%
          paste(collapse = ";")

        BANDS_clusters$n_features[i] <- chromoObject@data %>%
          filter(!!sym(chromosome) == BANDS_clusters$chr[i],
                 !!sym(avg_position) < BANDS_clusters$baseEnd[i],
                 !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
          nrow()

        BANDS_clusters$DEGs[i] <- chromoObject@data %>%
          filter(!!sym(DEG) == k) %>%
          filter(!!sym(chromosome) == BANDS_clusters$chr[i],
                 !!sym(avg_position) < BANDS_clusters$baseEnd[i],
                 !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
          pull(!!sym(gene_col)) %>%
          paste(collapse = ";")

        BANDS_clusters$n_DEG[i] <- chromoObject@data %>%
          filter(!!sym(DEG) == k) %>%
          filter(!!sym(chromosome) == BANDS_clusters$chr[i],
                 !!sym(avg_position) < BANDS_clusters$baseEnd[i],
                 !!sym(avg_position) > BANDS_clusters$baseStart[i]) %>%
          nrow()
      }

      BANDS_clusters <- BANDS_clusters %>%
        mutate(
          n_not_DEG = n_features - n_DEG
        ) %>%
        as.data.frame()

      # Fischer exact test
      total_DEG <- chromoObject@data %>% filter(!!sym(DEG) == k) %>% nrow()
      total_no_DEG <- chromoObject@data %>% filter(!!sym(DEG) != k) %>% nrow()
      BANDS_clusters$pval <- apply(BANDS_clusters, 1, function(row){
        matrix_data <- matrix(c(row["n_DEG"], row["n_not_DEG"], total_DEG - as.numeric(row["n_DEG"]), total_no_DEG - as.numeric(row["n_not_DEG"])), nrow = 2)
        matrix_data <- apply(matrix_data, 2, as.numeric)
        fish_res <- fisher.test(matrix_data)
        return(fish_res[["p.value"]])
      })

      # padj and score
      BANDS_clusters$padj <- p.adjust(BANDS_clusters$pval, method = "BH")
      BANDS_clusters$score <- -log10(BANDS_clusters$padj)

      # saving df
      BANDS_clusters_list[[k]] <- BANDS_clusters

      # scaling
      if(scale_density){ ####### scaling should change density and not density_list!!! FIX!
        for(i in 1:nlevels(ALL_density[[chromosome]])){
          # aux vars
          max_y_all_density <- max(ALL_density[ALL_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
          if(nrow(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i],]) == 0){
            max_y_deg_density <- 0
          }else{
            max_y_deg_density <- max(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
          }
          max_density <- max(max_y_deg_density, max_y_all_density)

          # all clusters
          clusters_names_all <- ALL_clusters[ALL_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]
          max_chr_all <- bind_rows(ALL_density_list[clusters_names_all])
          max_chr_all <- max(max_chr_all$y)
          for (j in clusters_names_all) {
            ALL_density_list[[j]]$y = (ALL_density_list[[j]]$y) * (max_density) / max_chr_all
          }

          # DEGs
          clusters_names_degs <- DEG_clusters[DEG_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]
          if(length(clusters_names_degs) > 0){
            max_chr_deg <- bind_rows(DEG_density_list[clusters_names_degs])
            max_chr_deg <- max(max_chr_deg$y)
            for (j in clusters_names_degs) {
              DEG_density_list[[j]]$y = (DEG_density_list[[j]]$y) * (max_density) / max_chr_deg
            }
          }
        }
      }

      DEG_density_aux[[k]] <- DEG_density
      ALL_density_aux[[k]] <- ALL_density

    }

    chromoObject@density$separated = list(
      DEG_clusters = list(UP = DEG_clusters_list[["UP"]], DOWN = DEG_clusters_list[["DOWN"]]),
      BANDS_clusters = list(UP = BANDS_clusters_list[["UP"]], DOWN = BANDS_clusters_list[["DOWN"]]),
      DEG_density = list(UP = DEG_density_aux[["UP"]], DOWN = DEG_density_aux[["DOWN"]]),
      ALL_clusters = ALL_clusters,
      ALL_density = list(UP = ALL_density_aux[["UP"]], DOWN = ALL_density_aux[["DOWN"]]),
      parameters = list(bandwidth = bandwidth, threshold = cluster_threshold, weight_by = weight_by, scaled = scale_density)
    )
  }

  return(chromoObject)
}

####################################
#####    chromoPlotDensity    ######
####################################

chromoPlotDensity <- function(
    chromoObject,
    separate_up_down = F, #c("notSeparated", "separated")
    n_top_clusters = 10,
    n_top_bands = 7,
    color_enrich = "#990099",
    color_enrich_up = "#dd2200",
    color_enrich_down = "#0022dd"
){

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  chromosome <- chromoObject@columns$chromosome
  start_position <- chromoObject@columns$start_position
  end_position <- chromoObject@columns$end_position
  avg_position <- chromoObject@columns$avg_position
  gene_length <- chromoObject@columns$gene_length
  DEG <- chromoObject@columns$DEG

  if(!separate_up_down){
    # importing
    ALL_density <- chromoObject@density$notSeparated$ALL_density
    DEG_density <- chromoObject@density$notSeparated$DEG_density

    # creating ALL_density_list
    ALL_clusters <- ALL_density[ALL_density$y != 0 | (ALL_density$y == 0 & (c(TRUE, ALL_density$y[-length(ALL_density$y)] != 0) | c(ALL_density$y[-1] != 0, TRUE))), , drop = FALSE]
    ALL_clusters <- as.data.frame(ALL_clusters)
    if(ALL_clusters$y[1] == 0 & ALL_clusters$y[2] == 0){ALL_clusters <- ALL_clusters[-1,]}
    if(ALL_clusters$y[nrow(ALL_clusters)] == 0 & ALL_clusters$y[nrow(ALL_clusters)-1] == 0){ALL_clusters <- ALL_clusters[-nrow(ALL_clusters),]}
    ALL_density_list <- list()
    boundaries <- which(ALL_clusters$y == 0)
    for (i in seq(1, length(boundaries), by=2)) {
      ALL_density_list[[(i+1)/2]] <- ALL_clusters[boundaries[i]:boundaries[i + 1], ]
    }

    # creating DEG_density_list
    DEG_clusters <- DEG_density[DEG_density$y != 0 | (DEG_density$y == 0 & (c(TRUE, DEG_density$y[-length(DEG_density$y)] != 0) | c(DEG_density$y[-1] != 0, TRUE))), , drop = FALSE]
    DEG_clusters <- as.data.frame(DEG_clusters)
    if(DEG_clusters$y[1] == 0 & DEG_clusters$y[2] == 0){DEG_clusters <- DEG_clusters[-1,]}
    if(DEG_clusters$y[nrow(DEG_clusters)] == 0 & DEG_clusters$y[nrow(DEG_clusters)-1] == 0){DEG_clusters <- DEG_clusters[-nrow(DEG_clusters),]}
    DEG_density_list <- list()
    boundaries <- which(DEG_clusters$y == 0)
    for (i in seq(1, length(boundaries), by=2)) {
      DEG_density_list[[(i+1)/2]] <- DEG_clusters[boundaries[i]:boundaries[i + 1], ]
    }

    # importing
    DEG_clusters <- chromoObject@density$notSeparated[["DEG_clusters"]]
    BANDS_clusters <- chromoObject@density$notSeparated[["BANDS_clusters"]]
    ALL_clusters <- chromoObject@density$notSeparated[["ALL_clusters"]]

    # top clusters to include name
    top_clusters <- DEG_clusters %>% arrange(-score) %>% head(n_top_clusters)
    top_bands <- BANDS_clusters %>% arrange(-score) %>% head(n_top_bands) %>% pull(band)

    # plot
    plots_list <- c()
    for(i in 1:nlevels(ALL_density[[chromosome]])){

      # aux vars
      max_y_all_density <- max(ALL_density[ALL_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
      if(nrow(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i],]) == 0){
        max_y_deg_density <- 0
      }else{
        max_y_deg_density <- max(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
      }
      max_density <- max(max_y_deg_density, max_y_all_density)

      # plotting
      plots_list[[i]] <- ggplot() +
        geom_point( # not DEGs
          data = chromoObject@data %>% filter(!!sym(chromosome) == levels(ALL_density[[chromosome]])[i], !!sym(DEG) == "NO"),
          aes(x = !!sym(avg_position)), y = max_density * -0.1, color = "#00000022"
        ) +
        geom_point( # DEGs
          data = chromoObject@data %>% filter(!!sym(chromosome) == levels(ALL_density[[chromosome]])[i], !!sym(DEG) != "NO"),
          aes(x = !!sym(avg_position)), y = max_density * -0.1, color = paste0(color_enrich, "aa")
        ) +
        labs(
          title = NULL,
          x = NULL,
          y = paste0("chr", levels(ALL_density[[chromosome]])[i])) +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0)) + # remove padding to x axis
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(BANDS_clusters[BANDS_clusters$band == "chr1_q44","baseEnd"]))) + # remove padding to y axis and IMPORTANT fix size for all chr
        theme(
          plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
          panel.border = element_blank(),  # Remove panel borders
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
          #axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
          axis.title.y = element_text(face = "bold", size = 20, angle = 0, vjust = 0.5),  # Make axis titles bold
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.line.x = element_blank(),  # Remove x-axis line
          axis.line.y = element_blank()  # Remove y-axis line
          #axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
          #axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
          #axis.ticks.x = element_blank(),  # Remove x-axis ticks
          #axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
          #legend.position = "none" # No legend
        )

      # all featured polygons
      for (j in ALL_clusters[ALL_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
        plots_list[[i]] <- plots_list[[i]] +
          geom_polygon( # DEGs
            data = ALL_density_list[[j]],
            aes(x = x, y = y),
            color = "#333333dd",
            linewidth = 0.5,
            fill = "#bbbbbb77"
          )
      }

      # DEGs polygons
      for (j in DEG_clusters[DEG_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
        plots_list[[i]] <- plots_list[[i]] +
          geom_polygon( # DEGs
            data = DEG_density_list[[j]],
            aes(x = x, y = y),
            color = paste0(color_enrich, "dd"),
            linewidth = 0.5,
            fill = paste0(substr(colorRampPalette(c("white", color_enrich))(256)[as.integer(rescale(DEG_clusters[DEG_clusters$cluster_num == j, "score"], to = c(0, 1), from = c(min(DEG_clusters$score), max(DEG_clusters$score))) * 255) + 1], 1, 7), "cc")
            #paste0(substr(inferno(256)[as.integer(rescale(DEG_clusters[DEG_clusters$cluster_num == j,"score"], to = c(0, 1), from = c(min(DEG_clusters$score), max(DEG_clusters$score))) * 255) + 1],1,7), "cc")
          )
      }

      # cytogenetic bands
      for (j in BANDS_clusters[BANDS_clusters$chr == levels(ALL_density[[chromosome]])[i], "band"]) {
        plots_list[[i]] <- plots_list[[i]] +
          annotate(
            "rect",
            xmin = BANDS_clusters[BANDS_clusters$band == j, "baseStart"],
            xmax = BANDS_clusters[BANDS_clusters$band == j, "baseEnd"],
            ymin = ifelse(sub("^[^_]*_", "", j) %in% c("p11.1", "p11", "q11.1", "q11"), max_density * -0.5, max_density * -0.6),
            ymax = ifelse(sub("^[^_]*_", "", j) %in% c("p11.1", "p11", "q11.1", "q11"), max_density * -0.3, max_density * -0.2),
            fill = paste0(substr(colorRampPalette(c("white", color_enrich))(256)[as.integer(rescale(BANDS_clusters[BANDS_clusters$band == j, "score"], to = c(0, 1), from = c(min(BANDS_clusters$score), max(BANDS_clusters$score))) * 255) + 1], 1, 7), "cc"),
            color = color_enrich
          )
      }

      # cytogenetic bands top labels
      for (j in BANDS_clusters %>% filter(chr == levels(ALL_density[[chromosome]])[i], band %in% top_bands) %>% pull(band)) {
        plots_list[[i]] <- plots_list[[i]] +
          annotate(
            "text",
            x = (BANDS_clusters[BANDS_clusters$band == j, "baseStart"] + BANDS_clusters[BANDS_clusters$band == j, "baseEnd"])/2,
            y = (max_density * -0.6 + max_density * -0.2)/2,
            label = sub("^[^_]*_", "", j),
            color = "black",
            size = 1.5,
            angle = 90,
            fontface = "bold"
          )
      }

      # cluster name
      for (j in top_clusters[top_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
        plots_list[[i]] <- plots_list[[i]] +
          annotate(
            "text",
            x = (max(DEG_density_list[[j]]$x) + min(DEG_density_list[[j]]$x))/2,
            y = max(DEG_density_list[[j]]$y)*0.75,
            label = j,
            color = "black",
            size = 8,
            #angle = 45,
            fontface = "bold"
          )
      }

    }

    combined_plot <- purrr::reduce(plots_list, `/`)

  }else if(separate_up_down){

    plots_list_both <- list()

    for(k in c("UP", "DOWN")){

      # importing
      ALL_density <- chromoObject@density$separated$ALL_density[[k]]
      DEG_density <- chromoObject@density$separated$DEG_density[[k]]

      # creating ALL_density_list
      ALL_clusters <- ALL_density[ALL_density$y != 0 | (ALL_density$y == 0 & (c(TRUE, ALL_density$y[-length(ALL_density$y)] != 0) | c(ALL_density$y[-1] != 0, TRUE))), , drop = FALSE]
      ALL_clusters <- as.data.frame(ALL_clusters)
      if(ALL_clusters$y[1] == 0 & ALL_clusters$y[2] == 0){ALL_clusters <- ALL_clusters[-1,]}
      if(ALL_clusters$y[nrow(ALL_clusters)] == 0 & ALL_clusters$y[nrow(ALL_clusters)-1] == 0){ALL_clusters <- ALL_clusters[-nrow(ALL_clusters),]}
      ALL_density_list <- list()
      boundaries <- which(ALL_clusters$y == 0)
      for (i in seq(1, length(boundaries), by=2)) {
        ALL_density_list[[(i+1)/2]] <- ALL_clusters[boundaries[i]:boundaries[i + 1], ]
      }

      # creating DEG_density_list
      DEG_clusters <- DEG_density[DEG_density$y != 0 | (DEG_density$y == 0 & (c(TRUE, DEG_density$y[-length(DEG_density$y)] != 0) | c(DEG_density$y[-1] != 0, TRUE))), , drop = FALSE]
      DEG_clusters <- as.data.frame(DEG_clusters)
      if(DEG_clusters$y[1] == 0 & DEG_clusters$y[2] == 0){DEG_clusters <- DEG_clusters[-1,]}
      if(DEG_clusters$y[nrow(DEG_clusters)] == 0 & DEG_clusters$y[nrow(DEG_clusters)-1] == 0){DEG_clusters <- DEG_clusters[-nrow(DEG_clusters),]}
      DEG_density_list <- list()
      boundaries <- which(DEG_clusters$y == 0)
      for (i in seq(1, length(boundaries), by=2)) {
        DEG_density_list[[(i+1)/2]] <- DEG_clusters[boundaries[i]:boundaries[i + 1], ]
      }

      # importing
      ALL_clusters <- chromoObject@density$separated$ALL_clusters
      DEG_clusters <- chromoObject@density$separated$DEG_clusters[[k]]
      BANDS_clusters <- chromoObject@density$separated$BANDS_clusters[[k]]

      # top clusters to include name
      top_clusters <- DEG_clusters %>% arrange(-score) %>% head(n_top_clusters)
      top_bands <- BANDS_clusters %>% arrange(-score) %>% head(n_top_bands) %>% pull(band)

      # plot
      plots_list <- c()
      for(i in 1:nlevels(ALL_density[[chromosome]])){

        # aux vars
        max_y_all_density <- max(ALL_density[ALL_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
        if(nrow(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i],]) == 0){
          max_y_deg_density <- 0
        }else{
          max_y_deg_density <- max(DEG_density[DEG_density[[chromosome]] == levels(ALL_density[[chromosome]])[i], "y"])
        }
        max_density <- max(max_y_deg_density, max_y_all_density)

        # plotting
        plots_list[[i]] <- ggplot() +
          geom_point( # not DEGs
            data = chromoObject@data %>% filter(!!sym(chromosome) == levels(ALL_density[[chromosome]])[i], !!sym(DEG) != k),
            aes(x = !!sym(avg_position)),
            y = max_density * -0.1,
            color = "#00000022"
          ) +
          geom_point( # DEGs
            data = chromoObject@data %>% filter(!!sym(chromosome) == levels(ALL_density[[chromosome]])[i], !!sym(DEG) == k),
            aes(x = !!sym(avg_position)),
            y = max_density * -0.1,
            color = ifelse(k == "UP", paste0(color_enrich_up, "aa"), paste0(color_enrich_down, "aa"))
          ) +
          labs(
            title = NULL,
            x = NULL,
            y = paste0("chr", levels(ALL_density[[chromosome]])[i])) +
          theme_minimal() +
          scale_y_continuous(expand = c(0, 0)) + # remove padding to x axis
          scale_x_continuous(expand = c(0, 0), limits = c(0, max(BANDS_clusters[BANDS_clusters$band == "chr1_q44","baseEnd"]))) + # remove padding to y axis and IMPORTANT fix size for all chr
          theme(
            plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
            panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
            panel.border = element_blank(),  # Remove panel borders
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
            #axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
            axis.title.y = element_text(face = "bold", size = 20, angle = 0, vjust = 0.5),  # Make axis titles bold
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(),  # Remove x-axis ticks
            axis.title.x = element_blank(),  # Remove x-axis title
            axis.line.x = element_blank(),  # Remove x-axis line
            axis.line.y = element_blank()  # Remove y-axis line
            #axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
            #axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
            #axis.ticks.x = element_blank(),  # Remove x-axis ticks
            #axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
            #legend.position = "none" # No legend
          )

        # all featured polygons
        for (j in ALL_clusters[ALL_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
          plots_list[[i]] <- plots_list[[i]] +
            geom_polygon( # DEGs
              data = ALL_density_list[[j]],
              aes(x = x, y = y),
              color = "#333333dd",
              linewidth = 0.5,
              fill = "#bbbbbb77"
            )
        }

        # DEGs polygons
        for (j in DEG_clusters[DEG_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
          plots_list[[i]] <- plots_list[[i]] +
            geom_polygon(
              data = DEG_density_list[[j]],
              aes(x = x, y = y),
              color = paste0(ifelse(k == "UP", color_enrich_up, color_enrich_down), "dd"),
              linewidth = 0.5,
              fill = paste0(substr(colorRampPalette(c("white", ifelse(k == "UP", color_enrich_up, color_enrich_down)))(256)[as.integer(rescale(DEG_clusters[DEG_clusters$cluster_num == j, "score"], to = c(0, 1), from = c(min(DEG_clusters$score), max(DEG_clusters$score))) * 255) + 1], 1, 7), "cc")
              #paste0(substr(inferno(256)[as.integer(rescale(DEG_clusters[DEG_clusters$cluster_num == j,"score"], to = c(0, 1), from = c(min(DEG_clusters$score), max(DEG_clusters$score))) * 255) + 1],1,7), "cc")
            )
        }

        # cytogenetic bands
        for (j in BANDS_clusters[BANDS_clusters$chr == levels(ALL_density[[chromosome]])[i], "band"]) {
          plots_list[[i]] <- plots_list[[i]] +
            annotate(
              "rect",
              xmin = BANDS_clusters[BANDS_clusters$band == j, "baseStart"],
              xmax = BANDS_clusters[BANDS_clusters$band == j, "baseEnd"],
              ymin = ifelse(sub("^[^_]*_", "", j) %in% c("p11.1", "p11", "q11.1", "q11"), max_density * -0.5, max_density * -0.6),
              ymax = ifelse(sub("^[^_]*_", "", j) %in% c("p11.1", "p11", "q11.1", "q11"), max_density * -0.3, max_density * -0.2),
              fill = ifelse(k == "UP",
                            paste0(substr(colorRampPalette(c("white", color_enrich_up))(256)[as.integer(rescale(BANDS_clusters[BANDS_clusters$band == j, "score"], to = c(0, 1), from = c(min(BANDS_clusters$score), max(BANDS_clusters$score))) * 255) + 1], 1, 7), "cc"),
                            paste0(substr(colorRampPalette(c("white", color_enrich_down))(256)[as.integer(rescale(BANDS_clusters[BANDS_clusters$band == j, "score"], to = c(0, 1), from = c(min(BANDS_clusters$score), max(BANDS_clusters$score))) * 255) + 1], 1, 7), "cc")
              ),
              color = ifelse(k == "UP", color_enrich_up, color_enrich_down)
            )
        }

        # cytogenetic bands top labels
        for (j in BANDS_clusters %>% filter(chr == levels(ALL_density[[chromosome]])[i], band %in% top_bands) %>% pull(band)) {
          plots_list[[i]] <- plots_list[[i]] +
            annotate(
              "text",
              x = (BANDS_clusters[BANDS_clusters$band == j, "baseStart"] + BANDS_clusters[BANDS_clusters$band == j, "baseEnd"])/2,
              y = (max_density * -0.6 + max_density * -0.2)/2,
              label = sub("^[^_]*_", "", j),
              color = "black",
              size = 1.5,
              angle = 90,
              fontface = "bold"
            )
        }

        # cluster name
        for (j in top_clusters[top_clusters$chromosome == levels(ALL_density[[chromosome]])[i],"cluster_num"]) {
          plots_list[[i]] <- plots_list[[i]] +
            annotate(
              "text",
              x = (max(DEG_density_list[[j]]$x) + min(DEG_density_list[[j]]$x))/2,
              y = max(DEG_density_list[[j]]$y)*0.75,
              label = j,
              color = "black",
              size = 8,
              #angle = 45,
              fontface = "bold"
            )
        }

      }

      plots_list_both[[k]] <- plots_list

    }

    combined_plot <- purrr::reduce(plots_list_both[["UP"]], `/`) | purrr::reduce(plots_list_both[["DOWN"]], `/`)
  }

  return(combined_plot)
}

###############################
#####     chromoZoom     ######
############################### ###### add whole chr plot and plot by base start and end defined by user

chromoZoom <- function(
    chromoObject,
    density_slot = "notSeparated", #c("notSeparated","separated_up","separated_down")
    cluster,
    whole_chr = NULL
){

  custom_labels <- function(x) {
    ifelse(x >= 1e9, paste0(round(x / 1e9, 2), " Gb"),
           ifelse(x >= 1e6, paste0(round(x / 1e6, 2), " Mb"),
                  ifelse(x >= 1e3, paste0(round(x / 1e3, 2), " kb"), as.character(x))))
  }

  # full chromosome
  if (!is.null(whole_chr)){ # usar geom point instead
    zoom_plot <- ggplot()

    return(zoom_plot)
  }

  zoom_plot <- ggplot() +
    labs(
      title = paste0("Cluster ", cluster, ", chr", df[df$cluster_num == cluster, "chromosome"], ", ", custom_labels(df[df$cluster_num == cluster, "start_position"]), " - ", custom_labels(df[df$cluster_num == cluster, "end_position"])),
      x = NULL,
      y = expression(log[2]*FC)
    ) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) + # remove padding to x axis
    scale_x_continuous(expand = c(0, 0), labels = custom_labels) + # remove padding to y axis and IMPORTANT fix size for all chr
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
      panel.border = element_blank(),  # Remove panel borders
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.title = element_text( hjust = 0.5, face = "bold", size = 20),
      axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
      #axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
      axis.title.y = element_text(face = "bold", size = 20),  # Make axis titles bold
      # axis.text.y = element_blank(),
      # axis.text.x = element_blank(),  # Remove x-axis text
      # axis.ticks.x = element_blank(),  # Remove x-axis ticks
      # axis.title.x = element_blank(),  # Remove x-axis title
      # axis.line.x = element_blank()  # Remove x-axis line
      #axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
      #axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
      #axis.ticks.x = element_blank(),  # Remove x-axis ticks
      #axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
      legend.position = "none" # No legend
    )

  # min and max log2fc
  fc_vector <- DEdf %>% filter(Symbol %in% unlist(strsplit(df[df$cluster_num == cluster,"all_features"], ";"))) %>% pull(log2FoldChange)
  max_fc <- max(fc_vector)
  min_fc <- min(fc_vector)

  not_DEGs <- setdiff(unlist(strsplit(df[df$cluster_num == cluster,"all_features"], ";")), unlist(strsplit(df[df$cluster_num == cluster,"DEGs"], ";")))
  # not DEGs
  for (i in not_DEGs){
    zoom_plot <- zoom_plot +
      annotate("rect",
               xmin = DEdf[DEdf$Symbol == i,"start_position"],
               xmax = DEdf[DEdf$Symbol == i,"end_position"],
               ymin = DEdf[DEdf$Symbol == i,"log2FoldChange"] - (max_fc - min_fc)/30,
               ymax = DEdf[DEdf$Symbol == i,"log2FoldChange"],
               fill = '#77777777'
      )
  }

  # DEGs
  for (i in unlist(strsplit(df[df$cluster_num == cluster,"DEGs"], ";"))){
    zoom_plot <- zoom_plot +
      annotate("rect",
               xmin = DEdf[DEdf$Symbol == i,"start_position"],
               xmax = DEdf[DEdf$Symbol == i,"end_position"],
               ymin = DEdf[DEdf$Symbol == i,"log2FoldChange"] - (max_fc - min_fc)/30,
               ymax = DEdf[DEdf$Symbol == i,"log2FoldChange"],
               fill = ifelse(DEdf[DEdf$Symbol == i,"DEG"] == "DOWN", "#0033ffaa", "#ff3300aa")
      )
  }

  # gene name
  zoom_plot <- zoom_plot +
    geom_text_repel(
      data = DEdf %>% filter(Symbol %in% unlist(strsplit(df[df$cluster_num == cluster,"DEGs"], ";"))),
      aes(x = avg_position, y = log2FoldChange, label = Symbol, color = DEG),
      hjust = 0.5,
      vjust = 0.5,
      size = 3,
      fontface = "bold",
      inherit.aes = F
    )+
    scale_color_manual(values = c("DOWN" = "#002277", "UP" = "#772200"))

  return(zoom_plot)
}

####################################
#####    chromoInteractions    #####
####################################






#############################
###   Bar plot function   ###
#############################

barfunc <- function(
    df,
    x_axis,
    y_axis,
    fill_groups = F,
    annot_loc = "center", # "center" or "repel" or "top" or F
    fill_col = F,
    title_name = NULL,
    annot_text = y_axis,
    size_annot = 4,
    size_xaxis_text = 15,
    size_xaxis_title = 15,
    size_yaxis_text = 15,
    size_yaxis_title = 15,
    rotate_x = 0,
    x_title = F,
    y_title = F,
    extra_annot = NULL
){
  # Base plot
  bar_plot <- ggplot(df, aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
    geom_bar(stat = "identity") +
    labs(
      title = title_name,
      x = x_axis,
      y = y_axis
    ) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) + # remove extra padding between bars and x axis
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
      panel.border = element_blank(),  # Remove panel borders
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
      axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
      axis.title.y = element_text(face = "bold", size = size_yaxis_title),  # Make axis titles bold
      axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
      axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
      #legend.position = "none" # No legend
    )
  # y axis ticks
  if(size_yaxis_text == 0){
    bar_plot <- bar_plot + theme(axis.ticks.y = element_blank())
  }

  # Conditionally add fill aesthetic and scale_fill_manual
  if (fill_groups != F) {
    bar_plot <- bar_plot +
      aes(fill = !!sym(fill_groups)) +
      labs(fill = fill_groups)
  }

  # coloring bars
  if(is.list(fill_col)){
    bar_plot <- bar_plot +
      scale_fill_manual(values = fill_col)
  } else if (fill_col != F) {
    bar_plot <- bar_plot +
      geom_bar(stat = "identity", fill = fill_col) +
      theme(legend.position = "none")
  }

  # Annotations location
  if(!is.null(extra_annot)){
    bar_plot <- bar_plot +
      annotate(
        geom = "text",
        x = extra_annot[[x_axis]],
        y = extra_annot$y,
        label = extra_annot$annotation_text,
        color = extra_annot$color,
        size = 3,
        angle = 90,
        fontface = "bold"
      )
  }else if (annot_loc == "center") {
    bar_plot <- bar_plot +
      geom_text(
        aes(label = !!sym(annot_text)),
        position = position_stack(vjust = 0.5),
        size = size_annot
      )
  } else if (annot_loc == "repel") {
    bar_plot <- bar_plot +
      geom_text_repel(
        aes(label = !!sym(annot_text)),
        direction = "y",
        size = size_annot
      )
  } else if(annot_loc == "top"){
    bar_plot <- bar_plot +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.1*max(ifelse(fill_groups != F, df %>% group_by(!!sym(fill_groups)) %>% pull(!!sym(y_axis)), df[[y_axis]])))) +
      geom_text(
        aes(label = !!sym(annot_text)),
        vjust = -0.2,
        size = size_annot
      )
  }
  # x label title
  if(x_title != F){
    bar_plot <- bar_plot +
      xlab(x_title)
  }
  # y label title
  if(y_title != F){
    bar_plot <- bar_plot +
      ylab(y_title)
  }

  return(bar_plot)
}



####################################
###     Markers plot function    ###
####################################



####################################
###   Violin/Box plot function   ###
####################################



###############################
###      Volcano plot       ###
###############################



############################################
###      Functional enrich bar Plot      ###
############################################

funcenrichplot <- function(df, p_col = "p.adjust", highlight, color_list, text_size = 4, x_text_size = 13){

  df <- df %>%
    arrange(-!!sym(p_col)) %>%
    mutate(
      score = -log10(!!sym(p_col)),
      Description = factor(Description, levels = unique(df %>% arrange(-!!sym(p_col)) %>% pull(Description)))
    )

  bar_plot <- ggplot(df) +
    geom_col( # fill = ONTOLOGY ######### colorir de acordo com BP, CC, MF
      aes(x = score, y = Description, fill = !!sym(highlight)),
      #fill = "#4422aa77",
      width = 0.6
    ) +
    labs(
      title = NULL,
      y = NULL,
      x = expression("-log"[10]*"(padjust)")
    ) +
    scale_x_continuous(
      limits = c(0, max(df$score)),
      breaks = seq(0, max(df$score), by = 1),
      expand = c(0, 0), # The horizontal axis does not extend to either side
      position = "top"  # Labels are located on the top
    ) +
    # The vertical axis only extends upwards
    #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
    theme(
      # Set background color to white
      panel.background = element_rect(fill = "white"),
      # Set the color and the width of the grid lines for the horizontal axis
      panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
      # Remove tick marks by setting their length to 0
      axis.ticks.length = unit(0, "mm"),
      # Remove the title for both axes
      #axis.title = element_blank(),
      # Only left line of the vertical axis is painted in black
      axis.line.y.left = element_line(color = "black"),
      # Remove labels from the vertical axis
      axis.text.y = element_blank(),
      # But customize labels for the horizontal axis
      axis.text.x = element_text(family = "Econ Sans Cnd", size = x_text_size),
      # padj axis title
      axis.title.x = element_text(size = 13)
    ) +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      color = "#dd2222"
    ) +
    scale_fill_manual(
      values = color_list
    ) +
    geom_text( # descriptions
      aes(0, y = Description, label = Description),
      hjust = 0,
      nudge_x = 0.05,
      colour = "black",
      family = "Econ Sans Cnd",
      size = text_size
    )

  return(bar_plot)
}

##############################
#####    metaHeatmap    ######
##############################

sampDistri_func <- function(myMeta, xaxis, yaxis){
  if (yaxis == "alone"){
    myMeta$alone <- NA
  }

  # Creates long format dataframe with number of samples per xaxis and yaxis to do a heatmap
  long_df <- myMeta %>%
    group_by(!!sym(yaxis), !!sym(xaxis)) %>%
    summarise(num_Samples = n())

  # Reorder rows from most samples to least samples per row
  ordered_regions <- long_df %>%
    group_by(!!sym(yaxis)) %>%
    summarise(total_samples = sum(num_Samples)) %>%
    arrange(total_samples) %>%
    pull(!!sym(yaxis))

  long_df[[yaxis]] <- factor(long_df[[yaxis]], levels = ordered_regions)

  # Plots heatmap
  plot <- ggplot(long_df, aes(!!sym(xaxis), !!sym(yaxis), fill= num_Samples)) +
    geom_tile() +
    scale_fill_gradient2(high = "red") +
    xlab(xaxis) +
    ylab(yaxis) +
    geom_text(aes(label = num_Samples), vjust = +0.5, size = 12) +
    theme(
      text = element_text(size = 50),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.grid = element_line(color = "gray", linetype = "dotted"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(plot)
}

######################################
#####    featureDistribution    ######
######################################

######### adapt using ggdensity!!!!!!

# featureDistribution <- function(
    #     data,
#     feature,
#     bandwidth = "nrd0",
#     separate_by = NULL,
#     color = "#227722"
#   ){
#   data <- data %>%
#     filter(rownames(.) == feature) %>%
#     t() %>%
#     as.data.frame()
#
#   # density calculation
#   if(!is.null(separate_by)){
#     for(i in unique(separate_by)){
#       dens <- density(
#         data %>% dplyr::select(which(colnames(.) == i)) %>% as.vector(),
#         #n = length(data),
#         kernel = "gaussian",
#         from = 0,
#         to = max(data[1,]),
#         bw = bandwidth,
#       )
#       testing <- data.frame(x = dens$x, y = dens$y)
#     }
#   }else{
#     dens <- density(
#       data[,1],
#       #n = length(data),
#       kernel = "gaussian",
#       from = 0,
#       to = max(data[,1]),
#       bw = bandwidth,
#     )
#     dens <- data.frame(x = dens$x, y = dens$y) %>%
#       add_row(x = min(dens$x) - 1, y = 0) %>%
#       add_row(x = max(dens$x) + 1, y = 0) %>%
#       arrange(x)
#
#     max_density <- max(dens$y)
#   }
#
#   distPlot <- ggplot() +
#   geom_point(
#     data = data,
#     aes(x = !!sym(feature)), y = max_density * -0.1, color = paste0(color, "aa")
#   ) +
#   geom_polygon(
#     data = dens,
#     aes(x = x, y = y),
#     color = paste0(color, "aa"),
#     linewidth = 0.5,
#     fill = paste0(color, "88")
#   ) +
#   labs(
#     title = NULL,
#     x = NULL,
#     y = NULL
#   ) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) + # remove padding to x axis
#   scale_x_continuous(expand = c(0, 0)) + # remove padding to y axis and IMPORTANT fix size for all chr
#   theme(
#     plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
#     panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
#     panel.border = element_blank(),  # Remove panel borders
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
#     #axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
#     axis.title.y = element_text(face = "bold", size = 20, angle = 0, vjust = 0.5),  # Make axis titles bold
#     axis.text.y = element_blank(),
#     axis.text.x = element_blank(),  # Remove x-axis text
#     axis.ticks.x = element_blank(),  # Remove x-axis ticks
#     axis.title.x = element_blank(),  # Remove x-axis title
#     axis.line.x = element_blank(),  # Remove x-axis line
#     axis.line.y = element_blank()  # Remove y-axis line
#     #axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
#     #axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
#     #axis.ticks.x = element_blank(),  # Remove x-axis ticks
#     #axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
#     #legend.position = "none" # No legend
#   )
#
#   return(distPlot)
# }

