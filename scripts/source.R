# Find inflection points for selection of principle components
get_npcs <- function(seurat_object, create_plot = T){
  library(ggplot2)
  
  std_data = Stdev(object = seurat_object, reduction = "pca")
  ndims = length(x = std_data)
  elbow_data = data.frame(dims = 1:ndims, stdev = std_data[1:ndims])
  
  reference = elbow_data[,2]
  difference_vect = c(elbow_data[2:nrow(elbow_data),2],0)
  difference = which((reference - difference_vect) < 0.05)
  
  difference_consec = c(difference[2:length(difference)],0) - difference
  names(difference_consec) = difference
  
  npcs = as.numeric(names(which(difference_consec ==1)[1]))
  
  if(create_plot){
    
    plt = ggplot(elbow_data, aes(x = dims, y = stdev)) +
      geom_point() + 
      geom_vline(xintercept=npcs) +
      geom_text(aes(npcs+10, max(stdev), label=paste("Inferred NPCs:", npcs), vjust=0), size=6) +
      xlab("Number of principle components") +
      ylab("Standard deviation of principle components") +
      theme(text = element_text(size=17), axis.text=element_text(size=12))
    
  }
  
  out <- list(npcs=npcs, plot=plt)
  return(out)
  
}

# Assign cell identities using file of marker genes and likely cell types
assign.identity <- function(seurat_object, markers){
  
  #Find significant marker genes for each cluster
  all.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  #Create one row for each cell.type label and marker gene
  identity <- strsplit(markers$Markers, ",")
  markers <- data.frame(Cell.Type = rep(markers$Cell.Type, sapply(identity, length)), Markers = unlist(identity))
  
  #Make sure both gene names are in sentence case
  markers$Markers <- str_to_sentence(markers$Markers)
  all.markers$gene <- str_to_sentence(all.markers$gene)
  
  #Merge identity on if available
  add <- merge(markers, all.markers, by.x="Markers", by.y="gene", all.y=TRUE)
  
  #Subset identity and cluster numbers
  ids <- add %>% select("Cell.Type", "cluster") %>% unique()
  
  #Gather identities if multiple per cluster
  ids <- ids %>% group_by(cluster) %>% mutate(label=paste(Cell.Type, collapse=",")) %>% select(-Cell.Type) %>% unique()
  
  #Remove "NA," values 
  ids$label <- str_replace_all(ids$label, "NA,", "") 
  
  #Convert "NA" character to NA
  ids$label[ids$label == "NA"] <- NA
  
  #Reorder
  ids <- ids[order(ids$cluster),]
  
  #Return
  return(ids)
}

#Function for running progeny on object and plotting outputs
runprog <- function(x){
  
  #Create dataframe of clusters
  CellsClusters <- data.frame(Cell = names(Idents(x)),
                              CellType = as.character(Idents(x)),
                              stringsAsFactors = FALSE)
  
  #Run progeny
  x <- progeny(x, scale=FALSE, organism="Human", top=500, perm=1,
               return_assay = TRUE)
  
  ## We can now directly apply Seurat functions in our Progeny scores. 
  ## For instance, we scale the pathway activity scores. 
  x <- Seurat::ScaleData(x, assay = "progeny")
  
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <-
    as.data.frame(t(GetAssayData(x, slot = "scale.data",
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell)
  
  ## We match Progeny scores with the cell clusters.
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
  
  ## We summarize the Progeny scores by cellpopulation
  summarized_progeny_scores <- progeny_scores_df %>%
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  
  #Create dataframe for plotting
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  paletteLength = 100
  myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
  
  progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(summarized_progeny_scores_df)/paletteLength,
                        max(summarized_progeny_scores_df),
                        length.out=floor(paletteLength/2)))
  
  progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14,
                          fontsize_row = 10,
                          color=myColor, breaks = progenyBreaks,
                          main = "PROGENy (500)", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
  
  out <- list(scores=summarized_progeny_scores, heat=progeny_hmap)
  
  return(out)
}

#See Monica covid paper  
compute_stats = function(df, celltype, pathways, conditions, plot = FALSE){
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(viridis))
  
  p = NULL
  df_sub = df %>% filter(condition %in% conditions,
                         Pathway %in% pathways,
                         cell.type %in% celltype)
  
  # p.adjust for number of cells * number of pathways tested
  num_test = nrow(df_sub) * length(pathways)
  
  stats = df_sub %>% group_by(Pathway) %>% 
    nest() %>%
    mutate(wilcox = map(data, function(df){
      stest = wilcox.test(Activity ~ condition, data = df, alternative = 'two.sided')
      broom::tidy(stest) %>%
        mutate(corr_pvalue = p.adjust(p.value, method = 'BH', n = num_test))})) %>%
    select(-data) %>%
    unnest(wilcox) %>%
    ungroup() %>%
    arrange(corr_pvalue) %>% 
    as.data.frame
  
  if(plot){
    p = ggplot(df_sub, aes(x = condition, y = Activity, fill = condition)) +
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), colour = 'lightgrey') +
      theme_classic() + 
      scale_fill_viridis(discrete = TRUE, option = 'viridis') +
      theme(strip.background = element_rect(fill = 'lightgrey'),
            axis.ticks = element_blank()) +
      xlab('') +
      ylab('Pathway activity') +
      facet_wrap( ~ Pathway)
  }
  return(list('stats' = stats, 'p' = p))
}

#---- GSEA function

#Written by Monica Hannani - https://github.com/monicahannani/CovidEpiMap/blob/main/sc_source/sc_source.R

run_gsea = function(bg.genes, stats, category, plot.title = NULL, subcategory = NULL, out.dir = '.', file.prefix, n = 30){
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(fgsea))
  suppressPackageStartupMessages(library(dplyr))
  
  # Fetch geneset
  geneSets = msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
  geneSets = geneSets[geneSets$human_gene_symbol %in% bg.genes,]
  m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)
  
  # Run GSEA
  gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)
  order = order(gsea$padj, decreasing = FALSE)
  
  # Plot
  file.name = paste0(out.dir, '/', file.prefix, '_gsea_', paste0(c(category, subcategory), collapse = '_'))
  
  pdf(file = paste0(file.name, '.pdf'), width = 10, height = 9)
  print(plot_go(gsea.res = gsea,
                gsea.res.order = order, 
                n = n, 
                plot.title = plot.title))
  dev.off()
  
  # Write to file
  write.table(gsea[order, -8], 
              file = paste0(file.name, '.txt'), 
              sep = '\t', 
              row.names = FALSE, 
              quote = FALSE)
}

#---- Plot GO terms from GSEA

plot_go = function(gsea.res, gsea.res.order, plot.title = NULL, n = 20){
  suppressPackageStartupMessages(library(ggplot2))
  
  # Format GOs
  plot.table = head(gsea.res[gsea.res.order,], n = n)
  plot.table$pathway = sub('GO_', '', plot.table$pathway)
  plot.table$pathway = gsub('_', ' ', plot.table$pathway)
  
  p = ggplot(plot.table,
             aes(x = NES, y = pathway)) +
    geom_point(aes(size = NES, color = padj)) +
    theme_bw(base_size = 8) +
    ylab(NULL) +
    ggtitle(plot.title) +
    scale_colour_gradient2(low = 'red', 
                           mid = 'lightgrey', 
                           high = 'blue', 
                           midpoint = 0.05, 
                           limits = c(0,0.1), 
                           oob = scales::squish)
  
  return(p)
}

#Create scatter plot for set of interest
gsea.expr.plot <- function(cluster, genesets, group) {
  c <- metadata %>% filter(seurat_clusters == cluster)
  
  c.df <- data@assays$RNA@data %>% as.data.frame()
  
  c.df <- c.df[,colnames(c.df) %in% row.names(c),] %>% t()
  
  #Create empty dataframe with rownames for averages
  avgs <- data.frame(cellnames=row.names(c.df), row.names = row.names(c.df))
  
  n.sets <- names(genesets)
  
  for (set.n in seq(1, length(genesets))) {
    #select set
    set <- genesets[set.n] %>% unlist()
    setname <- n.sets[set.n]
    
    #Means of set
    means.set <- rowMeans(c.df[,colnames(c.df) %in% set], na.rm=TRUE) %>% as.data.frame()
    colnames(means.set) <- setname
    avgs <- cbind(avgs, means.set)
  }
  
  #Put in long format
  avgs.long <- avgs %>% gather(gset, Value, names(genesets))
  
  #merge metadata info back onto df
  c.merge <- merge(avgs.long, c, all=TRUE, by.x="cellnames", by.y="row.names")
  
  #Plot expression differences
  ggplot(c.merge, aes_string(group, "Value", colour=group)) +
    geom_jitter() +
    facet_wrap(~gset) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5, colour="#616A6B")
    
}

#Create scatter plot for genes of interest
expr.plot <- function(cluster, genes, group) {
  c <- metadata %>% filter(seurat_clusters == cluster)
  
  c.df <- data@assays$RNA@data %>% as.data.frame()
  
  c.df <- c.df[,colnames(c.df) %in% row.names(c),]
  
  #Select genes of interest
  df.goi <- c.df[row.names(c.df) %in% genes,] %>% t()
  
  #merge metadata info back onto df
  c.merge <- merge(df.goi, c, all=TRUE, by="row.names")
  
  #Create short vector of only gene names that were found
  genes.short <- colnames(df.goi)
  
  #Put in long format
  c.merge <- c.merge %>% gather(Gene, Value, genes.short)
  
  #Plot expression differences
  ggplot(c.merge, aes_string(group, "Value", colour=group)) +
    geom_jitter() +
    facet_wrap(~Gene) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5, colour="#616A6B")
  
}