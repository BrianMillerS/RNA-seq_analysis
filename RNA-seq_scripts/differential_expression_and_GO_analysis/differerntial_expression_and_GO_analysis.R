#  Dependencies
options(java.parameters = "- Xmx1024m")
require(Rsamtools)
require(DESeq2)
require(GenomicFeatures)
require(GenomicAlignments)
require(genefilter)
require(pheatmap)
require(RColorBrewer)
library(gplots)
require(ggplot2)
require(sva)
require(openxlsx)
require(GSEABase)
require(GOstats)
require(treeGO)
require(ggrepel)
require(factoextra)
require(gridExtra)
require(plyr)
require(openxlsx)
require(xlsx)


##  FUNCTIONS  ##-----------------------------------------------------------------------------------------------------------------
load_GO_annotations_information <-function(use_BLAST_information = TRUE){
  # load the GO annotations from each specific source
  setwd("/media/pollardlab/POLLARDLAB3/RNAi_Project_T_thermophila_RNA-seq/GO_annotations/output_files")
  
  arabidopsis <- read.csv("arabidopsis_BLAST_derived_GO_universe.csv")
  drosophila <- read.csv("drosophila_BLAST_derived_GO_universe.csv")
  plasmodium <- read.csv("plasmodium_BLAST_derived_GO_universe.csv")
  ciliate <- read.csv("ciliate_org_GO_universe.csv")
  tetramine <- read.csv("tetramine_GO_universe.csv")
  geneOntology <- read.csv("geneOntology_org_GO_universe.csv")
  
  if(use_BLAST_information == TRUE){
    GO_sources_list <- list(arabidopsis= arabidopsis ,drosophila= drosophila ,plasmodium= plasmodium ,ciliate= ciliate ,tetramine= tetramine ,geneOntology= geneOntology)
  } else{
    GO_sources_list <- list(ciliate= ciliate ,tetramine= tetramine ,geneOntology= geneOntology)
  }
  return(GO_sources_list)
}


get_ttherms_with_their_evidence <-function(TTHERM_GO_utilized_init, GO_term, use_BLAST_sources){
  GO_sources_list <- load_GO_annotations_information(use_BLAST_sources)
  try_catch_for_attempt <- function(x){
    return(tryCatch(get_source_str(x) ,error=function(e) NULL))
  }
  
  get_source_str <- function(x){
    row_with_annotation <- x[x$frame.gene_id == TTHERM & x$frame.go_id ==  GO_term, ]
    return(as.character(row_with_annotation$frame.Evidance))
  }
  
  get_evidence_strs <- function(x,index){
    list_to_return <- list()
    return(paste0(index,":",x))
  }
  
  list_to_return_for_each_GO_term <- list()
  list_to_return_for_each_TTHERM <- list()
  for (TTHERM in TTHERM_GO_utilized_init){
    evidence_for_TTHERM_all_sources <- lapply(GO_sources_list, try_catch_for_attempt)
    evidence_for_TTHERM_non_empty_sources <- evidence_for_TTHERM_all_sources[lapply(evidence_for_TTHERM_all_sources,length)>0]
    source_info <- mapply(get_evidence_strs, evidence_for_TTHERM_non_empty_sources, names(evidence_for_TTHERM_non_empty_sources))
    
    if (length(source_info) > 0){
      to_add_to_GO_list <- paste0(TTHERM,"(", paste(source_info, collapse = ",") ,")")
      list_to_return_for_each_GO_term <- c(list_to_return_for_each_GO_term, to_add_to_GO_list)
    }
  }
  if(length(list_to_return_for_each_GO_term) >0 ){
    return(paste(list_to_return_for_each_GO_term, collapse = ", "))
  } else{
    return(NA)
  }
}


get_min_pval_to_use <- function(contrast_input, direction_input, universe_input){
  # This function returns the smallest pvalue one can use when doing a GO analysis with the specific input given
  base::remove(test)
  # make list of pvals to test
  x <- 100
  test_pvals <- c()
  for (i in 1:x){
    zeros <- paste(rep(0,i),collapse="")
    pval_to_test <- as.double(paste0("0.",zeros,"5",collapse=""))
    test_pvals <- c(test_pvals, pval_to_test)
  }
  test_pvals <- base::rev(test_pvals)
  
  # for every pval, if the object is present after running, then stop for it has worked.
  for (pval in test_pvals){
    test_DESEQ2_output <- extract_DE_results(dds_new_with6,
                                             contrast_genotype = contrast_input,
                                             padj_cutoff = pval,
                                             L2FC_magnatude_cutoff = FALSE,
                                             direction = direction_input,
                                             only_TTHERMs = TRUE,
                                             only_names = TRUE)
    
    try(test <- GO_analysis_BH_01(test_DESEQ2_output, universe_input,  direction = "over"), silent = T)
    
    if (base::exists("test")){
      number_of_zeros_in_pval <- as.double(substring(as.character(pval),4,nchar(as.character(pval))))-1  # gets the number of zeros in that smallest pvalue returned, this is done so that make_GO_TTHERM_tables() can use it as input
      return(number_of_zeros_in_pval)
    }
  }
}


generate_tables_for_Cytoscape <- function(universe , dds_used = dds_new_with6, DESEQ2_contrast = "KO_RDF1", DESEQ2_direction = "UP", GO_type = "BP", GO_direction = "over", use_BLAST_sources=TRUE ,max_num_of_zeros = 20){
  # This function generates the node and network csv files that cytoscape will use to make the GO trees
  # the output is a list of two items, the first item is the network df, the second is the node df for Cytoscape to use
  
  # make a list of all the padj values that will be tested
  pvals_to_use <- c()
  for (num_of_zeros in 1:max_num_of_zeros){
    zeros <- paste0(rep(0,num_of_zeros), collapse = "")
    pval_str <- paste0("0.",zeros,"5",collapse = "")
    pvals_to_use <- c(pvals_to_use, as.double(pval_str))
  }
  
  
  ## PART 1: initilize node table with the enriched GO terms, the parent terms will be added later  ##
  node_df <- data.frame(GO_term = as.character(), smallest_L2FC_padj = as.double(), size = as.integer(), largest_L2FC_padj_count = as.integer(), largest_L2FC_padj_count_TTHERM = as.character(), largest_L2FC_padj_TTHERMs_with_annotations = as.character(), has_annotations = as.logical(), stringsAsFactors=FALSE)
  
  for (pval in pvals_to_use){
    TTHERMs_for_GO <- extract_DE_results(dds_used,
                                         contrast_genotype = DESEQ2_contrast,
                                         padj_cutoff = pval,
                                         L2FC_magnatude_cutoff = FALSE,
                                         direction = DESEQ2_direction,
                                         only_TTHERMs = TRUE,
                                         only_names = TRUE)
    
    # run GO
    GO_output <- GO_analysis_BH_01(TTHERMs_for_GO, universe, direction = GO_direction)
    
    # get the correct GO type
    if (GO_type == "BP"){
      GO_output_specific <- GO_output$BP
    } else if (GO_type == "MF"){
      GO_output_specific <- GO_output$MF
    } else if (GO_type == "CC"){
      GO_output_specific <- GO_output$CC
    }
    
    # if there were no results, go to the next pvalue
    if (nrow(GO_output_specific) == 0){
      next
    }
    
    # add GO terms to node table
    for (GO_term in GO_output_specific$GOBPID){
      if (GO_term %in% node_df$GO_term){  # if the GO term is already in the node df
        row_num_to_replace <- which(node_df$GO_term == GO_term)  # get row number in node table that needs to be edited
        
        # replace the padj, as the later padj-cutoffs will be smaller
        node_df[row_num_to_replace, 2] = pval  # replace old pval with smaller one
        
      } else{  # if there is no GO term already present in the node df
        TTHERM_GO_utilized_init <- get_TTHERMs_for_enriched_GO_term_and_children(universe, TTHERMs_for_GO, GO_term)
        TTHERMs_to_return_as_character_init <- paste(TTHERM_GO_utilized_init, collapse =", ")
        
        # if it does then gather those annotations
        GO_sources_list <- load_GO_annotations_information(use_BLAST_sources)
        ttherms_with_evidence <- get_ttherms_with_their_evidence(TTHERM_GO_utilized_init, GO_term, use_BLAST_sources)
        
        # determine whether the GO term has GO Universe annotations
        has_annotations <- ! is.na(ttherms_with_evidence)
        
        # make a row to add the the node df
        GO_term_row <- GO_output_specific[GO_output_specific$GOBPID == GO_term, ]  # get the row from the GO results
        to_add_df <- data.frame(GO_term = as.character(GO_term), smallest_L2FC_padj = as.double(pval), size = as.integer(GO_term_row$Size), largest_L2FC_padj_count = as.integer(GO_term_row$Count), largest_L2FC_padj_count_TTHERM = as.character(TTHERMs_to_return_as_character_init), largest_L2FC_padj_TTHERMs_with_annotations = as.character(ttherms_with_evidence), has_annotations = as.logical(has_annotations),stringsAsFactors=FALSE)  # make row to add
        node_df <- rbind(node_df, to_add_df)  # add row
      }
    }
  }
  
  
  ## PART 2: make the network table  ##
  # we will use ONLY the GO terms that are the furtherest out on the tree to make the network table, this will reduce the number of nudundant network rows that will be made and speed up the process
  most_edge_GO_terms <- c()
  
  for (go_term in node_df$GO_term){
    descendants <- GO_child_assosiations[[go_term]]  # get all descendants of that GO term
    common_GO_terms <- intersect(descendants, node_df$GO_term)  # get the list of enriched GO terms that are descendants to the GO term in question
    
    if (length(common_GO_terms) == 0){  # if are no GO terms further out on the tree
      most_edge_GO_terms <- c(most_edge_GO_terms, go_term)  # add the GO term to the list
    }
  }
  
  # get the most parent GO terms for the correct GO type
  if (GO_type == "BP"){
    GO_root <- "GO:0008150"
  } else if (GO_type == "MF"){
    GO_root <- "GO:0003674"
  } else if (GO_type == "CC"){
    GO_root <- "GO:0005575"
  }
  
  # initilize network df
  network_df <- data.frame(source = as.character(), target = as.character(), stringsAsFactors = FALSE)
  GO_terms_tested <- c()
  
  for (go_Term in most_edge_GO_terms){
    
    current_generation <- c(go_Term)  # define the GO term as the current working generation
    
    while(TRUE){
      do_any_have_parents <- FALSE
      
      # get the parent terms of the current generation, call it upper_generation
      upper_generation <- c()
      for (current_generation_GOterm in current_generation){
        
        # if we have created network rows for this branch, skip it
        if (current_generation_GOterm %in% GO_terms_tested){  
          next
        }
        
        # if go term is the most parent term
        if (current_generation_GOterm == GO_root){
          next
        }
        
        # get the parents of the GO term
        current_generation_GOterm_parents <- getGOParents(current_generation_GOterm)
        current_generation_GOterm_parents_list <- current_generation_GOterm_parents[[current_generation_GOterm]]$Parents  # get list of parent GO terms
        
        # check that here are still parent terms to look at
        if (length(current_generation_GOterm_parents_list) > 0){
          
          # write to network table
          for (parent in current_generation_GOterm_parents_list){
            to_add_df <- data.frame(source = current_generation_GOterm, target = parent, stringsAsFactors = FALSE)
            network_df <- rbind(network_df, to_add_df)
          }
          
          do_any_have_parents <- TRUE
          upper_generation <- c(upper_generation, current_generation_GOterm_parents_list)
        }
      }
      
      # if there are no more parents to search
      if (do_any_have_parents == FALSE){
        break
      }
      
      GO_terms_tested <- c(GO_terms_tested, current_generation)
      current_generation <- upper_generation  # replace current generation with the next generation
    }
  }
  
  # add interaction column to network_df
  network_df$interaction <- rep("pp", nrow(network_df))
  
  
  ## PART 3: add parents to the node table, and add other things to the node table  ##
  network_all_GO_terms <- c(network_df$source, network_df$target)  # get all GO terms from network table
  network_GO_terms_not_in_node_df <- unique(network_all_GO_terms[! network_all_GO_terms %in% node_df$GO_term])  # of those, get the ones that are not covered in the node table, e.i. are parents
  
  # get node size (n) for the parent GO terms
  sizes <- c()
  for (GO_node_term in network_GO_terms_not_in_node_df){
    # get all descendants
    desendants <- unique(GO_child_assosiations[[GO_node_term]])
    desendants <- c(desendants, GO_node_term)
    subsetted_uni_from_desendants <- universe[universe$frame.go_id %in% desendants, ]
    subsetted_uni_from_desendants <- subsetted_uni_from_desendants[! duplicated.data.frame(subsetted_uni_from_desendants), ]
    sizes <- c(sizes, nrow(subsetted_uni_from_desendants))
  }
  NAs <- rep("NA", length(network_GO_terms_not_in_node_df))
  df_network_GO_terms_not_in_node_df <- data.frame(GO_term = network_GO_terms_not_in_node_df, smallest_L2FC_padj = rep(1, length(network_GO_terms_not_in_node_df)), size = sizes, largest_L2FC_padj_count = rep(0, length(network_GO_terms_not_in_node_df)), largest_L2FC_padj_count_TTHERM = NAs, largest_L2FC_padj_TTHERMs_with_annotations = NAs, has_annotations = rep(FALSE, length(network_GO_terms_not_in_node_df)))  # make a df of parent nodes to rbind to the current node table
  node_df <- rbind(node_df, df_network_GO_terms_not_in_node_df)
  
  # add description to node table
  GO_names <- c()
  for (GO_term in node_df$GO_term){
    GO_name <- as.character(annotate::Term(GO_term))
    GO_names <- c(GO_names, GO_name)
  }
  node_df$description <- GO_names
  
  # add column of smallest L2FC in log10 scale
  node_df$smallest_L2FC_padj_log10 <- log10(node_df$smallest_L2FC_padj)
  
  # add column of largest_L2FC_padj_count in log10 scale
  node_df$largest_L2FC_padj_count_log10 <- log10(node_df$largest_L2FC_padj_count)
  node_df <- data.frame(node_df)
  node_df[node_df$largest_L2FC_padj_count_log10 == "-Inf", ]$largest_L2FC_padj_count_log10 <- -1  # replace the "-Inf" with -1
  
  # remove redundant network edges
  network_df <- network_df[! duplicated.data.frame(network_df), ]
  
  return(list(network_df,node_df))
}


compare_GO_tables <- function(GO_table_just_Tet, GO_table_Tet_and_BLASTp){
  # this function compares GO results tables (variable padj tables) and returns four items
  # df_1 is a table including the GO terms that were specific to the GO results from just using the Tetrahymena GO annotations
  # df_2 is a table including the GO terms that were common between df_1 and df_3
  # df_3 is a table including the GO terms that were specific to the GO results from just using the Tetrahymena GO annotations with the BLASTp added annotations
  # is a list of the number of unique GO terms from the three dfs described above as.list(length(unique(df_1$GO_term)), length(unique(df_2$GO_term)), length(unique(df_3$GO_term)))
  
  # get the list of GO terms for the three dfs to make
  GOterms_Tet <- unique(GO_table_just_Tet$GO_term)
  GOterms_TetwithBLAST <- unique(GO_table_Tet_and_BLASTp$GO_term)
  
  GOterms_intersection <- intersect(GOterms_Tet, GOterms_TetwithBLAST)
  
  GOterms_Tet_specific <- GOterms_Tet[! GOterms_Tet %in% GOterms_intersection]
  GOterms_TetwithBLAST_specific <- GOterms_TetwithBLAST[! GOterms_TetwithBLAST %in% GOterms_intersection]
  
  # make dfs to return
  df_intersection <- GO_table_just_Tet[GO_table_just_Tet$GO_term %in% GOterms_intersection, ]  # note how it does not matter which input df is used for this line
  df_Tet_specific <- GO_table_just_Tet[GO_table_just_Tet$GO_term %in% GOterms_Tet_specific, ]
  df_TetwithBLAST_specific <- GO_table_Tet_and_BLASTp[GO_table_Tet_and_BLASTp$GO_term %in% GOterms_TetwithBLAST_specific, ]
  
  # make the list to return
  list_to_return <- c(length(GOterms_Tet_specific), length(GOterms_intersection), length(GOterms_TetwithBLAST_specific))
  
  return(list(df_Tet_specific, df_intersection, df_TetwithBLAST_specific, list_to_return))
}


make_GO_TTHERM_tables <- function(universe = subsetted_GO_Uni_P1, dds_used = dds_new_with6, DESEQ2_contrast = "KO_RSP1", DESEQ2_direction = "UP", GO_type = "BP", GO_direction = "over", max_num_of_zeros = NA){
  # make a list of the DESeq2 pvalues desired to test the padj space by 1 order of magnatude
  pvals_to_use <- c()
  for (num_of_zeros in 1:max_num_of_zeros){
    zeros <- paste0(rep(0,num_of_zeros), collapse = "")
    pval_str <- paste0("0.",zeros,"5",collapse = "")
    pvals_to_use <- c(pvals_to_use, as.double(pval_str))
  }
  
  # take the first pval and initiate the output df
  output_df <- extract_DE_results(dds_used,
                                  contrast_genotype = DESEQ2_contrast,
                                  padj_cutoff = pvals_to_use[1],
                                  L2FC_magnatude_cutoff = FALSE,
                                  direction = DESEQ2_direction,
                                  only_TTHERMs = TRUE,
                                  only_names = FALSE)
  
  # add gene name column and description column
  output_df <- add_gene_name_and_description_columns_to_df(output_df, gene_names_and_descriptions)
  
  for (pval in pvals_to_use){
    
    TTHERMs_for_GO <- extract_DE_results(dds_used,
                                         contrast_genotype = DESEQ2_contrast,
                                         padj_cutoff = pval,
                                         L2FC_magnatude_cutoff = FALSE,
                                         direction = DESEQ2_direction,
                                         only_TTHERMs = TRUE,
                                         only_names = TRUE)
    
    # run GO
    GO_output <- GO_analysis_BH_01(TTHERMs_for_GO, universe, direction = GO_direction)
    
    # get the correct GO type
    if (GO_type == "BP"){
      GO_output_specific <- GO_output$BP
    } else if (GO_type == "MF"){
      GO_output_specific <- GO_output$MF
    } else if (GO_type == "CC"){
      GO_output_specific <- GO_output$CC
    }
    
    # if there were no results, go to the next pvalue
    if (nrow(GO_output_specific) == 0){
      next
    }
    
    for (GO_term in GO_output_specific$GOBPID){
      
      # get description for that GO term
      description <- GO_output_specific[GO_output_specific$GOBPID == GO_term,]$Term
      
      # make output column name for that GO term
      output_col_name <- paste(description, GO_term, sep = ",")
      
      if (output_col_name %in% colnames(output_df)){  # if the GO term already exists, skip it as you have already added all the TTHERMs possible (from a less stringent pvalue)
        next
      }
      
      # get TTHERMs associated with that GO terms enrichment (searching all its children)
      TTHERM_GO_utilized <- get_TTHERMs_for_enriched_GO_term_and_children(universe, TTHERMs_for_GO, GO_term)
      
      # make the column to add, with 1 for representation
      to_add_df <- data.frame(placeholder_name = rep(1, length(TTHERM_GO_utilized)))
      rownames(to_add_df) <- TTHERM_GO_utilized
      names(to_add_df)[names(to_add_df) == "placeholder_name"] <- output_col_name
      
      # add the new column information
      output_df <- merge.data.frame(output_df, to_add_df, by= "row.names", all = TRUE)
      # reset the rownames as the TTHERMs, clean up the column names
      rownames(output_df) <- output_df$Row.names
      output_df <- output_df[ ,2:ncol(output_df)]  # drop column with Row.names
    }
  }
  return(output_df)
}


get_DE_TTHERMs_from_GOterm <- function(GO_Universe, Deseq2_list_input, GOterm){
  # RETURNS a list of the TTHERMs associated with that GO term that are DE (intersected from Deseq2_list_input)
  
  # get all TTHERMs from Universe for that GO term
  all_GO_term_TTHERMs <- GO_Universe[GO_Universe$frame.go_id == GOterm, ]$frame.gene_id  # get list of TTHERMs that are associated with that significant_GO_term
  all_GO_term_TTHERMs <- droplevels.factor(all_GO_term_TTHERMs)  # drop unused levels
  
  # from those TTHERMs only get ones that were significantly differentially expressed
  GO_term_DE_TTHERMs <- intersect(Deseq2_list_input, all_GO_term_TTHERMs)
  
  return(unique(GO_term_DE_TTHERMs))
}


get_TTHERMs_for_enriched_GO_term_and_children <- function(universe, DESEQ2_TTHERM_input, GO_term){
  # get TTHERMs from GO term, input
  output <- get_DE_TTHERMs_from_GOterm(universe, DESEQ2_TTHERM_input, GO_term)
  
  # get TTHERMs from GO terms, children
  for (child_term in GO_child_assosiations[[GO_term]]){
    output <- c(output, get_DE_TTHERMs_from_GOterm(universe, DESEQ2_TTHERM_input, child_term))
  }
  return(unique(output))
}


add_overlapping_TTHERMs_to_custom_annotation_df <- function(input_df, input_counts_filename = "featurecounts_combined_newdat.txt"){
  # For a df with custom annotations as rownames, this adds a column of the TTHERMs that are within their genomic ranges, this uses the genomic
  # location information from the featurecounts file, which comes from the gff 
  # This function also prints out the overlapping TTHERMs and more specific information
  
  # load the counts data and get the TTHERMs
  loaded_count_df <- read.table(input_counts_filename, header=TRUE, row.names=1)
  loaded_count_df <- loaded_count_df[,c(1:4)]  # get the columns that are needed (scf, start, stop, strand)
  TTHERM_specific_count_df <- subset.data.frame(loaded_count_df, ! grepl(":",rownames(loaded_count_df)))  # get all counts rows that are not custom annotations
  
  # set up lists to fill
  lists_of_overlapping_TTHERMs <- c()
  
  for (anot in rownames(input_df)){
    
    anot_row <- loaded_count_df[rownames(loaded_count_df) == anot, ]  # get counts row for every custom annotation
    anot_scf <- as.character(anot_row$Chr)  # get scf for custom annotation
    
    anot_start <- as.numeric(as.character(droplevels(anot_row$Start)))
    anot_end <- as.numeric(as.character(anot_row$End))
    anot_strand <- as.character(anot_row$Strand)
    anot_strand <- paste0("\\",anot_strand)
    print(paste(anot,anot_scf,anot_start,anot_end,anot_strand))
    
    Overlapping_TTHERMs_list <- c()  # define list for this annotation
    TTHERM_samescf_countdf <- TTHERM_specific_count_df[grepl(anot_scf, as.character(TTHERM_specific_count_df$Chr)), ]  # get counts rows for TTHERMs that are on the same scf
    
    TTHERM_samescf_sameStrand_countdf <- TTHERM_samescf_countdf[grepl(anot_strand, as.character(TTHERM_samescf_countdf$Strand)), ]  # of those get the ones that are on the same strand
    
    TTHERM_samescf_countdf_singleexon <- TTHERM_samescf_sameStrand_countdf[! grepl(";", as.character(TTHERM_samescf_sameStrand_countdf$Chr)), ]  # if there is one exon
    for (single_exon_TTHERM in rownames(TTHERM_samescf_countdf_singleexon)){
      single_exon_TTHERM_row <- TTHERM_samescf_countdf_singleexon[single_exon_TTHERM, ]
      single_exon_TTHERM_scf <- as.character(single_exon_TTHERM_row$Chr)
      single_exon_TTHERM_start <- as.numeric(as.character(single_exon_TTHERM_row$Start))
      single_exon_TTHERM_end <- as.numeric(as.character(single_exon_TTHERM_row$End))
      single_exon_TTHERM_strand <- as.character(single_exon_TTHERM_row$Strand)
      
      if(does_it_overlap(anot_start, anot_end, single_exon_TTHERM_start, single_exon_TTHERM_end)){  # if the TTHERM overlaps the annotation at all
        Overlapping_TTHERMs_list <- c(Overlapping_TTHERMs_list, single_exon_TTHERM)  # add the TTHERM to the list
        
        print(paste("Overlapping TTHERM with a single exon:",single_exon_TTHERM, single_exon_TTHERM_scf, as.character(single_exon_TTHERM_start), as.character(single_exon_TTHERM_end),single_exon_TTHERM_strand))
      }
    }
    
    TTHERM_samescf_countdf_multiexon <- TTHERM_samescf_sameStrand_countdf[grepl(";", as.character(TTHERM_samescf_sameStrand_countdf$Chr)), ]  # if there are more than one exon
    
    for (multi_exon_TTHERM in rownames(TTHERM_samescf_countdf_multiexon)){
      multi_exon_TTHERM_row <- TTHERM_samescf_countdf_multiexon[multi_exon_TTHERM, ]
      multi_exon_TTHERM_scf <- as.character(multi_exon_TTHERM_row$Chr)
      multi_exon_TTHERM_start <- as.character(multi_exon_TTHERM_row$Start)
      multi_exon_TTHERM_end <- as.character(multi_exon_TTHERM_row$End)
      multi_exon_TTHERM_strand <- as.character(multi_exon_TTHERM_row$Strand)
      
      # split the multi exon lines into lists
      multi_exon_TTHERM_scf_split <- unlist(strsplit(multi_exon_TTHERM_scf,";"))
      multi_exon_TTHERM_start_split <- unlist(strsplit(multi_exon_TTHERM_start,";"))
      multi_exon_TTHERM_end_split <- unlist(strsplit(multi_exon_TTHERM_end,";"))
      multi_exon_TTHERM_strand_split <- unlist(strsplit(multi_exon_TTHERM_strand,";"))
      
      for (exon_num in 1:length(multi_exon_TTHERM_start_split)){  # for every exon
        exon_num_scf <- as.character(multi_exon_TTHERM_scf_split[[exon_num]])  # get exon scf (some exons from a given gene are on different scaffolds so every exon must be checked)
        
        if (exon_num_scf == anot_scf){  # if the exon is on the same scf
          exon_num_start <- as.numeric(multi_exon_TTHERM_start_split[[exon_num]])  # get exon start
          exon_num_end <- as.numeric(multi_exon_TTHERM_end_split[[exon_num]])  # get exon end
          exon_num_strand <- as.character(multi_exon_TTHERM_strand_split[[exon_num]])  # get exon strand
          
          if(does_it_overlap(anot_start, anot_end, exon_num_start, exon_num_end)){  # if the TTHERM exon overlaps the annotation at all
            Overlapping_TTHERMs_list <- c(Overlapping_TTHERMs_list, multi_exon_TTHERM)  # add the TTHERM to the list
            print(paste("Overlapping TTHERM with multiple exons:", multi_exon_TTHERM, as.character(exon_num), exon_num_scf, as.character(exon_num_start), as.character(exon_num_end), exon_num_strand))
          }
        }
      }
    }
    Overlapping_TTHERMs_list <- paste(unique(Overlapping_TTHERMs_list), collapse = ",")  # remove redundancies (e.g. from multi exon TTHERMs all overlapping the same annotation)
    
    if (Overlapping_TTHERMs_list == ""){  # if blank, replace with NA
      Overlapping_TTHERMs_list = NA
    }
    print(Overlapping_TTHERMs_list)
    print("")
    
    lists_of_overlapping_TTHERMs <- c(lists_of_overlapping_TTHERMs, Overlapping_TTHERMs_list)
  }
  
  input_df$overlapping_TTHERMs <-lists_of_overlapping_TTHERMs
  return(input_df)
}


does_it_overlap <- function(anot_start, anot_end, TTHERM_start, TTHERM_end){
  # returns TRUE or FALSE, TRUE if the THERM overlaps (on either side) or is internal to the annotation(anot)
  # if there is one nt overlap it will return TRUE
  if ((TTHERM_start < anot_start) & (TTHERM_end < anot_start)){  # is the TTHERM fully 5' of the annotation? if so, return FALSE
    return(FALSE)
  } else if ((TTHERM_start > anot_end) & (TTHERM_end > anot_end)){  # is the TTHERM fully 3' of the annotation? if so, return FALSE
    return(FALSE)
  } else {
    return(TRUE)
  }
}


add_genomic_loci_information <- function(input_df, input_counts_filename = "featurecounts_combined_newdat.txt"){
  # adds genomic loci information to df
  loaded_count_data <- read.table(input_counts_filename, header=TRUE, row.names=1)
  
  chr_s <- c()
  start_s <- c()
  end_s <- c()
  strand_s <-c()
  length_s <- c()
  for (anot in rownames(input_df)){
    Chr_to_add <- as.character(loaded_count_data[rownames(loaded_count_data) == anot, ]$Chr)
    Start_to_add <- as.character(loaded_count_data[rownames(loaded_count_data) == anot, ]$Start)
    End_to_add <- as.character(loaded_count_data[rownames(loaded_count_data) == anot, ]$End)
    Strand_to_add <- as.character(loaded_count_data[rownames(loaded_count_data) == anot, ]$Strand)
    length_to_add <- as.character(loaded_count_data[rownames(loaded_count_data) == anot, ]$Length)
    
    chr_s <- c(chr_s, Chr_to_add)
    start_s <- c(start_s, Start_to_add)
    end_s <- c(end_s, End_to_add)
    strand_s <-c(strand_s, Strand_to_add)
    length_s <- c(length_s, length_to_add)
  }
  input_df$Scaffold <- chr_s
  input_df$Start <- start_s
  input_df$End <- end_s
  input_df$Strand <- strand_s
  input_df$Length <- length_s
  
  return(input_df)
}


get_L2FC_col_from_genes_of_interest <- function(input_data_frame, column_name = "column name", genes_of_interest, custom_anotation_only=FALSE){
  # returns a df with the L2FC if the gene is significant (using a df from the /dropbox)
  
  if (custom_anotation_only == TRUE){
    res_GOI <- subset(data.frame(input_data_frame), grepl(":", rownames(input_data_frame)))  # subset results df input to only have custom annotations
  } else {
    res_GOI <- subset(data.frame(input_data_frame), rownames(data.frame(input_data_frame)) %in% genes_of_interest$Gene)  # subset results input df to only contain the genes of interest
  }
  
  is.na(res_GOI$log2FoldChange) <- res_GOI$padj > 0.05
  is.na(res_GOI$log2FoldChange) <- is.na(res_GOI$padj)
  colnames(res_GOI)[2] <- column_name
  return(res_GOI[2])
}


add_has_reads_column_to_df <- function(input_df, input_counts_df, column_name = "has_reads_", genotypes_to_keep = c("P1","F2","N2","F1","SB")){
  # this function adds a column to a data frame (with annotations as rownames) that shows if that annotation had reads mapped to it, TRUE means that there was at least one read mappe to that annotaion
  # genotypes_to_keep can contain values of "P1", "N2", "F2", "F1", or "SB"
  
  # get columns from counts table, subset of input_counts_df for genotypes_to_keep
  counts_col_subset <- input_counts_df[ ,grepl(paste(genotypes_to_keep, collapse="|"), colnames(input_counts_df))]  #get colnames if contain any items in genotypes_to_keep
  
  hasCounts <- c()
  for (anot in rownames(input_df)){
    counts_anot_subset <- counts_col_subset[rownames(counts_col_subset) == anot, ]  # get row from counts
    
    if(sum(counts_anot_subset) == 0){
      hasCounts <-c(hasCounts, FALSE)
    } else {
      hasCounts <-c(hasCounts, NA)
    }
  }
  input_df$placeholder <- hasCounts
  
  names(input_df)[names(input_df) == "placeholder"] <- column_name
  return(input_df)
}


make_L2FC_df <-function(dds_used, genes_of_interest, custom_anotations = FALSE){
  # get constrast (res all results, even if padj = NA (excluded by HB correction))
  DESeq2_res_full_P1 <- extract_DE_results(dds_new_with6, contrast_genotype = "KO_RSP1")
  DESeq2_res_full_N2 <- extract_DE_results(dds_new_with6, contrast_genotype = "KO_RDN2")
  DESeq2_res_full_F2 <- extract_DE_results(dds_new_with6, contrast_genotype = "KO_RDF2")
  DESeq2_res_full_F1 <- extract_DE_results(dds_new_with6, contrast_genotype = "KO_RDF1")
  
  if(custom_anotations == TRUE){
    GOI_df_col_P1 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_P1, column_name = "KO_RSP1", custom_anotation_only=TRUE)
    GOI_df_col_N2 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_N2, column_name = "KO_RDN2", custom_anotation_only=TRUE)
    GOI_df_col_F2 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_F2, column_name = "KO_RDF2", custom_anotation_only=TRUE)
    GOI_df_col_F1 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_F1, column_name = "KO_RDF1", custom_anotation_only=TRUE)
  } else{
    GOI_df_col_P1 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_P1, column_name = "KO_RSP1", genes_of_interest, custom_anotation_only=FALSE)
    GOI_df_col_N2 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_N2, column_name = "KO_RDN2", genes_of_interest, custom_anotation_only=FALSE)
    GOI_df_col_F2 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_F2, column_name = "KO_RDF2", genes_of_interest, custom_anotation_only=FALSE)
    GOI_df_col_F1 <- get_L2FC_col_from_genes_of_interest(DESeq2_res_full_F1, column_name = "KO_RDF1", genes_of_interest, custom_anotation_only=FALSE)
  }
  
  # merge P1 and N2
  combined_GOI_L2FC_df <- merge.data.frame(GOI_df_col_P1, GOI_df_col_N2, by="row.names")
  rownames(combined_GOI_L2FC_df) <- combined_GOI_L2FC_df$Row.names
  combined_GOI_L2FC_df <- combined_GOI_L2FC_df[ ,c(2:ncol(combined_GOI_L2FC_df))]
  # merge, add F2
  combined_GOI_L2FC_df <- merge.data.frame(combined_GOI_L2FC_df, GOI_df_col_F2, by="row.names")
  rownames(combined_GOI_L2FC_df) <- combined_GOI_L2FC_df$Row.names
  combined_GOI_L2FC_df <- combined_GOI_L2FC_df[ ,c(2:ncol(combined_GOI_L2FC_df))]
  # merge, add F1
  combined_GOI_L2FC_df <- merge.data.frame(combined_GOI_L2FC_df, GOI_df_col_F1, by="row.names")
  rownames(combined_GOI_L2FC_df) <- combined_GOI_L2FC_df$Row.names
  combined_GOI_L2FC_df <- combined_GOI_L2FC_df[ ,c(2:ncol(combined_GOI_L2FC_df))]
  
  return(combined_GOI_L2FC_df)
}



get_TTHERM_info_df <- function(){
  # this function is only utilized by add_gene_name_and_description_columns_to_df()
  
  # ensure in the right directory
  setwd("/media/pollardlab/POLLARDLAB3/RNAi_Project_T_thermophila_RNA-seq/differential_expression_and_geneontology_steps")
  
  # import gene names, from ciliage.org scraped file (http://ciliate.org/index.php/show/namedgenes)
  gene_descriptions_df <- read.csv("TGD_List_of_Named_Genes_20180807.txt", header=TRUE, sep="\t")
  gene_descriptions_df <- gene_descriptions_df[ ,c(1,2)]  # remove the gene descriptions from the df, the gff is more complete
  
  # import gene descriptions, from T_thermophila_June2014.gff3 from (http://ciliate.org/index.php/home/downloads)
  gff_file <- read.delim("T_thermophila_June2014.gff3", header=F, comment.char="#")
  
  # get gff lines corresponding to genes
  gff_file_genes <- gff_file[gff_file$V3 == "gene",]
  gff_file_genes_descriptions <- gff_file_genes$V9
  
  # initilize new gff output df
  gff_annotations_df <- NULL
  gff_annotations_df$TTHERM <- c()
  gff_annotations_df$gene_description <- c()
  
  # extract gene descriptions from gff lines
  for (str in gff_file_genes_descriptions){
    splitme <- strsplit(str,";")
    
    #get gene description
    note <- splitme[[1]][3]
    note_split <- strsplit(note, "=")
    function_str <- note_split[[1]][2]
    gff_annotations_df$gene_description <- c(gff_annotations_df$gene_description, function_str)
    
    #get gene name
    name_str <- splitme[[1]][2]
    note_split <- strsplit(name_str, "=")
    name_str <- note_split[[1]][2]
    gff_annotations_df$TTHERM <- c(gff_annotations_df$TTHERM, name_str)
  }
  gff_gene_name_df <- merge.data.frame(gff_annotations_df, gene_descriptions_df, by = "TTHERM", all = TRUE)
  return(gff_gene_name_df)
}


add_gene_name_and_description_columns_to_df <- function(input, gene_info_df){
  # adds a gene description column to all TTHERMs, input is a df with rownames as TTHERMs
  # gene description and name infromation came from get_TTHERM_info_df()
  input_df <- as.data.frame(input)
  
  gene_name_description_df <- gene_info_df
  
  descriptions_to_add <- c()
  names_to_add <- c()
  
  for (row_num in 1:nrow(input_df)){  # for every row in input_df
    row_annotation <- rownames(input_df)[row_num]  # get annotation name
    
    if (grepl(":", row_annotation)){
      descriptions_to_add <- c(descriptions_to_add, NA)
      names_to_add <- c(names_to_add, NA)
      
    } else {
      row_in_gene_name_description_df <- which(gene_name_description_df$TTHERM == row_annotation)  # get row in gene_name_description_df with that TTHERM
      description <- as.character(gene_name_description_df[row_in_gene_name_description_df, 2])  # get the description from that row
      name <- as.character(gene_name_description_df[row_in_gene_name_description_df, 3])  # get the name from that row
      descriptions_to_add <- c(descriptions_to_add, description)
      names_to_add <- c(names_to_add, name)
    }
  }
  
  input_df$description <- descriptions_to_add
  input_df$gene_name <- names_to_add
  return(input_df)
}


GO_analysis_BH_01 <- function(list_of_gene_descriptions_for_analysis, GO_Tet_universe, direction = "over"){
  output <- c()
  
  fdr_alpha=0.1
  go_results_full <- GOanalysis(list_of_gene_descriptions_for_analysis, universe = GO_Tet_universe, pv=1, organism = "Tetrahymena thermophila SB210", testDirection = direction, conditional = FALSE)
  
  # #MF BH correction
  go_results_summary_MF <- summary(go_results_full$MF)
  go_results_summary_MF$Pvalue_BH_adj <- p.adjust(c(go_results_summary_MF$Pvalue), method = "BH")  # add another column for the BH corrected p-value
  output$MF <- subset(go_results_summary_MF, Pvalue_BH_adj<fdr_alpha)
  
  #BP BH correction
  go_results_summary_BP <- summary(go_results_full$BP)
  go_results_summary_BP$Pvalue_BH_adj <- p.adjust(c(go_results_summary_BP$Pvalue), method = "BH")  # add another column for the BH corrected p-value
  output$BP <- subset(go_results_summary_BP, Pvalue_BH_adj<fdr_alpha)
  
  # #CC BH correction
  go_results_summary_CC <- summary(go_results_full$CC)
  go_results_summary_CC$Pvalue_BH_adj <- p.adjust(c(go_results_summary_CC$Pvalue), method = "BH")  # add another column for the BH corrected p-value
  output$CC <- subset(go_results_summary_CC, Pvalue_BH_adj<fdr_alpha)
  return(output)
}


extract_DE_results <- function(dds_input, contrast_genotype = "KO_RSP1", padj_cutoff = FALSE, L2FC_magnatude_cutoff = FALSE, direction = FALSE, only_TTHERMs = FALSE, only_names = FALSE){
  # contrast_genotype can be: "KO_RSP1", "KO_RDF2", "KO_RDF1", or "KO_RDN2"
  # direction can be: "UP", or "DOWN"
  
  # make contrast
  res_subsetted <- as.data.frame(results(dds_input, contrast=c("Strain",contrast_genotype ,"wt")))  # BH FDR = 0.1
  
  # subset results by padj
  if (! padj_cutoff == FALSE){
    res_subsetted <- subset(res_subsetted, padj < padj_cutoff)
  }
  
  # subset results by L2FC magnatude cutoff
  if (! L2FC_magnatude_cutoff == FALSE){
    res_subsetted <- subset(res_subsetted, abs(log2FoldChange) > L2FC_magnatude_cutoff)
  }
  
  # get up or down regulated genes
  if (! direction == FALSE){
    if (direction == "UP"){
      res_subsetted <- subset(res_subsetted, log2FoldChange > 0)
    } else if (direction == "DOWN"){
      res_subsetted <- subset(res_subsetted, log2FoldChange < 0)
    }
  }
  
  # get just TTHERMs (exclude custom annotations)
  if (only_TTHERMs == TRUE){
    res_subsetted <- subset(res_subsetted, ! grepl(":",rownames(res_subsetted)))  # all custom annotations have a ":", so get all that do not have that
  }
  
  # get gene names
  if (only_names == TRUE){
    res_subsetted <- row.names(res_subsetted)
  }
  
  return(res_subsetted)
}


parse_counts_file <- function(counts_filename = "txt", sampleInfo_filename = "csv", use_new_RNA_samples = FALSE){
  # Load sample information csv
  csvfile <- file.path(sampleInfo_filename)
  sampleTable <- read.csv(csvfile) #  read csv
  
  # Load RNAi count data (from featureCounts) and convert to matrix
  countdata <- read.table(counts_filename, header=TRUE, row.names=1)
  countdata <- countdata[,6:ncol(countdata)]  # Remove first five columns (chr, start, end, strand, length)
  
  # define column names for count data
  colnames(countdata) <- sampleTable$Sample
  
  if (use_new_RNA_samples == FALSE){
    #### remove samples from new RNA prep <>
    countdata <- countdata[,c(4,5,7:17)]  # remove col from count date
    sampleTable <- sampleTable[c(4,5,7:17), ]  # remove rows from sample table
    sampleTable <- droplevels(sampleTable)  # drop unused levels, doesn't matter really, but it prevents DESeqq2 from throwing a warning
    print(sampleTable$Sample  == colnames(countdata))  # double check
  }
  
  #convert to matrix
  countdata <- as.matrix(countdata)
  return(list(countdata,sampleTable))
}

make_dds <- function(countdata, sampleTable){
  
  # run Deseq2
  coldata <- data.frame(row.names=colnames(countdata), factor(sampleTable$Sample_Match), sampleTable$Strain)
  colnames(coldata) <- c("Sample_Match", "Strain")
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~ Sample_Match + Strain)
  dds <- DESeq(dds)
  
  return(dds)
}


subset_GO_universe <- function(input_GO_universe, input_counts_df, genotypes_to_keep = c("P1", "SB")){
  # this function returns a subsetted GO universe that only contains TTHERMs associations where the TTHERM is expressed in the genotypes_to_keep genotypes
  # genotypes_to_keep can contain values of "P1", "N2", "F2", "F1", or "SB"
  
  # get columns, subset of input_counts_df for genotypes_to_keep
  counts_subset <- input_counts_df[ ,grepl(paste(genotypes_to_keep, collapse="|"), colnames(input_counts_df))]  #get colnames if contain any items in genotypes_to_keep
  
  # get rows that have at least one read across all samples in counts_subset
  counts_subset_expressedGenes <- counts_subset[rowSums(counts_subset) != 0, ]
  
  # only take the TTHERMs (exclude custom annotations (things that have a ":"))
  counts_subset_expressedGenes_TTHERMs <- subset(counts_subset_expressedGenes, ! grepl(":",rownames(counts_subset_expressedGenes)))  # all custom annotations have a ":", so get all that do not have that
  
  # get the names of the TTHERMs
  TTHERMs_to_subset_GOuni_with <- rownames(counts_subset_expressedGenes_TTHERMs)
  
  # subset GO universe to onoly have contain mappings with those TTHERMs
  output_GO_universe <- input_GO_universe[input_GO_universe$frame.gene_id %in% TTHERMs_to_subset_GOuni_with, ]
  
  # print outputs
  print(paste("From your genotype selections:",as.character(length(unique(input_GO_universe$frame.gene_id)) - length(unique(output_GO_universe$frame.gene_id))),"of the", as.character(length(unique(input_GO_universe$frame.gene_id))), "TTHERMs had no expression and were dropped form the GO Universe"))
  print(paste("From your gentype selections:", as.character(nrow(input_GO_universe) - nrow(output_GO_universe)), "of the", as.character(nrow(input_GO_universe)), "GO mappings corresponded to TTHERMs with no expression, and were dropped from the GO Universe"))
  
  return(output_GO_universe)
}


make_GO_Universe_file_for_cytoscape <- function(input_GO_Universe, output_filename = "custom_GO_annotations_RSP1vsSB210.txt"){
  # initiate file
  outfile <- file(output_filename)
  
  # get lines to write
  lines_to_write <- paste(input_GO_Universe$frame.gene_id, substr(input_GO_Universe$frame.go_id, 4,10), sep = " = ")
  writeLines(c("(species=Tetrahymena thermophila SB210)(type=All)(curator=GO)", lines_to_write), outfile)
  
  close(outfile)
}


make_GO_BH_padj_table <-function(universe, dds_used, DESEQ2_contrast = "KO_RSP1", DESEQ2_direction = "UP", GO_type = "BP", GO_direction = "over", max_num_of_zeros = FALSE){
  # make a list of the DESeq2 pvalues desire
  pvals_to_use <- c()
  for (num_of_zeros in 1:max_num_of_zeros){
    zeros <- paste0(rep(0,num_of_zeros), collapse = "")
    pval_str <- paste0("0.",zeros,"5",collapse = "")
    pvals_to_use <- c(pvals_to_use, as.double(pval_str))
  }
  
  # take the first pval and initiate the output df
  first_pval_DESEQ2_output <- extract_DE_results(dds_used,
                                                 contrast_genotype = DESEQ2_contrast, 
                                                 padj_cutoff = pvals_to_use[1],
                                                 L2FC_magnatude_cutoff = FALSE,
                                                 direction = DESEQ2_direction,
                                                 only_TTHERMs = TRUE,
                                                 only_names = TRUE)
  
  
  first_pval_GO_output <- GO_analysis_BH_01(first_pval_DESEQ2_output, universe,  direction = "over")
  
  first_pval_GO_output$BP$toReturn <- paste(first_pval_GO_output$BP$Count, formatC(first_pval_GO_output$BP$Pvalue_BH_adj, format = "e", digits = 2), sep = ", ")
  output_df <- data.frame(GO_term=first_pval_GO_output$BP$GOBPID, GO_Description=first_pval_GO_output$BP$Term, pval0.05=first_pval_GO_output$BP$toReturn)
  
  
  for (pval in pvals_to_use[2:length(pvals_to_use)]){
    
    next_pval_DESEQ2_output <- extract_DE_results(dds_used,
                                                  contrast_genotype = DESEQ2_contrast, 
                                                  padj_cutoff = pval,
                                                  L2FC_magnatude_cutoff = FALSE,
                                                  direction = DESEQ2_direction,
                                                  only_TTHERMs = TRUE,
                                                  only_names = TRUE)
    
    next_pval_GO_output <- GO_analysis_BH_01(next_pval_DESEQ2_output, universe,  direction = "over")
    
    if (nrow(next_pval_GO_output$BP) == 0){
      next
    }
    
    new_GO_terms_table <- subset(next_pval_GO_output$BP, ! GOBPID %in% output_df$GO_term)  # get subset of next output for the GO terms not currently in the output table
    add_to_output_df <- data.frame(GO_term=new_GO_terms_table$GOBPID, GO_Description=new_GO_terms_table$Term)  # make the subset into a df
    
    output_df <- rbind.fill(output_df, add_to_output_df)  # add new rows, GO terms and descriptions
    
    next_pval_GO_output$BP$toReturn <- paste(next_pval_GO_output$BP$Count, formatC(next_pval_GO_output$BP$Pvalue_BH_adj, format = "e", digits = 2), sep = ", ")  # add col to return
    
    new_col_to_add <- data.frame(GO_term=next_pval_GO_output$BP$GOBPID, placeholder_name=next_pval_GO_output$BP$toReturn)  # get new BH pvals to add
    
    output_df <- merge.data.frame(output_df, new_col_to_add, by="GO_term", all = TRUE) # add new columns, GO terms
    names(output_df)[names(output_df) == "placeholder_name"] <- pval
  }
  
  return(output_df)
}

