input_dir <- '/home/lisbet/Data/2021_GinoRoche/input/'
data_dir <- '/home/lisbet/Data/2021_GinoRoche/Data/'
results_dir <- '/home/lisbet/Data/2021_GinoRoche/results/'

# Data tweaking of CellPhoneDB input files
#------------------------------------------------------------------------------------------
# read in the data table
# symbols and entrez
phdb_input <- read.csv(paste0(input_dir, 'gene_input.csv'), header = T)
# read in the interaction and protein databases extracted from cellphone db
int_df <- read.csv(paste0(input_dir, 'interaction_input.csv'), header = T)
protein_df <- read.csv(paste0(input_dir, 'protein_curated.csv'), header = T)

# manipulate the hgnc symbols
protein_df <- protein_df %>%
  mutate(symbols_modified = gsub("\\_.*", "", protein_df$protein_name))

protein_df %>%
  filter(receptor == 'True') %>%
  dim()


receptor_vec <- protein_df %>%
  filter(receptor == 'True') %>%
  select(symbols_modified)
ligand_vec <- protein_df %>%
  filter(receptor == 'False') %>%
  select(symbols_modified)

# extract protein_name_a and _b from int_df -> why actually

int_df_reduced <- int_df[, c('partner_a', 'partner_b')]

# set up of the final df: Ligand -> A, Receptor -> B
# manipulate the int_df_reduced such that we identify receptors in partner a
# if so then add receptor_in_a as TRUE, else false
# do the call
receptors <- protein_df %>% filter(receptor == 'True')

int_df_reduced <- int_df_reduced %>%
  mutate(HGNC_partner_a = protein_df[match(int_df_reduced$partner_a, protein_df$uniprot), 'symbols_modified']) %>%
  mutate(HGNC_partner_b = protein_df[match(int_df_reduced$partner_b, protein_df$uniprot), 'symbols_modified']) %>%
  mutate(partner_a_hgnc = coalesce(HGNC_partner_a, partner_a)) %>%
            mutate(partner_b_hgnc = coalesce(HGNC_partner_b, partner_b))

int_df_reduced <- int_df_reduced %>%
  plyr::mutate(receptor_in_a = ifelse((grepl('receptor', int_df_reduced$partner_a_hgnc) == TRUE), TRUE, FALSE)) %>%
  mutate(is_receptor_from_proteinDB = ifelse((grepl(paste(receptors$symbols_modified, collapse = '|'), int_df_reduced$partner_a_hgnc) == TRUE), TRUE, FALSE))


# 87 entries must be reversed
wrong_receptors <- int_df_reduced %>%
  filter(receptor_in_a == TRUE | is_receptor_from_proteinDB == TRUE) %>%
  select('partner_b_hgnc', 'partner_a_hgnc') %>%
  set_colnames(c('partner_a_hgnc', 'partner_b_hgnc'))

# first call the pairs that are correct
# then add the reversed wrong entries
int_df_reordered <- int_df_reduced %>%
  filter(receptor_in_a == FALSE & is_receptor_from_proteinDB == FALSE) %>%
  select('partner_a_hgnc', 'partner_b_hgnc') %>%
  rbind(wrong_receptors)

all_ligands <- data.frame(uniprot = protein_df[match(int_df_reordered$partner_a, protein_df$uniprot), 'symbols_modified'],
                          hgnc = int_df_reordered$partner_a) %>%
  dplyr::mutate(combined = coalesce(uniprot, hgnc)) %>%
  dplyr::pull(combined)

all_receptors <- data.frame(uniprot = protein_df[match(int_df_reordered$partner_b, protein_df$uniprot), 'symbols_modified'],
                            hgnc = int_df_reordered$partner_b) %>%
  dplyr::mutate(combined = coalesce(uniprot, hgnc)) %>%
  dplyr::pull(combined)

get_uniprot <- function(hgnc_input, reference_df) {
  uniprot <- protein_df[reference_df$symbols_modified %in% hgnc_input,] %>% select(c(symbols_modified, uniprot))
  # uniprot <- subset_df[match(hgnc_input, subset_df$symbols_modified),]
  return(uniprot)
}

run_CellPhoneDB_across_all <- function(regulation) {
  sender_names <- c('CD4', 'CD8')
  receiver_names <- paste0('CD', c(4, 8, 16, 19))
  disease_states <- c('HD', 'pSS', 'SLE')
  
  for (i in disease_states) {
    for (k in sender_names) {
      
      
      sender_name_chr <- (paste0(k, '_', i))
      cat(paste('The sender data is:', sender_name_chr, '\n'))
      
      sender_ids = df_gmt %>%
        dplyr::filter(Condition == sender_name_chr) %>%
        dplyr::pull(Sample_ID) %>%
        unique()
      
      
      for (j in receiver_names) {
        
        receiver_name_chr <- (paste0(j, '_', i))
        cat(paste('The receiver data is:', receiver_name_chr, '\n'))
        
        receiver_ids = df_gmt %>%
          dplyr::filter(Condition == receiver_name_chr) %>%
          dplyr::pull(Sample_ID) %>%
          unique()
        
        
        if (regulation == 'up') {
          
          sender_median = df_gmt %>%
            dplyr::filter(Condition == sender_name_chr) %>%
            #  dplyr::filter(Condition == 'CD4_HD') %>%
            dplyr::summarise(th = quantile(ExpressionValue, 0.9))
          
          expressed_genes_sender = t_gmt[sender_ids, colnames(t_gmt) %in% all_ligands] %>%
            apply(2, function(x) { mean(x) }) %>%
            .[. > sender_median$th] %>%
            names()
          
          receiver_median = df_gmt %>%
            dplyr::filter(Condition == receiver_name_chr) %>%
            dplyr::summarise(th = quantile(ExpressionValue, 0.1))
          
          expressed_genes_receiver = t_gmt[receiver_ids, colnames(t_gmt) %in% all_receptors] %>%
            apply(2, function(x) { mean(x) }) %>%
            .[. > receiver_median$th] %>%
            names()
          
          
          cat(paste('Number of genes in sender:', length(expressed_genes_sender), '\n'))
          cat(paste('Number of genes in receiver:', length(expressed_genes_receiver), '\n'))
        } else if (regulation == 'down') {
          
          sender_median = df_gmt %>%
            dplyr::filter(Condition == sender_name_chr) %>%
            #  dplyr::filter(Condition == 'CD4_HD') %>%
            dplyr::summarise(th = quantile(ExpressionValue, 0.1))
          
          expressed_genes_sender = t_gmt[sender_ids, colnames(t_gmt) %in% all_ligands] %>%
            apply(2, function(x) { mean(x) }) %>%
            .[. < sender_median$th] %>%
            names()
          
          receiver_median = df_gmt %>%
            dplyr::filter(Condition == receiver_name_chr) %>%
            dplyr::summarise(th = quantile(ExpressionValue, 0.9))
          
          expressed_genes_receiver = t_gmt[receiver_ids, colnames(t_gmt) %in% all_receptors] %>%
            apply(2, function(x) { mean(x) }) %>%
            .[. < receiver_median$th] %>%
            names()
          
          cat(paste('Number of genes in sender:', length(expressed_genes_sender), '\n'))
          cat(paste('Number of genes in receiver:', length(expressed_genes_receiver), '\n'))
          
        }
        
        
        # Start of cellphoneDB workflow
        # Query for ligands
        uniprot_sender <- get_uniprot(hgnc_input = expressed_genes_sender, reference_df = protein_df)
        # uniprot_sender <- uniprot_sender %>% rownames_to_column() %>% set_colnames(c('rank', 'symbols_modified', 'uniprot'))
        
        df_ligands <- int_df_reordered[
          grepl(pattern = paste(uniprot_sender$uniprot, collapse = '|'), int_df_reordered$partner_a) |
            grepl(pattern = paste(uniprot_sender$symbols_modified, collapse = '|'), int_df_reordered$partner_a),]
        
        # Query for receptors
        
        uniprot_receiver <- get_uniprot(hgnc_input = expressed_genes_receiver, reference_df = protein_df)
        
        df_final <- df_ligands[
          grepl(pattern = paste(uniprot_receiver$uniprot, collapse = '|'), df_ligands$partner_b) |
            grepl(pattern = paste(uniprot_receiver$symbols_modified, collapse = '|'), df_ligands$partner_b),]
        
        df_final <- df_final %>%
          mutate(HGNC_partner_a = protein_df[match(df_final$partner_a, protein_df$uniprot), 'symbols_modified']) %>%
          mutate(HGNC_partner_b = protein_df[match(df_final$partner_b, protein_df$uniprot), 'symbols_modified'])
        
        if (any(is.na(df_final)) != FALSE) {
          print('Combining HGCN symbols first')
          df_final <- df_final %>%
            mutate(combine_a = coalesce(HGNC_partner_a, partner_a)) %>%
            mutate(combine_b = coalesce(HGNC_partner_b, partner_b)) %>%
            #mutate(rank_a = sender_df[match(df_final$HGNC_partner_a, sender_df$HGNC), 'rank']) %>%
            #mutate(rank_b = receiver_df[match(df_final$HGNC_partner_b, receiver_df$HGNC), 'rank']) %>%
            select(c('partner_a',
                     'partner_b', 'combine_a', 'combine_b')) %>%
            set_colnames(c('partner_a', 'partner_b', 'HGNC_partner_a', 'HGNC_partner_b'))
        }
        df_final %>%
          select(c('HGNC_partner_a', 'HGNC_partner_b')) %>%
          write.csv(file = paste0('/home/lisbet/Data/2021_GinoRoche/Data/M_001_CellPhoneDB_ExprData_', sender_name_chr, '_in_', receiver_name_chr, '.csv'))
      }
    }
  }
}


generate_count_info <- function(which_sender, which_receiver) {
  df_big <- com_df %>%
    filter(Sender_cell == which_sender & Receiver_cell == which_receiver) %>%
    group_by(Ligand, Receptor) %>%
    droplevels() %>%
    dplyr::mutate(count = factor(n(), levels = c('1', '2', '3'))) %>%
    dplyr::mutate(Ligand = factor(Ligand, levels = sort(levels(Ligand)))) %>%
    dplyr::mutate(Receptor = factor(Receptor, levels = sort(levels(Receptor)))) %>%
    dplyr::mutate(combination = paste0(Ligand, '_', Receptor))
  
  df_twos <- com_df %>%
    dplyr::filter(Sender_cell == which_sender & Receiver_cell == which_receiver) %>%
    dplyr::group_by(Ligand, Receptor) %>%
    droplevels() %>%
    dplyr::mutate(count = factor(n(), levels = c('1', '2', '3'))) %>%
    dplyr::mutate(Ligand = factor(Ligand, levels = sort(levels(Ligand)))) %>%
    dplyr::mutate(Receptor = factor(Receptor, levels = sort(levels(Receptor)))) %>%
    dplyr::filter(count == 2) %>%
    dplyr::mutate(combination = paste0(Ligand, '_', Receptor))
  
  fill_vect <- c()
  for (i in 1:dim(df_big)[1]) {
    if ((df_big$count[i] == 2) == TRUE) {
      df_sub <- df_twos[df_twos$combination == df_big$combination[i],]
      cul_comb <- paste(df_sub$Disease_state, collapse = "--")
      fill_vect <- c(fill_vect, cul_comb)
    } else if ((df_big$count[i] == 1) == TRUE) {
      cul_comb <- df_big$Disease_state[i]
      fill_vect <- c(fill_vect, cul_comb)
    } else if ((df_big$count[i] == 3) == TRUE) {
      cul_comb <- 'HD--SLE--pSS'
      fill_vect <- c(fill_vect, cul_comb)
    } else {
      cul_comb <- NA
      fill_vect <- c(fill_vect, cul_comb)
    }
  }
  
  df_big$Visualization <- factor(fill_vect)
  return(df_big)
}

plot_heatmap <- function(df_visualization, which_sender, which_receiver) {
  p <- ggplot(df_visualization, aes(Receptor, Ligand, fill = Visualization)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 23)) +
    ggtitle(paste('Ligand receptor interactions from', which_sender, 'to', which_receiver, 'cells'))
  return(p)
}

run_visualization <- function(sender_name, receiver_name, levels_disease) {
  data_df <- generate_count_info(which_sender = sender_name,
                                 which_receiver = receiver_name)
  
  
  data_df$Visualization <- factor(data_df$Visualization,
                                  levels = levels_disease)
  
  p <- plot_heatmap(df_visualization = data_df,
                    which_sender = sender_name,
                    which_receiver = receiver_name)
  
  
  p <- p + scale_fill_manual('Occurence', values = c("HD" = "#12B12A",
                                                     "SLE" = "#1212B1",
                                                     "pSS" = "#9212B1",
                                                     "HD--SLE" = "#07C9CF",
                                                     "HD--pSS" = "#F976F9",
                                                     "pSS--SLE" = "#CF0707",
                                                     "HD--SLE--pSS" = "#989898"))
  return(p)
}

# Function to run the ligand-receptor analysis per donor (for upregulated ligands)
run_CellPhoneDB_SampleWise <- function() {
  
  sender_names <- c('CD4', 'CD8')
  receiver_names <- paste0('CD', c(4, 8, 16, 19))
  
  
  all_donors <- unique(df_gmt$Donor)
  
  all_conditions <- unique(df_gmt$Condition)
  ################################################################################################
  # The threshold information on ligand site #####################################################
  ################################################################################################
  
  threshold_info <- c()
  for (cond in all_conditions) {
    sender_median = df_gmt %>%
      dplyr::filter(Condition == cond) %>%
      #  dplyr::filter(Condition == 'CD4_HD') %>%
      dplyr::summarise(th = quantile(ExpressionValue, 0.9))
    
    threshold_info <- c(threshold_info, sender_median$th)
  }
  names(threshold_info) <- all_conditions
  
  ################################################################################################
  
  for (n in all_donors) {
    for (s in sender_names) {
      which_donor_sender <- paste0(s, '_', n)
      which_th_name <- gsub(which_donor_sender, pattern = '_\\d+$', replacement = '')
      threshold <- threshold_info[which_th_name]
      cat(paste('used threshold:', threshold, '\n'))
      # expressed_genes_sender = t_gmt[which_donor_sender, colnames(t_gmt) %in% all_ligands] %>%
      #   #apply(2, function(x) { mean(x) }) %>%
      #   .[. >= threshold] %>%
      #   names()
      
      expressed_genes_sender = t_gmt[which_donor_sender,] %>%
        #apply(2, function(x) { mean(x) }) %>%
        .[. > threshold] %>%
        names()
      
      for (r in receiver_names) {
        which_donor_receiver <- paste0(r, '_', n)
        cat(paste('Comparison', which_donor_sender, 'vs', which_donor_receiver, '\n'))
        
        # expressed_genes_receiver = t_gmt[which_donor_receiver, colnames(t_gmt) %in% all_receptors] %>%
        #   #apply(2, function(x) { mean(x) }) %>%
        #   .[. != 0] %>%
        #   names()
        
        expressed_genes_receiver = t_gmt[which_donor_receiver,] %>%
          #apply(2, function(x) { mean(x) }) %>%
          .[. != 0] %>%
          names()
        
        cat(paste('Number of genes in sender:', length(expressed_genes_sender), '\n'))
        cat(paste('Number of genes in receiver:', length(expressed_genes_receiver), '\n'))
        
        
        # Start of cellphoneDB workflow
        # Query for ligands
        uniprot_sender <- get_uniprot(hgnc_input = expressed_genes_sender, reference_df = protein_df)
        # uniprot_sender <- uniprot_sender %>% rownames_to_column() %>% set_colnames(c('rank', 'symbols_modified', 'uniprot'))
        
        df_ligands <- int_df_reordered[
          grepl(pattern = paste(uniprot_sender$uniprot, collapse = '|'), int_df_reordered$partner_a_hgnc) |
            grepl(pattern = paste(uniprot_sender$symbols_modified, collapse = '|'), int_df_reordered$partner_b_hgnc),]
        
        # Query for receptors
        
        uniprot_receiver <- get_uniprot(hgnc_input = expressed_genes_receiver, reference_df = protein_df)
        
        df_final <- df_ligands[
          grepl(pattern = paste(uniprot_receiver$uniprot, collapse = '|'), df_ligands$partner_a_hgnc) |
            grepl(pattern = paste(uniprot_receiver$symbols_modified, collapse = '|'), df_ligands$partner_b_hgnc),]
        
        df_final <- df_final %>%
          mutate(HGNC_partner_a = protein_df[match(df_final$partner_a_hgnc, protein_df$uniprot), 'symbols_modified']) %>%
          mutate(HGNC_partner_b = protein_df[match(df_final$partner_b_hgnc, protein_df$uniprot), 'symbols_modified'])
        
        if (any(is.na(df_final)) != FALSE) {
          print('Combining HGCN symbols first')
          df_final <- df_final %>%
            mutate(combine_a = coalesce(HGNC_partner_a, partner_a_hgnc)) %>%
            mutate(combine_b = coalesce(HGNC_partner_b, partner_b_hgnc)) %>%
            #mutate(rank_a = sender_df[match(df_final$HGNC_partner_a, sender_df$HGNC), 'rank']) %>%
            #mutate(rank_b = receiver_df[match(df_final$HGNC_partner_b, receiver_df$HGNC), 'rank']) %>%
            select(c('partner_a_hgnc',
                     'partner_b_hgnc', 'combine_a', 'combine_b')) %>%
            set_colnames(c('partner_a_hgnc', 'partner_b_hgnc', 'HGNC_partner_a', 'HGNC_partner_b'))
        }
        df_final %>%
          select(c('HGNC_partner_a', 'HGNC_partner_b')) %>%
          write.csv(file = paste0('/home/lisbet/Data/2021_GinoRoche/Data/M_stSt_001_CellPhoneDB_ExprData_', which_donor_sender, '-in-', which_donor_receiver, '.csv'))
      }
    }
  }
}

# These functions were used to build the dataframe, calculate the fischer`s exact test and visualize the result
#-----------------------------------------------------------------------------------------------------------------


# Create the input used for the fischer`s exact test
create_dataframe <- function(sender_contrast, receiver_contrast, disease_contrast) {
  cat(paste('Building count table for', sender_contrast, 'agains', receiver_contrast, 'in', disease_contrast[1], 'vs', disease_contrast[2], '.----------\n'))
  df_a <- com_df %>%
    filter(Sender_cell == sender_contrast &
             Receiver_cell == receiver_contrast &
             Disease_state == disease_contrast[1]) %>%
    select(c('Ligand', 'Receptor')) %>%
    table() %>%
    melt()
  
  df_b <- com_df %>%
    filter(Sender_cell == sender_contrast &
             Receiver_cell == receiver_contrast &
             Disease_state == disease_contrast[2]) %>%
    select(c('Ligand', 'Receptor')) %>%
    table() %>%
    melt()
  
  comp_df <- merge(x = df_a, y = df_b, by = c('Ligand', 'Receptor'))
  
  comp_df <- comp_df %>%
    set_rownames(paste0(comp_df$Ligand, '_', comp_df$Receptor)) %>%
    select(c('value.x', 'value.y')) %>%
    set_colnames(c('Cond1', 'Cond2'))
  row_sub <- apply(comp_df, 1, function(row) !all(row == 0))
  comp_df_f <- comp_df[row_sub,]
  return(comp_df_f)
}

# Actually do the fischer`s exact test
run_Fishers_contrast <- function(comp_df, sender_contrast, receiver_contrast, disease_contrast) {
  # preprocessing the data such that we can generate the contigency table for the Fishers test
  # -----------------------------------------------------------------------------------------------
  cat(paste('The sender cell is', sender_contrast, ', the receiver', receiver_contrast, '.----------\n'))
  cat(paste('The contrast to be investigated is', disease_contrast[1], 'vs', disease_contrast[2], '.----------\n'))
  
  how_many_samples <- table(gsub(colnames(gmt), pattern = "\\_[0-9]*$", replacement = ''))
  
  contrast_sender_a <- paste0(sender_contrast, '_', disease_contrast[1])
  contrast_receiver_a <- paste0(receiver_contrast, '_', disease_contrast[1])
  
  contrast_sender_b <- paste0(sender_contrast, '_', disease_contrast[2])
  contrast_receiver_b <- paste0(receiver_contrast, '_', disease_contrast[2])
  
  n_sender_a <- how_many_samples[contrast_sender_a]
  n_receiver_a <- how_many_samples[contrast_receiver_a]
  
  n_sender_b <- how_many_samples[contrast_sender_b]
  n_receiver_b <- how_many_samples[contrast_receiver_b]
  # -----------------------------------------------------------------------------------------------
  
  p_value_df <- data.frame('Interaction' = character(),
                           'pvalue' = numeric(),
                           'adj.pvalue' = numeric(),
                           'direction' = character(),
                           'Sender' = character(),
                           'Receiver' = character(),
                           'Contrast' = character())
  # calculate the statistics
  interactions <- rownames(comp_df)
  for (i in 1:dim(comp_df)[1]) {
    contrast_a <- c(comp_df$Cond1[i], ((n_sender_a) - comp_df$Cond1[i]))
    contrast_b <- c(comp_df$Cond2[i], ((n_sender_b) - comp_df$Cond2[i]))
    
    df <- data.frame(Cond1 = contrast_a, Cond2 = contrast_b)
    fischer_result <- pairwise_fisher_test(t(df), alternative = 'greater', p.adjust.method = 'BH')
    p_value_vect_greater <- data.frame('Interaction' = interactions[i],
                                       'pvalue' = fischer_result$p,
                                       'adj.pvalue' = fischer_result$p.adj,
                                       'direction' = 'greater',
                                       'Sender' = sender_contrast,
                                       'Receiver' = receiver_contrast,
                                       'Contrast' = paste0(disease_contrast, collapse = '_'))
    
    fischer_result <- pairwise_fisher_test(t(df), alternative = 'less', p.adjust.method = 'BH')
    p_value_vect_less <- data.frame('Interaction' = interactions[i],
                                    'pvalue' = fischer_result$p,
                                    'adj.pvalue' = fischer_result$p.adj,
                                    'direction' = 'less',
                                    'Sender' = sender_contrast,
                                    'Receiver' = receiver_contrast,
                                    'Contrast' = paste0(disease_contrast, collapse = '_'))
    p_value_df <- p_value_df %>%
      bind_rows(p_value_vect_greater) %>%
      bind_rows(p_value_vect_less)
    
  }
  return(p_value_df)
}

# Data visualization of up and down regulated signals in one ggplot object
plot_contrast_high_low <- function(results_df, which_disease, threshold_p) {
  
  results_df <- results_df %>%
    dplyr::mutate(Receiver = factor(results_df$Receiver, levels =c('CD19', 'CD16', 'CD8', 'CD4'))) %>%
    dplyr::mutate(Sender = factor(results_df$Sender, levels =c('CD4', 'CD8')))
  #plot for sign. lower signals
  # -----------------------------------------------------------------------------------------------
  p <- results_df %>%
    dplyr::mutate(Interaction = factor(Interaction)) %>%
    #dplyr::mutate(Sender = factor(results_df$Sender)) %>%
    #dplyr::mutate(Receiver = factor(results_df$Receiver, levels =c('CD4', 'CD8', 'CD16', 'CD19'))) %>%
    dplyr::mutate(direction = factor(results_df$direction)) %>%
    dplyr::mutate(Interaction = factor(Interaction, levels = sort(levels(Interaction)))) %>%
    mutate('pvalue filtered' = ifelse(adj.pvalue <= threshold_p, adj.pvalue, NA)) %>%
    drop_na() %>%
    filter(direction == 'less') %>%
    ggplot(aes(Receiver, Interaction, fill = `pvalue filtered`)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(color = "black", fill = NA)) +
    geom_tile() +
    scale_fill_stepsn("DN in disease", colours = c('#0000FF','#9191FF', '#DADAFF' , '#FFFFFF'), breaks = c(0,0.001,0.05,0.1), limits = c(0, 0.11), na.value = 'white') +
    #labs(subtitle = paste('CCC signals sig. different in ', which_disease, 'compared to healthy'), caption = 'Fishers test \n If adjusted p-value > 0.1 then not visualized') +
    facet_grid(rows = vars(Sender), scales = "free") +
    coord_flip() +
    theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1))
  
  # plot sign higher interactions
  # -----------------------------------------------------------------------------------------------
  # first filter the data frame for higher interactions
  greater_df <- results_df %>%
    dplyr::mutate(Interaction = factor(Interaction)) %>%
    #dplyr::mutate(Sender = factor(results_df$Sender)) %>%
    #dplyr::mutate(Receiver = factor(results_df$Receiver, levels =c('CD4', 'CD8', 'CD16', 'CD19'))) %>%
    dplyr::mutate(direction = factor(results_df$direction)) %>%
    dplyr::mutate(Interaction = factor(Interaction, levels = sort(levels(Interaction)))) %>%
    mutate('pvalue filtered' = ifelse(adj.pvalue <= threshold_p, adj.pvalue, NA)) %>%
    drop_na() %>%
    filter(direction == 'greater')
  
  # then add to the previous object
  q <- p +
    new_scale_fill() +
    geom_tile(aes(fct_inorder(greater_df$Receiver) , Interaction, fill = `pvalue filtered`), data = greater_df) +
    scale_fill_stepsn("UP in disease", colours = c('#FF0000','#FF9191', '#FFDADA', '#FFFFFF' ), breaks = c(0, 0.001,0.05,0.1), limits = c(0, 0.11), na.value = 'white') +
    facet_grid(rows = vars(factor(Sender))) +
    coord_flip() +
    theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = "bottom")
  return(q)
}

# Wrapper funtion for all three funtions above
full_run_sampleWise <- function(disease_contrast) {
  sender_cells <- c('CD4', 'CD8')
  receiver_cells <- c('CD4', 'CD8', 'CD16', 'CD19')
  
  results_df <- data.frame('Interaction' = character(),
                           'pvalue' = numeric(),
                           'adj.pvalue' = numeric(),
                           'direction' = character(),
                           'Sender' = character(),
                           'Receiver' = character(),
                           'Contrast' = factor())
  
  for (j in sender_cells) {
    for (i in receiver_cells) {
      comp_df <- create_dataframe(sender_contrast = j, receiver_contrast = i, disease_contrast = disease_contrast)
      
      contrast_df <- run_Fishers_contrast(comp_df = comp_df,
                                          sender_contrast = j,
                                          receiver_contrast = i,
                                          disease_contrast = disease_contrast)
      results_df <- results_df %>% bind_rows(contrast_df)
    }
  }
  save(results_df,file = '/home/lisbet/Data/2021_GinoRoche/Data/DF_stat_CellPhoneDB_output.RData')
  p_heatmap <- plot_contrast_high_low(results_df = results_df, which_disease = disease_contrast[1], threshold_p = 0.1)
  return(p_heatmap)
}
          
