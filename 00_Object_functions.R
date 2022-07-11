### 'gino.kwon@drfz.de' or 'gino.kwon@gmail.com' corresponded ###

# qunatile tranformation
func_quantile <- function(x){
  df_rank <- apply(x,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(x, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)  
  index_to_mean <- function(my_index, my_mean){ return(my_mean[my_index]) }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(x)
  return(df_final)
}

# compute raw sum of reads
tmp_func_count <- function(x1,x2,x3,x4,c1,c2,tb=tb_4,lab) {
  y1 <- names(x1)[which( (x1 >= quantile(x1,c1)) & (x1 <= quantile(x1,c2) ) ) ]
  y2 <- names(x2)[which( (x2 >= quantile(x2,c1)) & (x2 <= quantile(x2,c2) ) ) ]
  y3 <- names(x3)[which( (x3 >= quantile(x3,c1)) & (x3 <= quantile(x3,c2) ) ) ]
  y4 <- names(x4)[which( (x4 >= quantile(x4,c1)) & (x4 <= quantile(x4,c2) ) ) ]
  log <- unique( c(y1,y2,y3,y4) )
  pdf(file=paste(lab,'.pdf',sep=''))
  y5 <- tb[which(rownames(tb) %in% log),]
  y6 <- as.numeric(y5)
  plot(density(y6),main=length(log),xlab='Log2 RPKM',ylab='Freq')
  dev.off()
  return(log) 
}

# to compute Treated / Control
func_expression <- function(cnt,tr) {
  cnt <- as.numeric(cnt)
  tr <- as.numeric(tr)
  vals <- (tr - cnt) / cnt
  return(vals)
}

# to import Genesets GMT format to R list
func_imp_path <- function(path,set,min,max) {
  min <- as.numeric(min)
  max <- as.numeric(max)
  raw_path <- readLines(path)
  ls_path <- as.list(rep(NA,length(raw_path)) )
  ls_path_len <- rep(0,length(raw_path))
  ls_path_org <- ls_path
  ls_path_org_len <- ls_path_len
  for(i in 1:length(raw_path)) {
    tmp <- unlist(strsplit(raw_path[i],split='\t',fixed=T) )
    names(ls_path)[i] <- tmp[1]
    names(ls_path_org)[i] <- tmp[1]
    ls_path_org[[i]] <- tmp[-(1:2)]
    ls_path_org_len[i] <- length(tmp)-2
    tmp <- tmp[-(1:2)]
    if( length( which(!tmp %in% set) ) >= 1 ) { tmp <- tmp[ which(tmp %in% set) ] }
    ls_path[[i]] <- tmp
    ls_path_len[i] <- length(tmp)
  }

  # removed more than 300 and 3
  del_ls <- which( (ls_path_len > max) | ( ls_path_len < min) )
  if( length(del_ls) >= 1) {
  ls_path <- ls_path[-del_ls]
  ls_path_len <- ls_path_len[-del_ls]
  ls_path_org <- ls_path_org
  ls_path_org_len <- ls_path_org_len
  }
  return(list(ls_path,ls_path_len,ls_path_org,ls_path_org_len ) )
}

# extracting( and )splitting) expression table to each Disease stage
func_extr <- function(v1,v2,v3,v4,imp00,imp18,index_labs) {
  v1 <- v1 #index
  v2 <- v2 #donor
  v3 <- v3 #celltype
  v4 <- v4 #Hour
  imp00 <- imp00
  imp18 <- imp18
  index_labs <- index_labs
  # column name 'up' is name of sample
  tar_lab <- index_labs$up[which( (index_labs$index %in% v1)&(index_labs$don == v2)&(index_labs$typ == v3)&(index_labs$hrs == v4) )]
  if(v4 == 'H18') {
    x1 <- imp18[,which(colnames(imp18) %in% tar_lab)]
    x2 <- imp18[,grep(paste(v3,'_HD',sep='' ),colnames(imp18)) ]
    } else {
      x1 <- imp00[,which(colnames(imp00) %in% tar_lab)]
      x2 <- imp00[,grep(paste(v3,'_HD',sep='' ),colnames(imp00)) ]
  }
  return(list(x1,x2))
}

# calculating cumilative distance of group
func_ks_avg_inv <- function(x1,sets) {
  x1 <- x1 #input gene list
  y <- unlist(sets)
  lab <- names(sets)
  gs <- rownames(x1)
  dists <- numeric()
  for(i in 1:ncol(x1)) {
    a1 <- x1[,i]
    names(a1) <- gs
    r1 <- sort(a1)
    l1 <- names(r1)
    bak <- 1:length(l1)
    tar <- which(l1 %in% y)
    # x=back,y=target . Oneside
    log1 <- ks.test(x=bak,y=tar,alternative='greater')
    log2 <- ks.test(x=bak,y=tar,alternative='less')
    d_plus <- log1[[1]]
    d_minus <- log2[[1]]
    p_plus <- log1[[2]]
    p_minus <- log2[[2]]
    if( d_plus > d_minus) {
      d_stat <- as.numeric( d_plus )
    } else {
      d_stat <- -1 * as.numeric( d_minus )
    }
    dists <- append(dists,d_stat)
  }
  med_dist <- median(dists)
  mu_dist <- mean(dists)
  sds_dist <- sd(dists)
  return(list(med_dist,mu_dist,sds_dist))
}

# compute p-value for the gene expression similarity
func_coexpr_perm <- function(inp,prob_num,cand_id,y,pops) {
  tmp_rst <- matrix(NA,nrow=nrow(inp),ncol=length(cand_id))
  rownames(tmp_rst) <- pops #N of population
  colnames(tmp_rst) <- y
  # random value for perturbation
  rand_y <- lapply(1:prob_num, function(i) sapply(1:ncol(inp), function(j) { inp[(sample(nrow(inp),1)),j] } ) )
  # Probability = Random < Real expression A-B / Number of permutation
  for(z in 1:ncol(tmp_rst)) {
    real_x <- inp[cand_id[z],]
    real_y <- lapply(1:nrow(tmp_rst), function(i) inp[i,] ) # ordered by Pops. list
    # Calc Dist
    real_dist <- unlist(lapply(1:length(real_y), function(i) sqrt(sum((real_x-real_y[[i]])^2))  )) #Num Vec. Perm to Each gene
    rand_dist_bak <- unlist(lapply(1:length(rand_y), function(i) sqrt(sum((real_x-rand_y[[i]])^2))  ))
    # calc Pvalue per gene
    pvals <- sapply(1:nrow(tmp_rst), function(i) length(which( rand_dist_bak < real_dist[i]))/prob_num )
    tmp_rst[,z] <- pvals
  }
  # If Input geneset has less than 5 significant purterbation. stop.
  if( sum(as.numeric(tmp_rst)) >= 5 ) {
    return(tmp_rst)
  } else {
    return("ERROR")
  }
}

# Filtering out, post-procedures. by PPI. But PPI wasn't applied
func_connect <- function(imp,rms,ppi) {
  # Filtering by Squre-root-mean. (=genes that have no differences between samples)
  #ids <- which( rownames(imp) %in% rms$sym[which(rms$pass == 'T')] )
  #anz <- imp[ids,]
  # Mute RMSE

  # remove bias + non-effect
  sumrow <- sapply(1:nrow(anz), function(i) sum(anz[i,])  )
  sumcol <- sapply(1:ncol(anz), function(i) sum(anz[,i])  )
  # for unique interaction between Gene-Gene. selected genes that 1. At least 3 co-funtional targets(potential candidates). 2. Not too many or co-expressed to all genes(= regarding as Normal signal such as Median+abs(Sigma))
  # Maximum interaction(conntection) threshold is 5% of Gene population.
  ids_col <- which( (sumcol >= 2) & (sumcol <= (nrow(anz)*0.05) ) )
  # Deleting Perturbation(member of geneset)
  if( length(ids_col) >= 3 ) {
    sub_anz <- anz[,ids_col]
    if((nrow(sub_anz) >= 5) && (length(sub_anz) >= 1) && (is.matrix(sub_anz))) {
      sumrow <- sapply(1:nrow(sub_anz), function(i) sum(sub_anz[i,])  )
      ids_row <- which( sumrow > 1)
      sub_anz <- sub_anz[ids_row,]
      if((nrow(sub_anz) >= 5) && (length(sub_anz) >= 1) && (is.matrix(sub_anz))) {
        sumcol <- sapply(1:ncol(sub_anz), function(i) sum(sub_anz[,i])  )
        ids_col <- which( sumcol > 1)
        sub_anz <- sub_anz[,ids_col]
        # PPI filtering (PPI >= 0.6), didn't used. Source : STRING DB
        if((nrow(sub_anz) >= 5) && (length(sub_anz) >= 1) && (is.matrix(sub_anz))) {
          samp <- rownames(sub_anz)
          ppi$p1 <- sapply(1:nrow(ppi), function(i) if(ppi$a[i] %in% samp) {1} else {0} )
          ppi$p2 <- sapply(1:nrow(ppi), function(i) if(ppi$b[i] %in% samp) {1} else {0} )
          ids <- (ppi$p1 + ppi$p2)
          ids <- which(ids == 2)
          syms <- unique(append(ppi$a[ids],ppi$b[ids]) )
          samp_ppi <- samp[which(samp %in% syms)]
          ppi_anz <- sub_anz[which(samp %in% syms),]
          return(list(sub_anz,ppi_anz))
        } else {
          return(list(matrix(NA),matrix(NA)))
        }
      } else {
        return(list(matrix(NA),matrix(NA)))
      }
    } else {
      return(list(matrix(NA),matrix(NA)))
    }

  } else {
    return(list(matrix(NA),matrix(NA)))
  }
}

# Geneset enrichment to other Genesets. Hypergeometric test
func_ann <- function(anz,sets,outfile) {
  pop <- 10835 #size of gene population
  tmp <- paste(unlist(strsplit(outfile,split='_',fixed=T))[-(1:2)],collapse='_')
  # Batch_number(=any), Origin (=Input Geneset), Option(=Dummy), Label(=Name of Geneset), DB(=Original source of DB), Count(=N Enriche gene), Size(=Original Geneset size)
  hds <- c('Batch','Origin','Option','Label','DB','Count','P_value','Genes','Size')
  form <- matrix(NA,nrow=length(sets),ncol=length(hds))
  colnames(form) <- hds
  form <- as.data.frame(form,stringsAsFactors=F)
  form$Batch <- as.character( strsplit(outfile,split='_',fixed=T)[[1]][1] )
  form$Origin <- strsplit(outfile,split='_',fixed=T)[[1]][2]
  form$Option <- tmp
  samp <- rownames(anz)
  for(i in 1:length(sets)) {
    succ <- sets[[i]]
    label <- names(sets)[i]
    DB <- strsplit(label,split='_',fixed=T)[[1]][1]
    if(! DB %in% c('BIOCARTA','HALLMARK','KEGG','PID','REACTOME','MB','REF') ) { DB <- 'PUB' }
    samp_succ <- samp[which(samp %in% succ)]
    # Minimum Enriched gene is 3. Smaller intersection than 3 regarded as P-value as 1.0
    if( length(samp_succ) >= 3) {
      p <- phyper((length(samp_succ)-1),length(succ),pop-length(succ),length(samp), lower.tail=FALSE)
      form$DB[i] <- DB
      form$Count[i] <- length(samp_succ)
      form$Genes[i] <- paste(samp_succ,collapse='|')
      form$Size[i] <- length(succ)
    } else {
      p <- 1
      samp_succ <- NA
    }
  form$P_value[i] <- p
  form$Label[i] <- label
  }
  # False-discovery rate
  form$FDR <- p.adjust(form$P_value,method='fdr')
  #filtering
  filter <- form[which( (form$FDR <= 0.10) & (form$DB %in% c('BIOCARTA','KEGG','REACTOME','PID','REF')) ),]
  return(list(filter,form))
}

# Main part of Co-functional candidates
func_coexpr_main <- function(gmt,cand,ind,typ='all',don='all',stg='all',hrs='all',decimal=5) {
  pops <- rownames(gmt) #number of genes
  x <- gmt #matrix format of numeric expression
  y <- cand #Character for mat; list of genes
  indtb <- ind #Disease Index or Stage ifnormation table. must have don,typ,stg (donor,celltype,stage)
  if( (typ == 'all') && (length(typ)==1) ) { typ <- c('CD4','CD8','CD16','CD19') }
  if( (don == 'all') && (length(don)==1) ) { don <- c('HD','pSS','SLE') }
  if( (stg == 'all') && (length(stg)==1) ) { stg <- c('Low','Mid','High') }
  prob_num <- 10^decimal #number of trier, decimal of p-value
  cand_id <- unlist(lapply(1:length(y), function(i) which(y[i] == pops) )) #getting input from genesets
  # Minimum number of gene is 5
  if( length(cand_id) >= 5 ) {
    # Collecting Stage information from Input table(Index). indtb
    ids <- which( (indtb$typ %in% typ) & (indtb$don %in% don) & (indtb$stg %in% stg) )
    x1 <- x[,ids]
    samp_inf <- indtb[ids,]
    # Minimum number of sample is 2. otherwise can't computeZ
    if( nrow(samp_inf) >= 2 ) {
      # get random permutatin exoression
      tb_pvals <- func_coexpr_perm(x1,prob_num,cand_id,y,pops) #compute P-value per each input genes
      if( (tb_pvals != "ERROR") && (is.matrix(tb_pvals)) ) {
        tb_fdrs <- tb_pvals
        # FDR
        for(i in 1:ncol(tb_pvals)) { tb_fdrs[,i] <- p.adjust(tb_pvals[,i],method='fdr') }   
        # remove self loop and initial nodes
        tb_fdrs <- tb_fdrs[-cand_id,]
        # Make a Binary matrix. 1 is interacted, 0 is not. Row is input, Col is target
        tb_bin <- as.matrix( tb_fdrs )
        tb_bin[which(tb_bin >= 0)] <- 0
        tb_bin[which(tb_fdrs < 0.1)] <- 1
        # FDR 10% as threshold
        tb_bin_take <- tb_bin
        # remaining only Significant connection. #1 All(without filtering) #2 PPI filter. PPI not used.
        data_mat <- func_connect(tb_bin_take,tb_rmse,ppi)
        anz <- data_mat[[1]]
        anz_ppi <- data_mat[[2]]
        # Co-funtional candidate genes TO Geneset enrichement P-value
        if( (length(anz) >= 1) && (nrow(anz) >= 3 ) ) {
          annotat <- func_ann(anz,path_set,outfile) #matching Geneset
          ann_filter <- annotat[[1]]
          write.table(ann_filter,file=paste(outfile,'.txt',sep=''),sep='\t',col.names=T)
          } else {
            cat('Anotation_FAIL','\n')
        }
      }
    }
  }
}


