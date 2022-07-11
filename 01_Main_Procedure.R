### 'gino.kwon@drfz.de' or 'gino.kwon@gmail.com'  ###

library(plyr)
library(EnsDb.Hsapiens.v86)

source("00_Object_functions.R")
##### Normalization and Builing analyze set #####
# import FPKM (=RPKM)
path_rpkm <- '/home/kwon/project/02_sle/data/gct_inf/rpkm_01.gct'#repository name 'Immunomics_PD1_ago_merged_1_2_3_rpkms_collapsed.gct'
raw_tb <- read.delim(file=path_rpkm,skip=2,stringsAsFactors=F,header=T,check.names=F,sep='\t')
tb_samp <- colnames(raw_tb)[-(1:2)] #sample list
tb_samp <- as.data.frame(tb_samp,stringsAsFactors=F)
#write.table(tb_samp,file='/home/kwon/project/02_sle/data/gct_inf/table_inf/tb_samp_00.txt',sep='\t',col.names=F,row.names=F)

# Load biotypes
tb_entrez <- raw_tb[,c(1:2)]
colnames(tb_entrez) <- c('entrez','symbol')
input_keys <- as.character( tb_entrez$entrez )
tb_biotype <- ensembldb::select(EnsDb.Hsapiens.v86, keys=input_keys,columns=c("SYMBOL", "GENEBIOTYPE"), keytype="ENTREZID")
tb_entrez$type <- NA
for(i in 1:nrow(tb_entrez)) {
	if( tb_entrez[i,1] %in% tb_biotype[,1] ) {
		log <- tb_biotype[ which( tb_entrez[i,1] == tb_biotype[,1]), ]
		if( nrow(log) >= 2) {
			if( "protein_coding" %in% log[,3] ) {
				log1 <- "protein_coding"
			} else {
				log1 <- log[1,3]
			}
		} else {
			log1 <- log[1,3]
		}
		tb_entrez$type[i] <- log1
	}
}

# prioroty genes for the data processing
base_genes <- c("STAT1","STAT2","STAT3","STAT4","BTLA","PDCD1","CD274","PDCD1LG2","HAVCR2","CD22","CD40","CD40LG",
"SYK","PLCG2","BTK","PIK3CA","AKT1","IRF4","IRF5")

# to select biotype
entrz_protein <- tb_entrez[which(tb_entrez$type == 'protein_coding'),1]
# remove non coding
raw_tb_0 <- raw_tb[which(raw_tb[,1] %in% entrz_protein),]

# import 1st re-named sample list
path_samp_lab1 <- '/home/kwon/project/02_sle/data/gct_inf/table_inf/samp_01_rename.txt'
raw_samp_lab <- read.delim(file=path_samp_lab1,sep='\t',header=T,stringsAsFactors=F)
# import sorted order of 1st re-named list
path_samp_sort <- '/home/kwon/project/02_sle/data/gct_inf/table_inf/samp_02_sort.txt'
raw_samp_sort <- read.delim(file=path_samp_sort,sep='\t',header=T,stringsAsFactors=F)
# change col names to re-named headers
colnames(raw_tb_0) <- c(colnames(raw_tb_0)[1:2],unlist(lapply(1:ncol(raw_tb_0), function(x) raw_samp_lab$New_Label[which(raw_samp_lab$Raw_Header==colnames(raw_tb_0)[x] )]  )))

# remove IgG samples
rm_ls <- raw_samp_lab$New_Label[ which(raw_samp_lab$Pass==F) ] #col_id of IgG
raw_tb_1 <- raw_tb_0[, -( which(colnames(raw_tb_0) %in% rm_ls) ) ]
tb_lab <- raw_tb_1[,(1:2)] #Table row Entrez and symbol
raw_tb_2 <- raw_tb_1[,-(1:2)]

# log2 scale
raw_tb_3 <- apply(raw_tb_2, 1:2, function(x) log((x+1),2) )
raw_log <- as.numeric(raw_tb_3)
rownames(raw_tb_3) <- as.character(raw_tb_row_lab[,1])

# remove Non expressed genes (non-mapped) and overexpression
tmp_sum <- unlist( lapply(1:nrow(raw_tb_3), function(x)	sum( raw_tb_3[x,]	)	) ) #Sum of value per gene
names(tmp_sum) <- tb_lab[,1]
tmp_sum <- sort(tmp_sum) #sorting sum of value
tmp_tot_all <- tmp_sum # row_sum_tot as total read
# to collect name of samples per donor on stimulated condition, and compute
ls_sle <- raw_samp_sort$New_Label[ which( (raw_samp_sort$V2_Class == 'SLE') & (raw_samp_sort$V3_Condition == 'aCD3') ) ]
ls_pss <- raw_samp_sort$New_Label[ which( (raw_samp_sort$V2_Class == 'pSS') & (raw_samp_sort$V3_Condition == 'aCD3') ) ]
ls_hd <- raw_samp_sort$New_Label[ which( (raw_samp_sort$V2_Class == 'HD') & (raw_samp_sort$V3_Condition == 'aCD3') ) ]
tmp_tot_sle <- unlist( lapply(1:nrow(raw_tb_3), function(x)	sum( raw_tb_3[x,which( colnames(raw_tb_3) %in% ls_sle ) ]	)	) )
tmp_tot_pss <- unlist( lapply(1:nrow(raw_tb_3), function(x)	sum( raw_tb_3[x,which( colnames(raw_tb_3) %in% ls_pss ) ]	)	) )
tmp_tot_hd <- unlist( lapply(1:nrow(raw_tb_3), function(x)	sum( raw_tb_3[x,which( colnames(raw_tb_3) %in% ls_hd ) ]	)	) )
names(tmp_tot_sle) <- tb_lab[,1]
names(tmp_tot_pss) <- tb_lab[,1]
names(tmp_tot_hd) <- tb_lab[,1]
# to collect gene symbols within range
case0352 <- tmp_func_count(x1=tmp_tot_sle,x2=tmp_tot_pss,x3=tmp_tot_hd,c1=0.45,c2=0.995,x4=tmp_tot_all,tb=raw_tb_3,lab='Case_0352_45Percto99Percent')
# final decision. Nonmapped+Low reads regions as 45%, Overexpressed as 99.5%, considering after the "base_genes"
raw_tb_4 <- raw_tb_3[which(rownames(raw_tb_3) %in% case0352),]

# to swap log2rpkm to quantile normalization
raw_tb_5 <- func_quantile(raw_tb_4)

# to switch Entrez to Symbol
rownames(raw_tb_5) <- sapply(1:nrow(raw_tb_5), function(x) raw_tb_row_lab$Description[ which( as.character( raw_tb_row_lab$Name ) == rownames(raw_tb_5)[x]) ]   )
bak <- raw_tb_5

# import sorted order of samples, manually sorted by Treated condition, Celltype, Donor
path_samp_ord <- '/home/kwon/project/02_sle/rimage/Sample_order.txt'
tb_samp_ord <- read.delim(file=path_samp_ord,sep='\t',stringsAsFactors=F,header=T)

# to split Steady-state matrix (Steadystate == H00, anz_00)
mat_00 <- raw_tb_5[,grep('BASE',colnames(raw_tb_5))]
# Labelling to Final Sample name and Sorting
ids <- sapply(1:ncol(mat_00), function(x) which(colnames(mat_00) == tb_samp_ord$old[x] ) )
anz_00 <- mat_00[,ids]
colnames(anz_00) <- tb_samp_ord$new[which(tb_samp_ord$hour == 'H00')]
tmp <- anz_00
for(i in 1:ncol(tmp)) {
    tmp[,i] <- rank(imp[,i])
}
rnk_00 <- tmp #rank table for Steady-state
# upper limit same as up-rpkm, than z-scaling and zero centralizing
tmp <- anz_00
tmp[which(tmp > 8)] <- 8
fit_00 <- tmp
for(i in 1:ncol(tmp)) {
        ms <- mean(tmp[,i] )
        sds <- sd(tmp[,i])
        ys <- ( (tmp[,i] - ms) / sds )
        fit_00[,i] <- ys
}

# final normalized analyze set
norm_00 <- anz_00 #unfitted
rank_00 <- rnk_00 #rank order
anz_00 <- fit_00 #normalized set

# Making up stimulated set
# three main variables
v1 <- c("HD","pSS","SLE")
v2 <- c("PBS","aCD3")
v3 <- c("CD8","CD4","CD19","CD16")
# to extract 18hours table
mat_18 <- raw_tb_5[,-(grep('BASE',colnames(raw_tb_5)))]
# remove not pair sample( absent of either PBS or aCD3) / Variable "Singular" 99 is Absent ID number of sample 
mat_18_pair <- mat_18[,which(!colnames(mat_18) %in% samp_sort_num$New_Label[which(samp_sort_num$Singluar == 99)] )]
# to identify column number of each PBS(as 0) and aCD3
tmp_id <- 1:ncol(mat_18_pair)
tmp_hds <- unlist(lapply(1:ncol(mat_18_pair), function(x) if( strsplit( colnames(mat_18_pair)[x],split='_',fixed=T ) [[1]][3] == "aCD3" ) {x} else {0} ))
tmp_tar <- sapply(1:length(tmp_hds), function(x) 
if( tmp_hds[x] == 0 ) { which( colnames(mat_18_pair) == gsub("_PBS_","_aCD3_",colnames(mat_18_pair)[x] ) ) } else {NA}
)
# make result table template / Row = Genes, Col = Binded samples (half of total) 
rst <- matrix(NA,nrow=nrow(mat_18_pair),ncol=length(which(tmp_hds != 0)) )
tb_lab <- matrix(NA,nrow=ncol(rst),ncol=7) #sample description for sorting
colnames(tb_lab) <- c('label','typ','don','id','loc','ord','sc')
tb_lab <- as.data.frame(tb_lab,stringsAsFactors=F)
# to get columnd id of PBS(control) and TR(aCD3). based on Sample ID matching.
col_id_cnt <- which(!is.na(tmp_tar))
col_id_tr <- unlist(lapply(1:length(col_id_cnt), function(x) tmp_tar[col_id_cnt[x]]     ))
rst_keep <- as.data.frame(matrix(NA,nrow=nrow(rst),ncol=1),stringsAsFactors=F)
# when Control value has non-expression (under 2% within the sample, it will changed to Zero. If Control is zero expression, taking expression of Tretaed(aCD3) 
percent_cut <- 0.05
limit_fc <- 2
for(i in 1:ncol(rst)) {
	vals_cnt <- mat_18_pair[,col_id_cnt[i]]
	vals_cnt[which( vals_cnt < as.numeric( quantile(vals_cnt,percent_cut) ) )] <- 0
	vals_tr <- mat_18_pair[,col_id_tr[i]]
	# compute expression value
	vals_fc <- func_expression(cnt=vals_cnt,tr=vals_tr)
	# swap value erro
	vals_fc[which(is.na(vals_fc) ) ] <- 0
	vals_fc[which( vals_fc > limit_fc) ] <- limit_fc
	rst[,i] <- vals_fc
	tmp <- as.data.frame(vals_cnt,strinsAsFactors=F)
	colnames(tmp) <- colnames(mat_18_pair)[col_id_cnt[i]]
	rst_keep <- cbind(rst_keep,tmp)
	tmp <- as.data.frame(vals_tr,strinsAsFactors=F)
	colnames(tmp) <- colnames(mat_18_pair)[col_id_tr[i]]
	rst_keep <- cbind(rst_keep,tmp)
	# Define name of samples, for now. but changing to the final name later
	raw_lab <- colnames(mat_18_pair)[col_id_cnt[i]]
	raw_del <- unlist( strsplit(raw_lab,split='_',fixed=T) )
	don <- raw_del[2]
	typ <- raw_del[4]
	ids <- raw_del[5]
	# Later Batch number experiment, add ID 50 number, during the NGS handled batch effect
	if( raw_del[1] == "01" ) {
		ids <- as.character( as.numeric(ids) + 50 )
	}
	new_lab <- paste(typ,don,ids,sep='_')
	tb_lab$label[i] <- new_lab
	tb_lab$typ[i] <- typ
	tb_lab$don[i] <- don
	tb_lab$id[i] <- ids
	tb_lab$loc[i] <- i
	sc <- 0
	# value for column sorting
	if( typ == "CD4" ) { sc <- sc+10 } else if( typ == "CD8" ) { sc <- sc+20 } else if( typ == "CD16" ) { sc <- sc+30 } else { sc <- sc+40}
	if( don == "HD" ) { sc <- sc+1 } else if( don == "pSS") { sc <- sc+2} else { sc <- sc+3 }
	tb_lab$sc[i] <- sc 
}
rownames(rst) <- rownames(mat_18_pair)
colnames(rst) <- tb_lab$label
rst_keep <- rst_keep[,-1] # remove null column
# makr order and arranging
tb_lab <- arrange(tb_lab,sc,id)
anz_18 <- rst
for( i in 1:nrow(tb_lab)) {
	anz_18[,i] <- rst[,tb_lab$loc[i]]
}
colnames(anz_18) <- tb_lab$label
# rank table
tmp <- anz_18
for(i in 1:ncol(tmp)) {
    tmp[,i] <- rank(anz_18[,i])
}
rnk_18 <- tmp #rank table for Steady-state
# scaling and maximum value fitting as absolute 1
fit_18 <- anz_18
fit_18[which(fit_18 > 1)] <- 1
fit_18[which(fit_18 < -1)] <- -1

# final normalized analyze set
norm_18 <- anz_18 #unfitted
rank_18 <- rnk_18 #rank order
anz_18 <- fit_18 #normalized set

##### EOF Normalization and Builing analyze set #####





##### Regulatory (Enrichment level) of samples #####
# import Geneset file
genes <- rownames(anz_18) # symbol list of population
inpathway='/home/kwon/project/global_data/pathway/c2_v71.gmt'
raw_c2 <- func_imp_path(inpathway,genes,min=3,max=500)
path_set_c2 <- raw_c2[[1]] #LIST
path_set <- path_set_c2[which( sapply(1:length(path_set_c2), function(x) strsplit(names(path_set_c2)[x],split='_',fixed=T)[[1]][1] ) %in% c('KEGG','REACTOME','BIOCARTA','PID')	) ]
tmp <- c('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING','REACTOME_INTERFERON_GAMMA_SIGNALING')
path_ifn <- path_set[unlist(lapply(1:length(tmp), function(x) which(tmp[x] == names(path_set) )   ))]
# import Disease index information per patient samples
index_path <- '/home/kwon/project/02_sle/data/clin_cytof_luisa/index_label_table_ext.txt'
index_labs <- read.delim(file=index_path,sep='\t',stringsAsFactors=F,header=T)

# Extract(split) matrix
# H18(aCD3) / SLE : 2 / 4+5
tmp <- func_extr( c(2), "SLEDAI", "CD4", "H18", anz_00,anz_18,index_labs )
mat_H18_CD4_SLE_low <- tmp[[1]]
mat_H18_CD4_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD4", "H18", anz_00,anz_18,index_labs )
mat_H18_CD4_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD8", "H18", anz_00,anz_18,index_labs )
mat_H18_CD8_SLE_low <- tmp[[1]]
mat_H18_CD8_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD8", "H18", anz_00,anz_18,index_labs )
mat_H18_CD8_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD16", "H18", anz_00,anz_18,index_labs )
mat_H18_CD16_SLE_low <- tmp[[1]]
mat_H18_CD16_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD16", "H18", anz_00,anz_18,index_labs )
mat_H18_CD16_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD19", "H18", anz_00,anz_18,index_labs )
mat_H18_CD19_SLE_low <- tmp[[1]]
mat_H18_CD19_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD19", "H18", anz_00,anz_18,index_labs )
mat_H18_CD19_SLE_high <- tmp[[1]]
# H18(aCD3) / pSS : 0/1+2+3 / 4+10
tmp <- func_extr( c(0), "ESSDAI", "CD4", "H18", anz_00,anz_18,index_labs )
mat_H18_CD4_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD4", "H18", anz_00,anz_18,index_labs )
mat_H18_CD4_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,10), "ESSDAI", "CD4", "H18", anz_00,anz_18,index_labs )
mat_H18_CD4_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD8", "H18", anz_00,anz_18,index_labs )
mat_H18_CD8_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD8", "H18", anz_00,anz_18,index_labs )
mat_H18_CD8_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,10), "ESSDAI", "CD8", "H18", anz_00,anz_18,index_labs )
mat_H18_CD8_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD16", "H18", anz_00,anz_18,index_labs )
mat_H18_CD16_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD16", "H18", anz_00,anz_18,index_labs )
mat_H18_CD16_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,10), "ESSDAI", "CD16", "H18", anz_00,anz_18,index_labs )
mat_H18_CD16_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD19", "H18", anz_00,anz_18,index_labs )
mat_H18_CD19_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD19", "H18", anz_00,anz_18,index_labs )
mat_H18_CD19_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,10), "ESSDAI", "CD19", "H18", anz_00,anz_18,index_labs )
mat_H18_CD19_pSS_high <- tmp[[1]]
# H00(Stadystate) / SLE : 2 / 4+5 / 8+10+12
tmp <- func_extr( c(2), "SLEDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_SLE_low <- tmp[[1]]
mat_H00_CD4_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_SLE_mid <- tmp[[1]]
tmp <- func_extr( c(8,10,12), "SLEDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_SLE_low <- tmp[[1]]
mat_H00_CD8_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_SLE_mid <- tmp[[1]]
tmp <- func_extr( c(8,10,12), "SLEDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_SLE_low <- tmp[[1]]
mat_H00_CD16_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_SLE_mid <- tmp[[1]]
tmp <- func_extr( c(8,10,12), "SLEDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_SLE_high <- tmp[[1]]
tmp <- func_extr( c(2), "SLEDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_SLE_low <- tmp[[1]]
mat_H00_CD19_HD <- tmp[[2]]
tmp <- func_extr( c(4,5), "SLEDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_SLE_mid <- tmp[[1]]
tmp <- func_extr( c(8,10,12), "SLEDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_SLE_high <- tmp[[1]]
# H00(SteadyState) / pSS : 0 / 1+2+3 / 4+5+6+7 / 9+10
tmp <- func_extr( c(0), "ESSDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,5,6,7), "ESSDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_pSS_mid <- tmp[[1]]
tmp <- func_extr( c(9,10), "ESSDAI", "CD4", "H00", anz_00,anz_18,index_labs )
mat_H00_CD4_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,5,6,7), "ESSDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_pSS_mid <- tmp[[1]]
tmp <- func_extr( c(9,10), "ESSDAI", "CD8", "H00", anz_00,anz_18,index_labs )
mat_H00_CD8_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,5,6,7), "ESSDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_pSS_mid <- tmp[[1]]
tmp <- func_extr( c(9,10), "ESSDAI", "CD16", "H00", anz_00,anz_18,index_labs )
mat_H00_CD16_pSS_high <- tmp[[1]]
tmp <- func_extr( c(0), "ESSDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_pSS_non <- tmp[[1]]
tmp <- func_extr( c(1,2,3), "ESSDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_pSS_low <- tmp[[1]]
tmp <- func_extr( c(4,5,6,7), "ESSDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_pSS_mid <- tmp[[1]]
tmp <- func_extr( c(9,10), "ESSDAI", "CD19", "H00", anz_00,anz_18,index_labs )
mat_H00_CD19_pSS_high <- tmp[[1]]
# All stage (= binding all index)
mat_H18_CD4_SLE_all <- cbind(mat_H18_CD4_SLE_low,mat_H18_CD4_SLE_high)
mat_H18_CD8_SLE_all <- cbind(mat_H18_CD8_SLE_low,mat_H18_CD8_SLE_high)
mat_H18_CD16_SLE_all <- cbind(mat_H18_CD16_SLE_low,mat_H18_CD16_SLE_high)
mat_H18_CD19_SLE_all <- cbind(mat_H18_CD19_SLE_low,mat_H18_CD19_SLE_high)
mat_H18_CD4_pSS_all <- cbind(mat_H18_CD4_pSS_low,mat_H18_CD4_pSS_high)
mat_H18_CD8_pSS_all <- cbind(mat_H18_CD8_pSS_low,mat_H18_CD8_pSS_high)
mat_H18_CD16_pSS_all <- cbind(mat_H18_CD16_pSS_low,mat_H18_CD16_pSS_high)
mat_H18_CD19_pSS_all <- cbind(mat_H18_CD19_pSS_low,mat_H18_CD19_pSS_high)
mat_H00_CD4_SLE_all <- cbind(mat_H00_CD4_SLE_low,mat_H00_CD4_SLE_mid,mat_H00_CD4_SLE_high)
mat_H00_CD8_SLE_all <- cbind(mat_H00_CD8_SLE_low,mat_H00_CD8_SLE_mid,mat_H00_CD8_SLE_high)
mat_H00_CD16_SLE_all <- cbind(mat_H00_CD16_SLE_low,mat_H00_CD16_SLE_mid,mat_H00_CD16_SLE_high)
mat_H00_CD19_SLE_all <- cbind(mat_H00_CD19_SLE_low,mat_H00_CD19_SLE_mid,mat_H00_CD19_SLE_high)
mat_H00_CD4_pSS_all <- cbind(mat_H00_CD4_pSS_low,mat_H00_CD4_pSS_mid,mat_H00_CD4_pSS_high)
mat_H00_CD8_pSS_all <- cbind(mat_H00_CD8_pSS_low,mat_H00_CD8_pSS_mid,mat_H00_CD8_pSS_high)
mat_H00_CD16_pSS_all <- cbind(mat_H00_CD16_pSS_low,mat_H00_CD16_pSS_mid,mat_H00_CD16_pSS_high)
mat_H00_CD19_pSS_all <- cbind(mat_H00_CD19_pSS_low,mat_H00_CD19_pSS_mid,mat_H00_CD19_pSS_high)

# calculate distance(=cumulative dist =Enrichment level)
typs <- c('CD4','CD8','CD16','CD19')
dons <- c('pSS','SLE')
path_cand <- path_set #input genesets (candidates)
outlog <- '/home/kwon/project/02_sle/work/64_add_reg/sds.txt' #writting variant

# loop for all genesets
for(zzz in 1:length(path_cand)) {
	print(zzz)
	setname <- names(path_cand[zzz])
	ksdist_H18_CD4_HD <- func_ks_avg_inv(mat_H18_CD4_HD,path_cand[zzz])
	ksdist_H18_CD4_pSS_Low <- func_ks_avg_inv(mat_H18_CD4_pSS_low,path_cand[zzz])
	ksdist_H18_CD4_pSS_High <- func_ks_avg_inv(mat_H18_CD4_pSS_high,path_cand[zzz])
	ksdist_H18_CD8_HD <- func_ks_avg_inv(mat_H18_CD8_HD,path_cand[zzz])
	ksdist_H18_CD8_pSS_Low <- func_ks_avg_inv(mat_H18_CD8_pSS_low,path_cand[zzz])
	ksdist_H18_CD8_pSS_High <- func_ks_avg_inv(mat_H18_CD8_pSS_high,path_cand[zzz])
	ksdist_H18_CD16_HD <- func_ks_avg_inv(mat_H18_CD16_HD,path_cand[zzz])
	ksdist_H18_CD16_pSS_Low <- func_ks_avg_inv(mat_H18_CD16_pSS_low,path_cand[zzz])
	ksdist_H18_CD16_pSS_High <- func_ks_avg_inv(mat_H18_CD16_pSS_high,path_cand[zzz])
	ksdist_H18_CD19_HD <- func_ks_avg_inv(mat_H18_CD19_HD,path_cand[zzz])
	ksdist_H18_CD19_pSS_Low <- func_ks_avg_inv(mat_H18_CD19_pSS_low,path_cand[zzz])
	ksdist_H18_CD19_pSS_High <- func_ks_avg_inv(mat_H18_CD19_pSS_high,path_cand[zzz])
	ksdist_H18_CD4_HD <- func_ks_avg_inv(mat_H18_CD4_HD,path_cand[zzz])
	ksdist_H18_CD4_SLE_Low <- func_ks_avg_inv(mat_H18_CD4_SLE_low,path_cand[zzz])
	ksdist_H18_CD4_SLE_High <- func_ks_avg_inv(mat_H18_CD4_SLE_high,path_cand[zzz])
	ksdist_H18_CD8_HD <- func_ks_avg_inv(mat_H18_CD8_HD,path_cand[zzz])
	ksdist_H18_CD8_SLE_Low <- func_ks_avg_inv(mat_H18_CD8_SLE_low,path_cand[zzz])
	ksdist_H18_CD8_SLE_High <- func_ks_avg_inv(mat_H18_CD8_SLE_high,path_cand[zzz])
	ksdist_H18_CD16_HD <- func_ks_avg_inv(mat_H18_CD16_HD,path_cand[zzz])
	ksdist_H18_CD16_SLE_Low <- func_ks_avg_inv(mat_H18_CD16_SLE_low,path_cand[zzz])
	ksdist_H18_CD16_SLE_High <- func_ks_avg_inv(mat_H18_CD16_SLE_high,path_cand[zzz])
	ksdist_H18_CD19_HD <- func_ks_avg_inv(mat_H18_CD19_HD,path_cand[zzz])
	ksdist_H18_CD19_SLE_Low <- func_ks_avg_inv(mat_H18_CD19_SLE_low,path_cand[zzz])
	ksdist_H18_CD19_SLE_High <- func_ks_avg_inv(mat_H18_CD19_SLE_high,path_cand[zzz])
	ksdist_H00_CD4_HD <- func_ks_avg_inv(mat_H00_CD4_HD,path_cand[zzz])
	ksdist_H00_CD4_pSS_Low <- func_ks_avg_inv(mat_H00_CD4_pSS_low,path_cand[zzz])
	ksdist_H00_CD4_pSS_Mid <- func_ks_avg_inv(mat_H00_CD4_pSS_mid,path_cand[zzz])
	ksdist_H00_CD4_pSS_High <- func_ks_avg_inv(mat_H00_CD4_pSS_high,path_cand[zzz])
	ksdist_H00_CD8_HD <- func_ks_avg_inv(mat_H00_CD8_HD,path_cand[zzz])
	ksdist_H00_CD8_pSS_Low <- func_ks_avg_inv(mat_H00_CD8_pSS_low,path_cand[zzz])
	ksdist_H00_CD8_pSS_Mid <- func_ks_avg_inv(mat_H00_CD8_pSS_mid,path_cand[zzz])
	ksdist_H00_CD8_pSS_High <- func_ks_avg_inv(mat_H00_CD8_pSS_high,path_cand[zzz])
	ksdist_H00_CD16_HD <- func_ks_avg_inv(mat_H00_CD16_HD,path_cand[zzz])
	ksdist_H00_CD16_pSS_Low <- func_ks_avg_inv(mat_H00_CD16_pSS_low,path_cand[zzz])
	ksdist_H00_CD16_pSS_Mid <- func_ks_avg_inv(mat_H00_CD16_pSS_mid,path_cand[zzz])
	ksdist_H00_CD16_pSS_High <- func_ks_avg_inv(mat_H00_CD16_pSS_high,path_cand[zzz])
	ksdist_H00_CD19_HD <- func_ks_avg_inv(mat_H00_CD19_HD,path_cand[zzz])
	ksdist_H00_CD19_pSS_Low <- func_ks_avg_inv(mat_H00_CD19_pSS_low,path_cand[zzz])
	ksdist_H00_CD19_pSS_Mid <- func_ks_avg_inv(mat_H00_CD19_pSS_mid,path_cand[zzz])
	ksdist_H00_CD19_pSS_High <- func_ks_avg_inv(mat_H00_CD19_pSS_high,path_cand[zzz])
	ksdist_H00_CD4_HD <- func_ks_avg_inv(mat_H00_CD4_HD,path_cand[zzz])
	ksdist_H00_CD4_SLE_Low <- func_ks_avg_inv(mat_H00_CD4_SLE_low,path_cand[zzz])
	ksdist_H00_CD4_SLE_Mid <- func_ks_avg_inv(mat_H00_CD4_SLE_mid,path_cand[zzz])
	ksdist_H00_CD4_SLE_High <- func_ks_avg_inv(mat_H00_CD4_SLE_high,path_cand[zzz])
	ksdist_H00_CD8_HD <- func_ks_avg_inv(mat_H00_CD8_HD,path_cand[zzz])
	ksdist_H00_CD8_SLE_Low <- func_ks_avg_inv(mat_H00_CD8_SLE_low,path_cand[zzz])
	ksdist_H00_CD8_SLE_Mid <- func_ks_avg_inv(mat_H00_CD8_SLE_mid,path_cand[zzz])
	ksdist_H00_CD8_SLE_High <- func_ks_avg_inv(mat_H00_CD8_SLE_high,path_cand[zzz])
	ksdist_H00_CD16_HD <- func_ks_avg_inv(mat_H00_CD16_HD,path_cand[zzz])
	ksdist_H00_CD16_SLE_Low <- func_ks_avg_inv(mat_H00_CD16_SLE_low,path_cand[zzz])
	ksdist_H00_CD16_SLE_Mid <- func_ks_avg_inv(mat_H00_CD16_SLE_mid,path_cand[zzz])
	ksdist_H00_CD16_SLE_High <- func_ks_avg_inv(mat_H00_CD16_SLE_high,path_cand[zzz])
	ksdist_H00_CD19_HD <- func_ks_avg_inv(mat_H00_CD19_HD,path_cand[zzz])
	ksdist_H00_CD19_SLE_Low <- func_ks_avg_inv(mat_H00_CD19_SLE_low,path_cand[zzz])
	ksdist_H00_CD19_SLE_Mid <- func_ks_avg_inv(mat_H00_CD19_SLE_mid,path_cand[zzz])
	ksdist_H00_CD19_SLE_High <- func_ks_avg_inv(mat_H00_CD19_SLE_high,path_cand[zzz])
	alldist_H00_CD4_SLE <- func_ks_avg_inv(mat_H00_CD4_SLE_all,path_cand[zzz])
	alldist_H00_CD8_SLE <- func_ks_avg_inv(mat_H00_CD8_SLE_all,path_cand[zzz])
	alldist_H00_CD16_SLE <- func_ks_avg_inv(mat_H00_CD16_SLE_all,path_cand[zzz])
	alldist_H00_CD19_SLE <- func_ks_avg_inv(mat_H00_CD19_SLE_all,path_cand[zzz])
	alldist_H00_CD4_pSS <- func_ks_avg_inv(mat_H00_CD4_pSS_all,path_cand[zzz])
	alldist_H00_CD8_pSS <- func_ks_avg_inv(mat_H00_CD8_pSS_all,path_cand[zzz])
	alldist_H00_CD16_pSS <- func_ks_avg_inv(mat_H00_CD16_pSS_all,path_cand[zzz])
	alldist_H00_CD19_pSS <- func_ks_avg_inv(mat_H00_CD19_pSS_all,path_cand[zzz])
	alldist_H18_CD4_SLE <- func_ks_avg_inv(mat_H18_CD4_SLE_all,path_cand[zzz])
	alldist_H18_CD8_SLE <- func_ks_avg_inv(mat_H18_CD8_SLE_all,path_cand[zzz])
	alldist_H18_CD16_SLE <- func_ks_avg_inv(mat_H18_CD16_SLE_all,path_cand[zzz])
	alldist_H18_CD19_SLE <- func_ks_avg_inv(mat_H18_CD19_SLE_all,path_cand[zzz])
	alldist_H18_CD4_pSS <- func_ks_avg_inv(mat_H18_CD4_pSS_all,path_cand[zzz])
	alldist_H18_CD8_pSS <- func_ks_avg_inv(mat_H18_CD8_pSS_all,path_cand[zzz])
	alldist_H18_CD16_pSS <- func_ks_avg_inv(mat_H18_CD16_pSS_all,path_cand[zzz])
	alldist_H18_CD19_pSS <- func_ks_avg_inv(mat_H18_CD19_pSS_all,path_cand[zzz])

	# save result in the same format of table
	files <- ls()[ grep('ksdist_H18',ls()) ]
	temp <- as.data.frame( matrix(NA,nrow=length(files),ncol=7), stringsAsFactors=F )
	colnames(temp) <- c('set','hr','typ','don','med','sd','mu')
	for(z in 1:length(files)) {
		tmp <- get(files[z])
		temp$set <- setname
		temp$med[z] <- tmp[[1]]
		temp$mu[z] <- tmp[[2]]
		temp$sd[z] <- tmp[[3]]
		tmp <- unlist( strsplit(files[z],split='_',fixed=T) )
		temp$hr[z] <- tmp[2]
		temp$typ[z] <- tmp[3]
		temp$don[z] <- tmp[4]
		temp$stg[z] <- tmp[5]
	}
	temp <- temp[,c(1,2,3,4,5,8)]
	colnames(temp) <- c('Geneset','Hour','Celltype','Donor','Regulatory','Stage')
	temp$Stage[which(is.na(temp$Stage))] <- 'HD'
	#write.table(temp,file=paste(setname,'_KS_.txt',sep=''),sep='\t',col.names=T,row.names=F)
	for(i in 1:length(dons)) {
		for(j in 1:length(typs)) {
			hrs <- 'H18'
			case <- 'Stg'
			ty <- typs[j]
			do <- dons[i]
			xs <- temp$med[which( (temp$don == dons[i]) & (temp$typ == typs[j]) ) ]
			xh <- temp$med[which( (temp$don == 'HD') & (temp$typ == typs[j]) ) ]
			xs <- append(xs,xh)
			sds <- sd(xs)
			log <- paste(c(setname,hrs,case,do,ty,sds),collapse='\t')
			write(log,file=outlog,append=T) #variance for futher visualization
		}
	}
	rm(list=ls()[ grep('ksdist_H18',ls()) ])

	files <- ls()[ grep('ksdist_H00',ls()) ]
	temp <- as.data.frame( matrix(NA,nrow=length(files),ncol=7), stringsAsFactors=F )
	colnames(temp) <- c('set','hr','typ','don','med','sd','mu')
	for(z in 1:length(files)) {
		tmp <- get(files[z])
		temp$set <- setname
		temp$med[z] <- tmp[[1]]
		temp$mu[z] <- tmp[[2]]
		temp$sd[z] <- tmp[[3]]
		tmp <- unlist( strsplit(files[z],split='_',fixed=T) )
		temp$hr[z] <- tmp[2]
		temp$typ[z] <- tmp[3]
		temp$don[z] <- tmp[4]
		temp$stg[z] <- tmp[5]
	}
	temp <- temp[,c(1,2,3,4,5,8)]
	colnames(temp) <- c('Geneset','Hour','Celltype','Donor','Regulatory','Stage')
	temp$Stage[which(is.na(temp$Stage))] <- 'HD'
	#write.table(temp,file=paste(setname,'_KS_.txt',sep=''),sep='\t',col.names=T,row.names=F)
	for(i in 1:length(dons)) {
		for(j in 1:length(typs)) {
			hrs <- 'H00'
			case <- 'Stg'
			ty <- typs[j]
			do <- dons[i]
			xs <- temp$med[which( (temp$don == dons[i]) & (temp$typ == typs[j]) ) ]
			xh <- temp$med[which( (temp$don == 'HD') & (temp$typ == typs[j]) ) ]
			xs <- append(xs,xh)
			sds <- sd(xs)
			log <- paste(c(setname,hrs,case,do,ty,sds),collapse='\t')
			write(log,file=outlog,append=T) #variance for futher visualization
		}
	}
	rm(list=ls()[ grep('ksdist_H00',ls()) ])
}
# outlog is "For visualization, the F3A". Selecting variable candidate which = highly ranked (sensitivity on stage change) sets
# picked 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY','KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION','REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_BETA_CELLS'

##### EOF Regulatory (Enrichment level) of samples #####





##### Dissecting Co-funtional candidates #####
# Batch process
#args <- commandArgs()
#batch_number=as.character(args[9]) #Batch num
#i0=as.character(args[10]) #set i0
#i1=as.character(args[11]) #set i1

library(plyr)
library(RColorBrewer)
library(pheatmap)

decimal <- 5 #size of permutation (=P-value decimal)
# Input information of the Index
ind <- table_disease_index #must includes variables; don=Donor,typ=Celltype,stg=Stage or Index
# Input value of the Gene Expression
gmt <- anz_00 #Numeric table of Expression
# Input genesets(=perturbation)
tmp <- c('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING','REACTOME_INTERFERON_GAMMA_SIGNALING')
path_ifn <- path_set[unlist(lapply(1:length(tmp), function(x) which(tmp[x] == names(path_set) )   ))]
cand <- path_ifn[[1]]

# Each specific target of samples. 
# Low-stage
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD4',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD4',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD8',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD8',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD16',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD16',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD19',stg='Low')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD19',stg='Low')
# Mid-Stage
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD4',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD4',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD8',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD8',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD16',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD16',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD19',stg='Mid')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD19',stg='Mid')
#High-Stage
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD4',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD4',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD8',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD8',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD16',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD16',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='pSS',typ='CD19',stg='High')
func_coexpr_main(gmt=gmt,cand=cand,ind=ind,don='SLE',typ='CD19',stg='High')

##### EOF Dissecting Co-funtional candidates #####











