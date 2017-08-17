#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stats)
if (length(args)==0) {
  stop("At least one argument must be given (input file).n", call.=FALSE)
} 

local_stats<-c()

tree<-read.tree(args[1])
if (file.exists(dm_filename)){
    tree<-read.tree(filename)
    dm_matrix<-cophenetic(tree)
    dm<-as.dist(dm_matrix, diag = FALSE, upper = FALSE)
    HC<-hclust(dm, method = "complete", members = NULL)
    
    for (cf in seq(from = 0.4, to = 20, by = 0.4)){
      myhcl <- cutree(HC, h=cf)
      clusters<-as.data.frame(table(myhcl))$Freq
      clusters<-sort(clusters)
      if (any(clusters==1)){count_single_tip=table(clusters)[names(table(clusters))==1][[1]]
      }else{count_single_tip=0}
      
      if ((count_single_tip<sum(clusters))&count_single_tip>0){
        clusters_no_single_tip<-clusters[-c(1:count_single_tip)]
      }else if (count_single_tip==0){
        clusters_no_single_tip<-clusters[1:length(clusters)]
      }
      if (count_single_tip<sum(clusters)){
        cluster_size=mean(clusters_no_single_tip)
        cluster_number=length(clusters_no_single_tip)
        cluster_var=var(clusters_no_single_tip)
        if (length(clusters_no_single_tip)==1){
          cluster_var=0
        }
        fraction_clustered_tips=sum(clusters_no_single_tip)/(sum(clusters_no_single_tip)+count_single_tip)
        ratio_max_min<-max(clusters_no_single_tip)/min(clusters_no_single_tip)
        local_stats[[(n-1)*50+cf/0.2]]<-c(n, cf, cluster_size, cluster_number, cluster_var, fraction_clustered_tips, ratio_max_min)
        
      }
      
    }
    
  }

statsMatrix<-do.call(rbind,local_stats)
colnames(statsMatrix) <-c("tree_index","cutoff","average_cluster_size","cluster_number","var_cluster_szie","fraction_clustered_tips","max_over_min_size")
write.table(statsMatrix, file=args[2],row.names=FALSE)
