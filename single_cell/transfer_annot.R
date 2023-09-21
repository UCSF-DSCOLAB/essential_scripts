# Helper function to transfer annotations from a previous set of labels to 
# a new cluster using max per cluster

library(Seurat)
library(tidyverse)

transfer_annot = function(sobj, new_cluster, prev_label, new_cluster_label, min.frac=0.4){
  meta0 = sobj@meta.data  %>% 
    as_tibble(rownames="cell_id") %>%
    dplyr::select(cell_id, {{new_cluster}}, {{prev_label}}) %>%
    dplyr::rename(cluster={{new_cluster}}, orig_label={{prev_label}})
  
  clus_to_lab = meta0 %>%
    filter(!is.na(cluster), !is.na(orig_label)) %>%
    group_by(cluster) %>%
    dplyr::mutate(tot=n()) %>%
    group_by(cluster, orig_label, tot) %>%
    dplyr::count() %>%
    arrange(cluster, desc(n)) %>%
    mutate(frac=n/tot)
  
  clus_to_lab2 = clus_to_lab %>%
    ungroup() %>%
    group_by(cluster) %>%
    slice_max(order_by=frac, n=1)
  
  clus_to_lab3 = clus_to_lab2 %>%
    mutate(orig_label=as.character(orig_label)) %>%
    mutate(assigned_label=ifelse(frac < min.frac | is.na(orig_label), sprintf("%s*", orig_label), 
                                 orig_label)) %>%
    arrange(assigned_label) %>%
    ungroup() %>%
    group_by(assigned_label) %>%
    dplyr::mutate(idx=1:n(),
           nlab=n()) %>%
    dplyr::mutate(new_label := ifelse(nlab ==1, assigned_label,
                               paste(assigned_label, idx))) %>%
    ungroup()
  
  clus_map = clus_to_lab3 %>% dplyr::select(cluster, new_label) %>% distinct()
  
  sobj = AddMetaData(sobj, meta0 %>% 
                       left_join(clus_map, by="cluster") %>% 
                       pull(new_label), new_cluster_label)
  return(sobj)
}
