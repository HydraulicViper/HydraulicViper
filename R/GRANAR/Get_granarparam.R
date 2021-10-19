

Get_granarparam <- function(Data_params){
  
  Data_params <- Data_params%>%
    dplyr::group_by(Image, tissue, dist, Root, Root_id)%>%
    dplyr::summarise(cell = mean(cell, na.rm = T),
                     layer = mean(layer, na.rm = T))%>%
    ungroup()%>%
    filter(!is.na(tissue))
  
  tiss_param <- NULL
  for(tiss in unique(na.omit(Data_params$tissue[Data_params$tissue != "proto"]))){
    
    fit<- lm(cell~dist, data = Data_params%>%filter(tissue == tiss))
    summary(fit)
    tis_temp <- Data_params%>%
      filter(tissue == tiss)%>%
      mutate(reg = cell-(fit$coefficients[1]+fit$coefficients[2]*dist))
    
    tmp <- tibble(one = fit$coefficients[1], 
                  slope = fit$coefficients[2], 
                  tiss, 
                  p = summary(fit)$coefficients[2,4],
                  va = var(tis_temp$reg, na.rm = T))
    if(summary(fit)$coefficients[2,4] > 0.05){
      tmp$slope <- 0
      tmp$one = median(Data_params$cell[Data_params$tissue == tiss])
      tmp$va = var(Data_params$cell[Data_params$tissue == tiss], na.rm = T)
    }
    tiss_param <- rbind(tiss_param, tmp)
  }
  layer_param <- NULL
  for(lay in c("stele", "cortex", "xylem", "proto")) {
    fit<- lm(layer~dist, data = Data_params%>%filter(tissue == lay))
    summary(fit)
    lay_temp <- Data_params%>%
      filter(tissue == lay)%>%
      mutate(reg = layer-(fit$coefficients[1]+fit$coefficients[2]*dist))
    
    tmp <- tibble(one = fit$coefficients[1], 
                  slope = fit$coefficients[2], 
                  lay, 
                  r = summary(fit)$r.squared,
                  va = var(lay_temp$reg, na.rm = T))
    
    
    ok <- summary(fit)$coefficients[2,4]
    if(is.na(summary(fit)$coefficients[2,4])){
      ok <- 1
      tmp$one <- median(Data_params$layer[Data_params$tissue == lay], na.rm = T)
      tmp$va <- var(Data_params$layer[Data_params$tissue == lay], na.rm = T)
      message(paste0("No enough data to make a regression for ",lay))
      
    }
    if(ok > 0.05){
      tmp$slope <- 0
      tmp$one = median(Data_params$layer[Data_params$tissue == lay], na.rm = T)
      tmp$va <- var(Data_params$layer[Data_params$tissue == lay], na.rm = T)
    }
    layer_param <- rbind(layer_param, tmp)
  }
  
  for (im in unique(Data_params$Image)) {
    # Get the total diameter of the root
    
    tmp <- Data_params%>%
      filter(Image == im)
    if(nrow(tmp) == 8){
      
      tot <- tmp$cell[tmp$tissue == "endo"]+
        tmp$cell[tmp$tissue == "pericycle"]+
        tmp$cell[tmp$tissue == "exo"]+
        tmp$cell[tmp$tissue == "epi"]+
        tmp$layer[tmp$tissue == "stele"]/2+
        tmp$layer[tmp$tissue == "cortex"]
      
      total <- tibble(Image = im, 
                      tissue = c("total"),
                      cell = NA, 
                      layer = tot, 
                      dist = unique(tmp$dist), 
                      Root = unique(tmp$Root),
                      Root_id = unique(tmp$Root_id))
      
      Data_params <- rbind(Data_params, total)
    }
  }
  
  param <- list(tiss_param, layer_param)
  return(param)
}


radial_corr <- function(RS_tot,x1,x2){ 
  # x1 where endodermis is fully suberized
  # x2 where exodermis has a Casparian strip 
  
  RS_tot$Kr_corr[RS_tot$apo_barriers == 1 & RS_tot$dist <= x1] <- RS_tot$kr[RS_tot$apo_barriers == 1 & RS_tot$dist <= x1]
  RS_tot$Kr_corr[RS_tot$apo_barriers == 2 & RS_tot$dist <= x2 & RS_tot$dist > x1 ] <- RS_tot$kr[RS_tot$apo_barriers == 2 & RS_tot$dist <= x2 & RS_tot$dist > x1]
  RS_tot$Kr_corr[RS_tot$apo_barriers == 4 & RS_tot$dist > x2] <- RS_tot$kr[RS_tot$apo_barriers == 4 & RS_tot$dist > x2]
  return(RS_tot)
}