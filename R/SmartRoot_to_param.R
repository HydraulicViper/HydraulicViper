

Root_data <- read_csv("Global_data.csv") %>%
  as.tibble(Root_data)%>%
  mutate(Id_plant = substr(image, 1, 2), # Get tag name of the plant
         Treatment = substr(image,9,9)) # Get tag treatment


age <- 18 # Age of the root system


is.na(Root_data) <- Root_data == "null"
Root_data[["insertion_first_child"]] <- as.numeric(Root_data[["insertion_first_child"]])
Root_data[["insertion_last_child"]] <- as.numeric(Root_data[["insertion_last_child"]])
is.na(Root_data$child_density) <- !Root_data$child_density

#-----------------------Compute the rparam files  ---------------------------------------------
Root_data <- Root_data%>%
  mutate(dist = NA) # Prepare column to get the information later
for(t in unique(Root_data$Treatment)){ # loop to works one treatment at the time
  temp_tre <- Root_data%>%
    filter(Treatment == t) # filter the Root_data to the selected treatment
  dady <- temp_tre%>% # dplyr:: is to call a precise function in a selected package
    dplyr::group_by(parent)%>% # Groups root that have kids 
    dplyr::summarise(n = n())%>% # How many kids
    ungroup()%>% # ungroup to select only the ones with high resolution 
    filter(n > 3,
           parent != "-1") # Axes cannot be parental root
  
  for (dad in dady$parent) {
    tmp_r <- Root_data%>%
      filter(parent == dad)%>%
      arrange(insertion_position) # order the lat root along the parental root
    tmp_r <- tmp_r%>%
      mutate(dist = insertion_position-c(insertion_position[1],insertion_position[1:(nrow(tmp_r)-1)]))%>% # compute the distance between the root and the previous one
      filter(dist < 1)%>% # Cluster of root close to each other
      select(root, dist) # simplify the data
    
    Root_data$dist[which(Root_data$root %in% tmp_r$root)] <- tmp_r$dist # Merge dist info into Root_data
    
  }
  
}
# ------ Ln ----------------
# Print the distance between laterals for every axial root type
Ln <- Root_data%>%
  dplyr::group_by(parent_name)%>%
  dplyr::summarise(ln = mean(dist, na.rm = TRUE),
                   ln_sd = sd(dist, na.rm = TRUE))%>%
  ungroup()%>%
  filter(parent_name != "-1")



# --------- Lb & La -------------------
# Axes must have attached root, Lb should be less than 4, La should be less than the half of the root length. 
L <- Root_data%>%
  filter(n_child > 0,
         insertion_first_child < 4,
         insertion_last_child > 0.3*length)%>%
  dplyr::group_by(root_name)%>%
  dplyr::summarise(lb = mean(insertion_first_child, na.rm = TRUE),
                   lb_std = sd(insertion_first_child, na.rm = TRUE),
                   la = mean((length - insertion_last_child), na.rm = TRUE),
                   la_std = sd((length - insertion_first_child), na.rm = TRUE))%>%
  ungroup()

# ------------- t ------------
# Insertion angle 
t <- Root_data%>%
  filter(root_name %in% c("Lat", "LongLat"))%>%
  dplyr::group_by(root_name)%>%
  dplyr::summarise(t = mean(insertion_angle, na.rm = TRUE),
                   t_std = sd(insertion_angle/180, na.rm = TRUE))

# ------------ a ------------
a <- Root_data%>%
  dplyr::group_by(root_name)%>%
  dplyr::summarise(a = mean(diameter, na.rm = TRUE),
                   a_std = sd(diameter, na.rm = TRUE))%>%
  ungroup()

# ----------- r ----------------
Axes <- Root_data%>%
  filter(root_name %!in% c("Lat", "LongLat"),
         n_child > 5 # Only works with good resolution data
  )
r_lat_rel <- NULL
for(i in unique(Axes$root)){
  #Li <- Root_data$length[which(Root_data$root == i)]
  temp <- Root_data%>%
    filter(parent == i)
  if(nrow(temp)>3){
    Reg <- lm(length ~ insertion_position, data = temp)
    R <- summary(Reg)$r.squared
    if(R > 0.7 && coefficients(aov(Reg))[2] < 0){
      m_rel = as.numeric(-coefficients(aov(Reg))[2])
      r_laty = cbind(Treatment = unique(temp$Treatment), root_name = unique(temp$parent_name), m = m_rel)
      r_lat_rel = rbind(r_lat_rel, r_laty)    
    }
  }
}
r_lat <-  as.tibble(r_lat_rel)%>%
  mutate(m = as.numeric(m))


# ---------- Lm --------------
Lm <- Root_data%>%
  dplyr::group_by(root_name)%>%
  dplyr::summarise(Lm = mean(length),
                   Lm_std = sd(length),
                   n = n())


