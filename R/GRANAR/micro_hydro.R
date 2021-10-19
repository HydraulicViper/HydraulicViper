

microhydro <- function(path = "MECHA/Projects/granar/in/Maize_Hydraulics.xml", 
                       kw = 2.4E-4, # hydraulic conductivity of standard walls
                       km = 3.0E-5, # Cell membrane permeability, with separate contribution of the biphospholipid layer (km) and AQP (kAQP)
                       kAQP = 4.3E-4, # cm/hPa/d
                       kpl = 5.3E-12){ # Individual plasmodesma conductance
  
  #Load the xml file with the hydraulics input for MECHA
  x <- read_xml(path)
  kw_tmp <- xml_children(x)[2]
  km_tmp <- xml_children(x)[4]
  kAQP_tmp <- xml_children(x)[5]
  kpl_tmp <- xml_children(x)[7]
  temp <- xml_find_all(kw_tmp, ".//kw")
  temp_km <- xml_find_all(km_tmp, ".//km")
  temp_kAQP <- xml_find_all(kAQP_tmp, ".//kAQP")
  temp_kpl <- xml_find_all(kpl_tmp, ".//kpl")
  xml_attr(temp_kpl, "value") <- kpl
  xml_attr(temp_kAQP, "value") <- kAQP
  xml_attr(temp, "value") <- kw
  xml_attr(temp_km, "value") <- km
  write_xml(x, file = path)
  print(paste0("the cell hydraulic parameters have been change"))
  
}

wallthick <- function(path = "MECHA/Projects/GRANAR/in/Maize_Geometry.xml", wall_thickness = 1.5){
  x <- read_xml(path)
  thick <- xml_children(x)[19]
  ver <- xml_name(thick)
  if(ver != "thickness"){
    message("thickness not found in xml")
    thickness
  }
  xml_attr(thick, "value") <- wall_thickness
  write_xml(x, file = path)
  
}
