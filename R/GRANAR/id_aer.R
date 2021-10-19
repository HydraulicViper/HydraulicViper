

# Function to get the id_aerenchyma from curent_root.xml

id_aer <- function(root_path = "./MECHA/cellsetdata/root_1.xml", geometry_path = "./MECHA/Projects/GRANAR/in/Maize_Geometry.xml"){
  require(xml2)  
  
  x <- read_xml(root_path)
  sh_path <- unlist(str_split(root_path, "_"))[1]
  root_id <- parse_number(unlist(str_split(root_path, "_"))[2])+2320
  
  cells <- xml_children(xml_find_all(x, "//cells"))
  id_aer <- which(xml_attr(cells, "group") == "4")-1
  
  # aer in geom xml

  if (is.null(geometry_path)) {
    warning("No path specified")
  }else{
    xml <- read_xml(geometry_path)
    
    aer <- xml_children(xml_find_all(xml, "//aerenchyma_range"))
  }
  if(length(id_aer) < 1){
    warning("No aerenchyma id specified")
  }else{
    id_aerenchyma <- id_aer
    
    # newbee <- 'aerenchyma id="0"'
    new_siblings <- paste0('aerenchyma id="',id_aerenchyma,'"')
    
    xml_add_sibling(aer, new_siblings)
  }
  
  xml_remove(aer[1])
  
  path <- paste0(c(unlist(str_split(geometry_path, ".xml"))[1]),"_aer_",root_id,".xml")
  
  write_xml(xml, path)
  fc <- file.copy(root_path, paste0(sh_path,"_",root_id,".xml"))
  
  return(TRUE)
  
  
}
  
  