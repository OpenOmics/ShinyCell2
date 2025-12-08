#' Generate data files required for shiny app (spatial data)
#'
#' Generate data files required for shiny app, specifically spatial data
#' A prefix is specified for each set of files to allow for multiple 
#' single-cell datasets in a single Shiny app.
#'
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param scConf shinycell config data.table
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate hdf5r
#'
#' @examples
#' makeShinyFilesSpatial(seu, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyFilesSpatial <- function(
  obj, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/"){
  ### Preprocessing and checks
  if(!file.exists(paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))){
    stop(paste0(shiny.prefix, "conf.rds file is missing! Have you ran makeShinyFilesGEX?"))
  }
  
  ### Start extraction
  sc1image = list()
  if(class(obj)[1] == "Seurat"){
    if(.hasSlot(obj, "images")){
      # coordinates - handle both VisiumV1 and VisiumV2 formats
      visium_image <- obj@images[[1]]
      if ("VisiumV2" %in% class(visium_image)) {
        # New VisiumV2 format - coordinates are in boundaries$centroids
        coords <- visium_image@boundaries$centroids@coords
        # Convert to expected format with imagecol and imagerow columns
        sc1image$coord <- data.frame(
          tissue = 1,  # Default tissue value
          row = NA,    # Not available in VisiumV2
          col = NA,    # Not available in VisiumV2  
          imagecol = coords[, "x"],
          imagerow = coords[, "y"]
        )
        rownames(sc1image$coord) <- rownames(coords)
      } else {
        # Traditional VisiumV1 format
        sc1image$coord <- visium_image@coordinates
      }
      
      # uncropped background image
      bg_image <- GetImage(obj, mode = "raster") # background image
      # bg_grob <- rasterGrob(bg_image, width=unit(1,"npc"), height=unit(1,"npc"), 
      #                       interpolate = FALSE) # background image grob
      sc1image$bg_image <- bg_image
      # sc1image$bg_grob  <- bg_grob
      
      # for cropped image generation in the server.R script
      sc1image$lowres   <- visium_image@scale.factors$lowres

    }
  } else {
    stop("Only Seurat objects are accepted!")
  }
  
  ### Saving objects
  saveRDS(sc1image, file = paste0(shiny.dir, "/", shiny.prefix, "image.rds"))
}


