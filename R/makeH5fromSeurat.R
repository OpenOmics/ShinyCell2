#' Make h5 object from Seurat assay data
#'
#' Make h5 object from Seurat assay data
#'
#' @param obj input Seurat (v3+) object
#' @param sc1meta data.table of cell metadata
#' @param filename filename of output h5 file
#' @param gex.assay assay in Seurat object to use
#' @param gex.slot slot in single-cell assay to use
#' @param chunkSize number of genes written to h5file at any one time
#'
#' @return h5 object
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate
#'
#' @export
makeH5fromSeurat <- function(obj, sc1meta, filename, 
                             gex.assay, gex.slot, chunkSize){
  # DEBUG: Print input parameters
  cat("=== DEBUG: makeH5fromSeurat parameters ===\n")
  cat("gex.assay:", gex.assay, "\n")
  cat("gex.slot:", gex.slot, "\n")
  cat("filename:", filename, "\n")
  cat("chunkSize:", chunkSize, "\n")
  cat("Number of cells in sc1meta:", nrow(sc1meta), "\n")
  cat("First few cellIDs in sc1meta:", head(sc1meta$cellID, 10), "\n")
  
  # Create h5 file and get ready
  if(class(obj@assays[[gex.assay]]) == "Assay5"){
    cat("Assay type: Assay5 (Seurat v5)\n")
    cat("Available layers in assay:", names(obj@assays[[gex.assay]]@layers), "\n")
    gex.matdim = dim(obj@assays[[gex.assay]]@layers[[gex.slot]])
    cat("Matrix dimensions:", gex.matdim, "\n")
    # Get available cells in this assay
    assay.cells = colnames(obj@assays[[gex.assay]]@layers[[gex.slot]])
  }else{
    cat("Assay type: Legacy Assay (Seurat v3/v4)\n")
    cat("Available slots in assay:", slotNames(obj@assays[[gex.assay]]), "\n")
    gex.matdim = dim(slot(obj@assays[[gex.assay]], gex.slot))
    cat("Matrix dimensions:", gex.matdim, "\n")
    # Get available cells in this assay
    assay.cells = colnames(slot(obj@assays[[gex.assay]], gex.slot))
  }
  
  cat("Number of cells in assay:", length(assay.cells), "\n")
  cat("First few cellIDs in assay:", head(assay.cells, 10), "\n")
  
  # Filter sc1meta to only include cells present in this assay
  sc1meta.filtered = sc1meta[sc1meta$cellID %in% assay.cells, ]
  cat("Number of cells after filtering:", nrow(sc1meta.filtered), "\n")
  
  # Check overlap
  overlap_count = sum(sc1meta$cellID %in% assay.cells)
  cat("Overlap between sc1meta and assay:", overlap_count, "\n")
  
  # ADDED: Verify all cellIDs in filtered metadata exist in assay
  missing_cells = setdiff(sc1meta.filtered$cellID, assay.cells)
  if(length(missing_cells) > 0){
    warning(paste("Removing", length(missing_cells), "cells not found in assay"))
    cat("First few missing cells:", head(missing_cells, 10), "\n")
    sc1meta.filtered = sc1meta.filtered[!sc1meta.filtered$cellID %in% missing_cells, ]
  }
  
  # ADDED: Check if we have any cells left
  if(nrow(sc1meta.filtered) == 0){
    cat("ERROR: No overlapping cells found!\n")
    cat("Sample cellIDs from metadata:", head(unique(sc1meta$cellID), 20), "\n")
    cat("Sample cellIDs from assay:", head(unique(assay.cells), 20), "\n")
    stop("No overlapping cells found between metadata and assay data")
  }
  
  cat("Final number of cells to process:", nrow(sc1meta.filtered), "\n")
  cat("==========================================\n")
  
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",  dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = c(gex.matdim[1], nrow(sc1meta.filtered)), maxdims = c(gex.matdim[1], nrow(sc1meta.filtered))),
    chunk_dims = c(1, nrow(sc1meta.filtered)))
  chk = chunkSize
  while(chk > (gex.matdim[1]-8)){
    chk = floor(chk / 2)     # Account for cases where nGene < chunkSize
  } 
  
  # Start writing to file
  nChunk = floor((gex.matdim[1]-8)/chk)
  if(class(obj@assays[[gex.assay]]) == "Assay5"){
    # Seurat v5
    for(i in 1:nChunk){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        obj[[gex.assay]][gex.slot][
          ((i-1)*chk+1):(i*chk), sc1meta.filtered$cellID])
    }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        obj[[gex.assay]][gex.slot][
          (i*chk+1):gex.matdim[1], sc1meta.filtered$cellID])
      gex.rownm = rownames(obj[[gex.assay]][gex.slot])
      
  } else {
    for(i in 1:nChunk){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        slot(obj@assays[[gex.assay]], gex.slot)[
          ((i-1)*chk+1):(i*chk), sc1meta.filtered$cellID])
    }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        slot(obj@assays[[gex.assay]], gex.slot)[
          (i*chk+1):gex.matdim[1], sc1meta.filtered$cellID])
      gex.rownm = rownames(slot(obj@assays[[gex.assay]], gex.slot))
  }
  sc1gexpr$close_all()
  # Output gene names
  return(gex.rownm)
}