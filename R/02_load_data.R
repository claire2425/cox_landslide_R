# 02_load_data.R
# Step 2: load the 5 raster layers + 3 shapefiles, line up the grids, and match the CRS.

cat("\n--- 02_load_data: read rasters and shapefiles ---\n")

suppressPackageStartupMessages({ library(terra); library(sf) })
if (!exists("align_to")) source("R/00_functions.R")

# ---- paths ----
# edit data_dir here if the folder name is different on your machine
data_dir<- "data/LSM Tutorial"
if (!dir.exists(data_dir)) {
  stop(sprintf(paste0(
    "Data directory not found: '%s'\n",
    "Unzip the LSM Tutorial data into ./data/ so this path exists.\n",
    "Tip: list what IS there with  list.dirs('data')"), data_dir))
}

# ---- load the 5 raster layers ----
rast_files<- c(
  Slope     = "Slope.tif",
  LandCover = "Land Cover.tif",
  NDVI      = "NDVI.tif",
  SoilType  = "Soil Type.tif",
  DistRoad  = "Distance to road (m).tif"
)

# check every file is there first -- gives a clearer error than letting rast() fail
missing_files<- rast_files[!file.exists(file.path(data_dir, rast_files))]
if (length(missing_files) > 0) {
  stop(sprintf("These raster files are missing in '%s':\n  %s\n%s",
               data_dir, paste(missing_files, collapse = "\n  "),
               "Check exact filenames with  list.files(data_dir)"))
}

cat("\nReading raster layers...\n")
rasters<- lapply(rast_files, function(f) rast(file.path(data_dir, f)))
names(rasters)<- names(rast_files)

cat("\nRaster summary:\n")
for (nm in names(rasters)) {
  r<- rasters[[nm]]
  v<- values(r, na.rm = TRUE)
  cat(sprintf("  %-10s  %5d x %5d  res=%.1f  CRS=%s  range=[%.2f, %.2f]\n",
              nm, nrow(r), ncol(r), res(r)[1],
              substr(crs(r, describe = TRUE)$name, 1, 22),
              min(v), max(v)))
}

# ---- align every raster to a common grid, using Slope as the reference ----
ref<- rasters$Slope
cat("\nChecking grid alignment against Slope...\n")
for (nm in setdiff(names(rasters), "Slope")) {
  categorical<- nm %in% c("LandCover", "SoilType")   # these two are categorical, so resample with 'near'
  before<- compareGeom(rasters[[nm]], ref, stopOnError = FALSE, messages = FALSE)
  rasters[[nm]]<- align_to(rasters[[nm]], ref, categorical = categorical)
  cat(sprintf("  %-10s  %s\n", nm,
              if (isTRUE(before)) "already aligned"
              else sprintf("RESAMPLED (%s)", if (categorical) "near" else "bilinear")))
}

# ---- load the shapefiles ----
cat("\nReading shapefiles...\n")
shp<- function(name) {
  p<- file.path(data_dir, name)
  if (!file.exists(p)) stop(sprintf("Shapefile not found: %s", p))
  st_read(p, quiet = TRUE)
}
boundary<- shp("Cox's_Bazar_Boundary.shp")
inventory<- shp("Landslide_Inventory.shp")
camp<- shp("Rohingya_camp.shp")

cat(sprintf("  Boundary:     %d feature(s)\n", nrow(boundary)))
cat(sprintf("  Inventory:    %d feature(s)  geom=%s\n",
            nrow(inventory), as.character(st_geometry_type(inventory))[1]))
cat(sprintf("  RohingyaCamp: %d feature(s)\n", nrow(camp)))

# ---- reproject the shapefiles to the raster CRS, but only if they're genuinely different ----
# (using the object-based same_crs() check, since comparing WKT strings directly is fragile)
ref_crs<- crs(ref)
to_fix<- c(
  boundary  = !same_crs(boundary,  ref_crs),
  inventory = !same_crs(inventory, ref_crs),
  camp      = !same_crs(camp,      ref_crs)
)
if (any(to_fix)) {
  cat(sprintf("\nReprojecting to raster CRS: %s\n",
              paste(names(to_fix)[to_fix], collapse = ", ")))
  if (to_fix["boundary"])  boundary<- st_transform(boundary,  ref_crs)
  if (to_fix["inventory"]) inventory<- st_transform(inventory, ref_crs)
  if (to_fix["camp"])      camp<- st_transform(camp,      ref_crs)
} else {
  cat("\nAll shapefiles already share the raster CRS (no reprojection).\n")
}

# ---- save (terra rasters have to be wrap()ped before they can be serialised) ----
if (!dir.exists("output")) dir.create("output", recursive = TRUE)
saveRDS(lapply(rasters, wrap), "output/rasters_raw.rds")
saveRDS(boundary,  "output/boundary.rds")
saveRDS(inventory, "output/inventory.rds")
saveRDS(camp,      "output/camp.rds")

cat("\nSaved: output/rasters_raw.rds, boundary.rds, inventory.rds, camp.rds\n\n")
