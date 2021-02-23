rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "boxr",
  "raster",
  "tigris",
  "sf",
  "exactextractr",
  "fmsb",
  "ggplot2",
  "ggalt"
)

install.packages(setdiff(packages, rownames(installed.packages()))) 
lapply(packages, require, character.only = TRUE)

#====================================================================================================
# Authenticate BOX connection using Client ID and Secret in renviron

source("box_auth.R")

#====================================================================================================

# load TNC tillage probability rasters
probList = list(
  "https://s3.amazonaws.com/TillageModel/tillage/model/Iowa/Iowa_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/Minnesota/Minnesota_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/Montana/Montana_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/Nebraska/Nebraska_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/North+Dakota/North+Dakota_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/South+Dakota/South+Dakota_probs.tif",
  "https://s3.amazonaws.com/TillageModel/tillage/model/Wyoming/Wyoming_probs.tif"
)

# if not present already, download the usxp s35 mtr raster from BOX
#if(!file.exists('s35_mtr.tif')){
  
  box_dl(file_id = '752345336503', overwrite = T, file_name = 's35_mtr.tif')
  
#}

# load the usxp s35 mtr raster (5 land use change classes)
r_mtr = raster("s35_mtr.tif")

# get a list of states for which probability rasters are available
states = lapply(probList, function(e)  gsub("[+]", " ", unlist(strsplit(e, '/'))[length(unlist(strsplit(e, '/')))-1]))

# get shapefiles of all counties within all states with data
counties = lapply(states, function(s) st_as_sf(counties(state = s, cb = FALSE, resolution = "500k", year = 2010)))
names(counties) = states
counties = do.call(rbind, counties)
counties['STATENAME'] = do.call(rbind, strsplit(row.names(counties), '[.]'))[,1]

# # calculating the number of pixels  of conversion (mtr = 3) in each county
# counties['nPixels_mtr3'] = exact_extract(
#   x = r_mtr,
#   y = st_transform(counties, st_crs(r_mtr)),
#   fun = function(value, cov_frac) sum(ifelse(value == 3, 1,0)*cov_frac)
# )
# 
# # calculate total count of conversion pixels in the region by summing those within each county.
# # (will be a non-integer because used coverage_fraction in county level calculation)
# total_nPixels_mtr3 = sum(counties$nPixels_mtr3, na.rm = T)

file = probList[[1]]
tp_rasters = lapply(probList, function(file){
  
  # create temporary directory to offload state's frequency tables.
  if(!dir.exists('temp_state_freqs')){
    dir.create('temp_state_freqs')
  }

  # get the name of the state's tillage probability raster
  tp_name = unlist(strsplit(file, '/'))[length(unlist(strsplit(file, '/')))]
  print(tp_name)

  # if file hasn'y yet been dowloaded
  #if(!file.exists(tp_name)){
    download.file(file, tp_name, mode="wb")
  #}

  # load  the tillage probability raster
  #return(raster(tp_name))
  r_tp = raster(tp_name)
  
  # get counties of state
  state_counties = counties[counties$STATENAME == gsub('_probs.tif', '', tp_name),]
  
  # iterate through each county (to reduce computational demand in comparison to a statewide calc)
  for(i in 1:nrow(state_counties)){
    
    # explicity cast the county's polygon as an sf multipolygon
    county = st_cast(state_counties[i,], "MULTIPOLYGON")
    
    # crop and mask tillage raster to the county's extent (ensuring projections match)
    r_tp_clip = crop(r_tp, st_transform(county, st_crs(r_tp)))
    r_tp_clip = mask(r_tp_clip, st_transform(county, st_crs(r_tp)))
    
    # crop mtr raster to the county's extent (ensuring projections match)
    r_mtr_clip = crop(r_mtr, st_transform(county, st_crs(r_mtr)))
    r_mtr_clip = mask(r_mtr_clip , st_transform(county, st_crs(r_mtr)))
    
    # reproject the cropped mtr raster to match the tillage raster's projection (using nearest neighbor resampling)
    r_mtr_clip_reproj = projectRaster(r_mtr_clip, r_tp_clip, method = 'ngb')
    
    # mask both rasters to each others data extent
    r_mtr_clip_reproj = mask(r_mtr_clip_reproj, is.na(r_tp_clip), maskvalue = TRUE)
    r_tp_clip = mask(r_tp_clip, is.na(r_mtr_clip_reproj), maskvalue = TRUE)
    
    # get mask of stable cropland in 2008 (i.e. mtr = 2 (stable cropland), mtr = 4 (abandoned cropland), and mtr = 5 (intermittent cropland))
    r_cropland_2008 = reclassify(r_mtr_clip_reproj, rcl = cbind(c(1,2,3,4,5), c(0,1,0,1,1)))
    
    # mask out 2008 cropland from the tillage probability raster
    r_tp_clip_noncropland = mask(r_tp_clip, r_cropland_2008, maskvalue = 1)
    r_tp_clip_noncropland[r_tp_clip_noncropland == 0] = NA
    
    # get frequency table of (scaled) tillage probability values
    freqs = as.data.frame(freq(r_tp_clip_noncropland*(10^40)))
    
    # count the number of remaining (after masking) conversion pixels in each county
    mtr_count = as.data.frame(freq(r_mtr_clip_reproj))
    mtr_count = mtr_count[complete.cases(mtr_count),]
    conversionPixels = mtr_count$count[mtr_count$value == 3]
    
    if(i == 1){
      all_freqs = freqs
      
      all_conversionPixels = conversionPixels
    }else{
      joined_freqs = merge(all_freqs, freqs, by = "value", all = T)
      
      all_freqs = data.frame('value' = joined_freqs$value, 
                             'count' = rowSums(joined_freqs[,grep('count', colnames(joined_freqs))], na.rm = T))
      
      all_conversionPixels = all_conversionPixels = conversionPixels
    }
    print(i)
  }
  
  write.csv(all_freqs, file = paste0('temp_state_freqs/',gsub('_probs.tif', '.csv', tp_name)), row.names = F)
  
  return(all_conversionPixels)
  
  rm(list = c("all_freqs", "freqs"))
  gc()
  
})





























r = tp_rasters[[1]]
r_reproj = projectRaster(from = r, to = r_mtr)
e = st_as_sfc(st_bbox(r))


# explicity cast the county's polygon as an sf multipolygon
e = st_cast(e, "MULTIPOLYGON")

# crop and mask tillage raster to the county's extent (ensuring projections match)
r_tp_clip = crop(r_tp, st_transform(county, st_crs(r_tp)))
r_tp_clip = mask(r_tp_clip, st_transform(county, st_crs(r_tp)))
# # align all rasters to same projection for mosaicing
# # (from: https://gis.stackexchange.com/questions/224781/merge-rasters-with-different-origins-in-r)
# tp_rasters = lapply(tp_rasters, function(r){
#   
#   template = projectRaster(from = r, to = r_mtr, alignOnly = T)
#   
#   return(projectRaster(from = r, to= template))
#   
# })


s = brick(tp_rasters)

s = stack(tp_rasters)
merge(tp_rasters[[1]], tp_rasters[[2]])
origin(tp_rasters[[1]])
origin(tp_rasters[[2]])
# mosaic rasters
rasterOptions(tolerance = 0.1)
tp_rasters$fun = max # if multiple estimates available for a pixel, choose max
r_tp = do.call(mosaic, tp_rasters)

# reproject r_tp to match projection of r_mtr
r_tp = projectRaster(r_tp, crs = crs(r_mtr))

# clip mtr to same extent as r_tp
r_mtr_clip = crop(r_mtr, extent(r_tp))
r_mtr_clip = mask(r_mtr_clip,extent(r_tp))


# mask out croplands in 2008


# for first state, extract all pixel values, rank in decending order, remove those in positions greater than nPixel
# for subsequent states, extract pixels, combine with previous pixels. Reorder, and remove those in positions greater than nPixel
file = probList[[1]]
for(file in probList){
  
  # get the name of the state's tillage probability raster
  tp_name = unlist(strsplit(file, '/'))[length(unlist(strsplit(file, '/')))]
  
  download.file(file, tp_name, mode="wb")
  
  # load  the tillage probability raster
  r_tp = raster(tp_name)
  
  # if tp projection does not match mtr, reproject tp
  if(as.character(crs(r_tp)) != as.character(crs(r_mtr))){
    r_tp = projectRaster(r_tp, crs = crs(r_mtr))
  }
  
  # MASK MTR RASTER TO ONLY INCLUDE PIXELS FOR WHICH TP RASTER HAS DATA.
  r_mtr_wData = mask(r_mtr, is.na(r_tp), maskvalue = TRUE)
    
    
  # get mask of stable cropland in 2008 (i.e. mtr = 2 (stable cropland), mtr = 4 (abandoned cropland), and mtr = 5 (itermittent cropland))
  r_cropland_2008 = reclassify(r_mtr_clip_reproj, rcl = cbind(c(1,2,3,4,5), c(0,1,0,1,1)))
  
  # mask out 2008 cropland from the tillage probability raster
  r_tp_clip_noncropland = mask(r_tp_clip, r_cropland_2008, maskvalue = 1)
  
  # extract and sort the remaining tillage probability pixel values in order of decreasing probability
  probs = getValues(r_tp_clip_noncropland)
  probs_sort = sort(probs, decreasing = T)
  
  
}







# iterate through each tillage probability raster
for(file in probList){
  
  # get the name of the state's tillage probability raster
  tp_name = unlist(strsplit(file, '/'))[length(unlist(strsplit(file, '/')))]
  
  download.file(file, tp_name, mode="wb")
  
  # load  the tillage probability raster
  r_tp = raster(tp_name)
  
  # get high res census county boundaries
  counties = st_as_sf(counties(state = gsub("[+]", " ",gsub("_probs.tif", "", tp_name)), cb = FALSE, resolution = "500k", year = 2010))
  
  # iterate through each of a state's counties
  countyList = list()
  
  for(i in 1:nrow(counties)){
    
    countyList[[i]] = tryCatch({
      
      # explicity cast the county's polygon as an sf multipolygon
      county = st_cast(counties[i,], "MULTIPOLYGON")
      
      # crop and mask tillage raster to the county's extent (ensuring projections match)
      r_tp_clip = crop(r_tp, st_transform(county, st_crs(r_tp)))
      r_tp_clip = mask(r_tp_clip, st_transform(county, st_crs(r_tp)))
      
      # crop mtr raster to the county's extent (ensuring projections match)
      r_mtr_clip = crop(r_mtr, st_transform(county, st_crs(r_mtr)))
      r_mtr_clip = mask(r_mtr_clip , st_transform(county, st_crs(r_mtr)))
      
      # reproject the cropped mtr raster to match the tillage raster's projection (using nearest neibhor resampling)
      r_mtr_clip_reproj = projectRaster(r_mtr_clip, r_tp_clip, method = 'ngb')
      
      # get the number of s35 conversion (2008-2016) pixels 
      f_mtr = as.data.frame(freq(r_mtr_clip_reproj))
      nPixels = as.numeric(na.omit(f_mtr$count[f_mtr$value == 3]))
      
      # if there's no conversion in a county, return just state and coutny fps
      if(length(nPixels) == 0){
        
        c("STATEFP" = county$STATEFP[1], "COUNTYFP" = county$COUNTYFP[1])
        
      }else{
        
        # get mask of stable cropland in 2008 (i.e. mtr = 2 (stable cropland), mtr = 4 (abandoned cropland), and mtr = 5 (itermittent cropland))
        r_cropland_2008 = reclassify(r_mtr_clip_reproj, rcl = cbind(c(1,2,3,4,5), c(0,1,0,1,1)))
        
        # mask out 2008 cropland from the tillage probability raster
        r_tp_clip_noncropland = mask(r_tp_clip, r_cropland_2008, maskvalue = 1)
        
        # extract and sort the remaining tillage probability pixel values in order of decreasing probability
        probs = getValues(r_tp_clip_noncropland)
        probs_sort = sort(probs, decreasing = T)
        
        # select values in rank order until their count == nPixels
        probs_selected = probs_sort[1:nPixels]
        
        # identify the minimum probability in the selected set
        minProb = min(probs_selected)
        
        # count the number of pixels having the minimum probability in the selected set
        count_minProb_selected = length(which(probs_selected == minProb))
        
        # count the total number of pixels (in the county's full set) having the minimum probability
        count_minProb_total = length(which(probs_sort == minProb))
        
        # calculate a weigth that will be applied to pixels in cases where minProb is both included and excluded from the selected set.
        minProb_weight = count_minProb_selected/count_minProb_total
        
        # get a binary 1/0 presence/absence raster of most probable conversion (based on minProb)
        r_probable_exp = mask(r_tp_clip_noncropland, r_tp_clip_noncropland < minProb, maskvalue = TRUE)
        r_probable_exp[r_probable_exp > 0] = 1
        r_probable_exp[is.na(r_probable_exp[])] = 0
        r_probable_exp = mask(r_probable_exp , st_transform(county, st_crs(r_probable_exp)))
        
        # get a binary 1/0 presence/absence raster of cropland expansion (2008-16)
        r_observed_exp = mask(r_mtr_clip_reproj, r_mtr_clip_reproj != 3, maskvalue = TRUE)
        r_observed_exp[r_observed_exp == 3] = 1
        r_observed_exp[is.na(r_observed_exp[])] = 0 
        r_observed_exp = mask(r_observed_exp , st_transform(county, st_crs(r_observed_exp)))
        
        
        # if the minProb_wheight == 1, no-need for weighting so calculate the comparison frequencies and stats
        if(minProb_weight == 1){
          
          # calculate the frequency of agreement/disagreement permutations
          f =  rComparionsonExp(r_observed_exp, r_probable_exp, county)
          
          # save corresponding stats to countyList
          c("STATEFP" = county$STATEFP[1], "COUNTYFP" = county$COUNTYFP[1], calcKappa(f))
          
          # if multiple pixels have the minProb value and the nPixel break is such that they could be both included or excluded from "most probable"
        }else{
          
          # calculate confusion matrix for probable expansion excluding pixels where prob == minProb
          r_probable_exp_rmv_minProb = mask(r_probable_exp, r_tp_clip_noncropland == minProb, maskvalue = TRUE)
          
          f_rmv_minProb = rComparionsonExp(r_observed_exp, r_probable_exp_rmv_minProb, county)
          
          # calculate confusion matrix for probable expansion including pixels where prob == minProb -- to be weighted
          r_probable_exp_only_minProb = mask(r_probable_exp, r_tp_clip_noncropland != minProb, maskvalue = TRUE)
          r_probable_exp_only_minProb = mask(r_probable_exp_only_minProb, r_probable_exp_only_minProb, maskvalue = 0)
          
          f_only_minProb = rComparionsonExp(r_observed_exp, r_probable_exp_only_minProb, county)
          
          # apply the minProb weight (and it's inverse) to the corresponding comparisons of the minProb matrix
          f_weighted = data.frame(
            value = c(1, 11, 0, 10),
            count = c(
              ifelse(1 %in% f_only_minProb$value, f_only_minProb$count[f_only_minProb$value == 1], 0) * minProb_weight,
              ifelse(11 %in% f_only_minProb$value, f_only_minProb$count[f_only_minProb$value == 1], 0) * minProb_weight,
              ifelse(1 %in% f_only_minProb$value, f_only_minProb$count[f_only_minProb$value == 1], 0) * (1 - minProb_weight),
              ifelse(11 %in% f_only_minProb$value, f_only_minProb$count[f_only_minProb$value == 1], 0) * (1 - minProb_weight)),
            comparison = c("obsNO_predYES", "obsYES_predYES", "obsNO_predNO", "obsYES_predNO")
          )
          
          # join both comparison frequency data frames
          f_merged = merge(f_rmv_minProb, f_weighted, by = 'value')
          
          # add the minProb frequency values back to those calculated by removing minProv
          f = data.frame(
            value = f_merged$value,
            count = f_merged$count.x + f_merged$count.y,
            comparison = f_merged$comparison.x)
          
          # save corresponding stats to countyList
          c("STATEFP" = county$STATEFP[1], "COUNTYFP" = county$COUNTYFP[1], calcKappa(f))
          
        }
        
      }
      
    }, error = function(e) {
      message("Error for  ", tp_name, " in County ", county$COUNTYFP[1], ":")
      message(e)
      message("\n")
      return(NA)
      
    })
    
    cat(county$COUNTYFP[1], ',')
    
  }