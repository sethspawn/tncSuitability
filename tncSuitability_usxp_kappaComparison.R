rm(list = ls())

#====================================================================================================
# load required packages, installing any that have not yet been installed

packages = c(
  "boxr",
  "raster",
  "tigris",
  "sf",
  "fmsb",
  "ggplot2",
  "ggalt"
)

install.packages(setdiff(packages, rownames(installed.packages()))) 
lapply(packages, require, character.only = TRUE)

#====================================================================================================

# Authenticate BOX connection using Client ID and Secret in renviron

# Sys.setenv(BOX_CLIENT_ID = "")
# Sys.setenv(BOX_CLIENT_SECRET = "")

box_auth()

#====================================================================================================

rComparionsonExp = function(r_observed_exp, r_probable_exp, county){
  
  # create comparison raster identifying agreement and disagreement by class (see included key below)
  r_comparison_exp = (r_observed_exp + (r_probable_exp/10))*10
  r_comparison_exp = mask(r_comparison_exp, st_transform(county, st_crs(r_comparison_exp)))
  
  comparison_key = data.frame(value = c(0,1,10,11), comparison = c('obsNO_predNO', 'obsNO_predYES', 'obsYES_predNO', 'obsYES_predYES'))
  
  # calculate the frequency of each comparison outcome
  f = as.data.frame(freq(r_comparison_exp))
  f = merge(f, comparison_key, by = 'value')
  
  return(f)
  
}

calcKappa = function(f){
  
  # create output vector with comparison counts
  comp = f$count
  names(comp) = f$comparison
  
  kappa = tryCatch({
    
    # construct confusion matrix
    m = matrix(
      c(ifelse(length(f$count[f$value ==  0]) == 0, 0, f$count[f$value ==  0]),
        ifelse(length(f$count[f$value ==  1]) == 0, 0, f$count[f$value ==  1]),
        ifelse(length(f$count[f$value == 10]) == 0, 0, f$count[f$value == 10]),
        ifelse(length(f$count[f$value == 11]) == 0, 0, f$count[f$value == 11])),
      nrow = 2, ncol = 2, byrow = T)
    
    # calculate kappa statistic
    k = Kappa.test(m)
    k_est = k$Result$estimate
    k_pval = k$Result$p.value
    
    c('k_est' = k_est, 'k_pval' = k_pval)
    
  }, error = function(e) {
    
    return(NA)
    
  })
  
  return(c(comp, kappa))
  
}

#===================================================================================================

start.time =  Sys.time()

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

# download the usxp s35 mtr raster from BOX
box_dl(file_id = '752345336503', overwrite = T, file_name = 's35_mtr.tif')

# load the usxp s35 mtr raster (5 land use change classes)
r_mtr = raster("s35_mtr.tif")

# initialize a list in which to accumulate each state's county stats
state_county_stats = list()

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
  
  
  # get all possible unique column names in countyList
  all_names = unique(unlist(lapply(countyList, function(e) names(e))))
  
  # combine all list elements  by row into a data frame by matching column names
  df = do.call(rbind.data.frame, lapply(countyList, function(x) x[match(all_names, names(x))]))
  colnames(df) = all_names
  
  # convert all data columns of the resulting data frame to numeric
  df[which(!colnames(df) %in% c("STATEFP", "COUNTYFP"))] = lapply(df[which(!colnames(df) %in% c("STATEFP", "COUNTYFP"))], function(x) as.numeric(as.character(x)))
  
  # join county's statitics with it's geometry and native attributes
  counties = merge(counties, df, by = c("STATEFP", "COUNTYFP"))
  
  state_county_stats[[tp_name]] = st_as_sf(counties)
  
  cat(gsub("[+]", " ",gsub("_probs.tif", "", tp_name)), 'COMPLETE \n')
  
}


difftime(Sys.time(), start.time, units = "mins")


#===================================================================================================

# get all possible unique column names in countyList
all_names = lapply(state_county_stats, function(e) names(e))

state_county_stats = lapply(state_county_stats, function(df) df[,which(colnames(df) != 'Var.27')]) # <<<====== need to address


# combine all list elements  by row into a data frame by matching column names
all_county_stats = do.call(rbind.data.frame, lapply(state_county_stats, function(x) x[match(all_names, names(x))]))

colnames(df) = all_names
all_county_stats = st_as_sf(do.call(rbind.data.frame, state_county_stats))

#===================================================================================================

ggplot()+
  geom_sf(data = all_county_stats)+
  geom_sf(data = all_county_stats, aes(fill = k_est)) +
  geom_sf(data = st_union(all_county_stats[all_county_stats$k_pval <= 0.05,]), colour = 'red', fill = 'red', alpha = 0.1, size = 2)+
  coord_sf(crs = crs(r_mtr))+
  theme_bw()
