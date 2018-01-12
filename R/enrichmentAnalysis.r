# parallelPermFun = function(offsets, leftSets, rghtSets, obs_zscr){
parallelPermFun = function(offsets, lefttemp, rghttemp, obs_zscr){

    if( is.character(lefttemp) ) {
        leftSets = readRDS(file = lefttemp)
    } else {
        leftSets = lefttemp;
    }
    
    if( is.character(rghttemp) ) {
        rghtSets = readRDS(file = rghttemp)
    } else {
        rghtSets = rghttemp;
    }    
    
    # offsets = offsetList[[offsetsi]];
    # library(shiftR);
    obs_cntE = matrix(NA_real_, length(leftSets), length(rghtSets));
    obs_cntD = matrix(NA_real_, length(leftSets), length(rghtSets));
    obs_cntO = matrix(NA_real_, length(leftSets), length(rghtSets));
    all_zmin =  Inf;
    all_zmax = -Inf;
    # offsets = offsetList[[1]];
    for(i in seq_along(leftSets)){ # i=1
        for(j in seq_along(rghtSets)){ # j=1
            z = shiftrPermBinary(
                left = leftSets[[i]],
                right = rghtSets[[j]],
                offsets = offsets,
                alsoDoFisher = FALSE,
                returnPermOverlaps = TRUE);
            all_z = cramerV(z$overlapsPerm,
                            z$lfeatures,
                            z$rfeatures,
                            z$nfeatures);
            
            obs_cntE[i, j] = sum(all_z >= obs_zscr[i,j]);
            obs_cntD[i, j] = sum(all_z <= obs_zscr[i,j]);
            obs_cntO[i, j] = sum(abs(all_z) >= abs(obs_zscr[i,j]));
            all_zmin = pmin.int(all_zmin, all_z);
            all_zmax = pmax.int(all_zmax, all_z);
            rm(z, all_z);
        } # j in seq_along(rghtSets)
    } # i in seq_along(leftSets)
    obs_zmax = max(obs_zscr);
    obs_zmin = min(obs_zscr);
    all_cntE = sum(all_zmax >= obs_zmax);
    all_cntD = sum(all_zmin <= obs_zmin);
    all_cntO = sum(pmax(all_zmax, -all_zmin) >= max(obs_zmax, -obs_zmin));
    
    result = list(
        npermute = length(offsets),
        obs_cntE = obs_cntE,
        obs_cntD = obs_cntD,
        obs_cntO = obs_cntO,
        all_cntE = all_cntE,
        all_cntD = all_cntD,
        all_cntO = all_cntO);
    return(result);
}

enrichmentAnalysis = function(
            pvstats1,
            pvstats2,
            percentiles1 = NULL,
            percentiles2 = NULL,
            npermute,
            margin = 0.05,
            threads = 1){
    
    # extract values for easier work
    stopifnot( length(pvstats1) == length(pvstats2) );
    n = length(pvstats1);
    
    if(margin > 1)
        margin = margin / length(pvstats1);
    
    
    # Prepare offsets
    maxperm = getNOffsetsMax(n = n, margin = margin);
    
    if(npermute >= maxperm) {
        npermute = maxperm;
        offsets = getOffsetsAll(n = n, margin = margin);
    } else {
        offsets = getOffsetsRandom(n = n, npermute = npermute, margin = margin);
    }
    
    # library(parallel)
    # Initiate multithreading parameters
    if(is.null(threads))
        threads = TRUE;
    if(is.logical(threads)){
        if(threads){
            threads = detectCores();
        } else {
            threads = 1L;
        }
    }
    
    # split offsets by thread
    # offsetList = lapply(splitIndices(length(offsets), threads), 
    #                     function(x){offsets[x]});
    # rm(offsets);
    
    # get thresholds
    if(sum(!duplicated(pvstats1)) > 2){
        thresholds1 = quantile(
                            x = pvstats1,
                            probs = percentiles1,
                            na.rm = TRUE,
                            names = FALSE);
    } else {
        thresholds1 = min(pvstats1);
    }
    if(sum(!duplicated(pvstats2)) > 2){
        thresholds2 = quantile(
                            x = pvstats2,
                            probs = percentiles2,
                            na.rm = TRUE,
                            names = FALSE);
    } else {
        thresholds2 = min(pvstats2);
    }
    
    # Prepare left binary sets
    leftSets = vector('list', length(thresholds1));
    for(i in seq_along(thresholds1)){ # i=1
        bool1 = (pvstats1 <= thresholds1[i]);
        if((!any(bool1)) || all(bool1))
            stop("No variaton in primary data after mapping and thresholding");
        leftSets[[i]] = shiftrPrepareLeft(bool1);
        rm(bool1);
    }
    
    # Prepare right binary left sets
    rghtSets = vector('list', length(thresholds2));
    for(j in seq_along(thresholds2)){ # j=1
        bool2 = (pvstats2 <= thresholds2[j]);
        if((!any(bool2)) || all(bool2))
            stop("No variaton in enrichment data after mapping and thresholding");
        rghtSets[[j]] = shiftrPrepareRight(bool2);
        rm(bool2);
    }
    
    obs_zscr = matrix(NA_real_, length(thresholds1), length(thresholds2));
    for(i in seq_along(thresholds1)){ # i=1
        for(j in seq_along(thresholds2)){ # j=1
            z = shiftrPermBinary(
                left = leftSets[[i]],
                right = rghtSets[[j]], 
                offsets = c(),
                alsoDoFisher = FALSE,
                returnPermOverlaps = FALSE);
            obs_zscr[i, j] = cramerV(
                z$overlap,
                z$lfeatures,
                z$rfeatures,
                z$nfeatures);
            rm(z);
        }
    }

    if( threads > 1){
        
        lefttemp = tempfile();
        rghttemp = tempfile();
        saveRDS(file = lefttemp, leftSets, compress = FALSE)
        saveRDS(file = rghttemp, rghtSets, compress = FALSE)

        cl = makeCluster(threads);
        clres = clusterApplyLB(
            cl,
            clusterSplit(cl, offsets),
            parallelPermFun,
            # leftSets = leftSets,
            # rghtSets = rghtSets,
            lefttemp = lefttemp,
            rghttemp = rghttemp,
            obs_zscr = obs_zscr);
        stopCluster(cl);
        
        file.remove(lefttemp);
        file.remove(rghttemp);
        
        # combine the cluster results
        sumlist = lapply(names(clres[[1]]), 
                         function(nm){ Reduce('+',lapply(clres,`[[`,nm)) });
        names(sumlist) = names(clres[[1]]);
    } else {
        sumlist = parallelPermFun(
            offsets = offsets,
            lefttemp = leftSets,
            rghttemp = rghtSets,
            obs_zscr = obs_zscr);
    }
    stopifnot(npermute == sumlist$npermute)
    
    result = list(
        overallPV = c( 
            TwoSided   = sumlist$all_cntO / npermute,
            Enrichment = sumlist$all_cntE / npermute,
            Depletion  = sumlist$all_cntD / npermute),
        byThresholdPV = list(
            TwoSided   = sumlist$obs_cntO / npermute, 
            Enrichment = sumlist$obs_cntE / npermute,
            Depletion  = sumlist$obs_cntD / npermute)
    );
    return(result);
}
