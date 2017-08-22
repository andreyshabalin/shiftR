simulateAR = function(n, cor){
    rn = rnorm(n);
    mult = sqrt(1 - cor^2);
    for( i in 2:n ){
        rn[i] = rn[i-1] * cor + rn[i] * mult;
    }
    return(rn);
}

simulateNumeric = function(n, corWithin, corAcross = 0){
    stopifnot(n>1)
    
    data1 = simulateAR(n, corWithin);
    data2 = simulateAR(n, corWithin);
  
    if(corAcross != 0){
        data2 = data1 * corAcross + data2 * sqrt(1 - corAcross^2);
    }
    return( list( data1 = data1, data2 = data2) );
}

simulateBinary = function(n, corWithin, corAcross = 0){
    tmp = simulateNumeric(n, corWithin, corAcross)
    return( list( data1 = as.integer(tmp$data1>0), data2 = as.integer(tmp$data2>0)) );
}

simulatePValues = function(n, corWithin, corAcross = 0){
    tmp = simulateNumeric(n, corWithin, corAcross)
    return( list( data1 = pnorm(tmp$data1), data2 = pnorm(tmp$data2)) );
}
