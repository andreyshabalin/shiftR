minOffset = function(n, margin){
    as.integer(ceiling(n*margin));
}

maxOffset = function(n, margin){
    as.integer(floor(n*(1-margin)));
}

getNOffsetsMax = function(n, margin){
    maxOffset(n, margin) - minOffset(n, margin) + 1L;
}

getOffsetsAll = function(n, margin){
    stopifnot( is.numeric(margin) );
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(n) );
    return( ceiling(n*margin):as.integer(n*(1-margin)) );
}

getOffsetsRandom = function(n, npermute, margin = 0.05){
    stopifnot( is.numeric(margin) );
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(n) );
    stopifnot( is.numeric(npermute) );
    offsets = sample.int(n = getNOffsetsMax(n, margin), size = floor(npermute), replace = FALSE) + (minOffset(n, margin)-1L);
    return( offsets );
}

getOffsetsUniform = function(n, npermute, margin = 0.05){
    stopifnot( is.numeric(margin) );
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(n) );
    stopifnot( is.numeric(npermute) );
    fr = minOffset(n, margin)
    to = maxOffset(n, margin)
    stopifnot( npermute <= to - fr + 1);
    offsets = as.integer(seq.int(from = fr, to = to, length.out = npermute));
    return( offsets );
}
