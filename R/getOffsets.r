getNOffsetsMax = function(n, margin){
    as.integer(n*(1-margin)) - as.integer(n*margin) + 1L;
}

getOffsetsAll = function(n, margin){
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(margin) );
    stopifnot( is.numeric(n) );
    return( as.integer(n*margin):as.integer(n*(1-margin)) );
}

getOffsetsRandom = function(n, npermute, margin = 0.05){
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(margin) );
    stopifnot( is.numeric(n) );
    stopifnot( is.numeric(npermute) );
    offsets = sample.int(n = floor(n*(1-2*margin))+1, size = floor(npermute), replace = FALSE) + as.integer(n*margin-1);
    return( offsets );
}

getOffsetsUniform = function(n, npermute, margin = 0.05){
    stopifnot( margin >= 0 );
    stopifnot( margin < 0.5 );
    stopifnot( is.numeric(margin) );
    stopifnot( is.numeric(n) );
    stopifnot( is.numeric(npermute) );
    from = as.integer(n*margin)
    to = as.integer(n*(1-margin))
    stopifnot( npermute <= to - from );
    offsets = as.integer(seq.int(from = from, to = to, length.out = npermute));
    return( offsets );
}
