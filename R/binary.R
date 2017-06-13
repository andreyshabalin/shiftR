shiftrPrepareLeft = function(set){
    if(typeof(set) == "double")
    	set = as.integer(set);
    zero = switch(typeof(set),
        "integer" = 0L,
        "logical" = FALSE,
        "raw" = as.raw(0));
    n = length(set);
    rez = packBits( c(set, set, rep(zero, (-2*n)%%32)), type = 'integer');
    class(rez) = "fcpLeft";
    attr(rez, "sum") = sum(set);
    attr(rez, "len") = length(set);
    return(rez);
}

shiftrPrepareRight = function(set){
	if(typeof(set) == "double")
		set = as.integer(set);
	zero = switch(typeof(set),
					  "integer" = 0L,
					  "logical" = FALSE,
					  "raw" = as.raw(0));
	n = length(set);
	rez = vector('list', 32);
	for( i in 0:31 ) {
		# cat('i',i,'\n')
		rez[[i+1]] = packBits( c(rep(zero, i), set, rep(zero, (-i-n)%%32)), type = 'integer');
	}
	class(rez) = "fcpRight";
	attr(rez, "sum") = sum(set);
	attr(rez, "len") = length(set);
	return(rez);
}

singlePermutation = function(left, right, offset){
	stopifnot( class(left) == "fcpLeft" );
	stopifnot( class(right) == "fcpRight" );
	stopifnot( attr(left, "len") == attr(right, "len") );
	stopifnot( offset >= 0L )
	stopifnot( offset < attr(right, "len") )
	.Call("CbitSumAndYinX", left, right[[1 + (offset) %% 32]], offset %/% 32, PACKAGE = "shiftR")
}

getNOffestsMax = function(n, margin){
    as.integer(n*(1-margin)) - as.integer(n*margin) + 1L;
}

# n = 100; margin = 0.05; npermute = 13
getOffsetsAll = function(n, margin){
	stopifnot( margin >= 0 );
	stopifnot( margin < 0.5 );
	stopifnot( is.numeric(margin) );
	stopifnot( is.numeric(n) );
	return( as.integer(n*margin):as.integer(n*(1-margin)) );
}
# typeof( shiftR:::getOffsetsAll(n,margin))
getOffsetsRandom = function(n, npermute, margin = 0.05){
	stopifnot( margin >= 0 );
	stopifnot( margin < 0.5 );
	stopifnot( is.numeric(margin) );
	stopifnot( is.numeric(n) );
	stopifnot( is.numeric(npermute) );
	offsets = sample.int(n = floor(n*(1-2*margin))+1, size = floor(npermute), replace = FALSE) + as.integer(n*margin-1);
	return( offsets );
}
# shiftR:::getOffsetsRandom(n=100, npermute=91, margin=0.05)
# range( shiftR:::getOffsetsRandom(n=100, npermute=91, margin=0.05) )
# typeof( shiftR:::getOffsetsRandom(n=100, npermute=91, margin=0.05))

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
# shiftR:::getOffsetsUniform(n=100, npermute=91, margin=0.05)
# range( shiftR:::getOffsetsUniform(n=100, npermute=91, margin=0.05) )
# typeof( shiftR:::getOffsetsUniform(n=100, npermute=91, margin=0.05))

shiftrPermBinary = function(left, right, offsets, alsoDoFisher = TRUE, returnPermOverlaps = FALSE){
	rez = list();
	stopifnot( typeof(offsets) %in% c("NULL","integer") );
	stopifnot( class(left) == "fcpLeft" );
	stopifnot( class(right) == "fcpRight" );
	stopifnot( attr(left, "len") == attr(right, "len") );
	sum1 = attr(left, "sum");
	sum2 = attr(right, "sum");
	sum12 = singlePermutation(left, right, 0);
	len =  attr(left, "len");
	if(alsoDoFisher) {
		fisherMat = matrix(c(sum12,sum1-sum12,sum2-sum12,len-sum1-sum2+sum12),2,2);
		fisherTest = fisher.test(fisherMat);
		rez = list(fisherTest = fisherTest,
					  fisherMat = fisherMat)
	}

	if(length(offsets)>0L){
    	overlapsPerm = integer(length(offsets));
    	for( i in seq_along(offsets) ) 
    		overlapsPerm[i] = .Call("CbitSumAndYinX", left, right[[1 + (offsets[i]) %% 32]], offsets[i] %/% 32, PACKAGE = "shiftR")
    	permPVenrich  = max(mean(overlapsPerm >= sum12), 0.5/length(offsets));
    	permPVdeplete = max(mean(overlapsPerm <= sum12), 0.5/length(offsets));
    	permPV = min(permPVenrich, permPVdeplete, 0.5)*2;
    	
    	meanO = mean(overlapsPerm);
    	stdO = sd(overlapsPerm);
    	
    	permZ = (sum12 - meanO) / stdO;
    	
    	rezPerm = list(
    	    permPVenrich = permPVenrich,
    	    permPVdeplete = permPVdeplete,
    	    permPV = permPV,
    	    permZ = permZ);
	} else {
	    rezPerm = NULL;
	}

	rez = c(rez, 
	        list(
        		nfeatures = len,
        		lfeatures = sum1,
        		rfeatures = sum2,
        		overlap = sum12,
        		overlapUnderNull = sum1 / len * sum2,
        		oddsRatio = sum12 * as.numeric(len - sum1 - sum2 + sum12) /
        		            (as.numeric(sum1-sum12) * (sum2-sum12))),
    		rezPerm);
	if(returnPermOverlaps){
		rez = c(list(overlapsPerm = overlapsPerm), rez);
	}
	return(rez);
}
