bitSum = function(x) {
	.Call("CbitSum", x, PACKAGE = "fastCircularPermutations")
}

bitSumAnd = function(x,y) {
	.Call("CbitSumAnd", x, y, PACKAGE = "fastCircularPermutations")
}

bitSumOr = function(x,y) {
	.Call("CbitSumOr", x, y, PACKAGE = "fastCircularPermutations")
}

bitSumAndShifted = function(x, y, yoffset) {
	.Call("CbitSumAndShifted", x, y, yoffset, PACKAGE = "fastCircularPermutations")
}

bitSumAndYinX = function(x, y, yoffset) {
	.Call("CbitSumAndYinX", x, y, yoffset, PACKAGE = "fastCircularPermutations")
}

prepareBinaryDataLeft = function(set) {
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

prepareBinaryDataRight = function(set) {
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

singleCircularPermutation = function(left, right, offset) {
	stopifnot( class(left) == "fcpLeft" );
	stopifnot( class(right) == "fcpRight" );
	stopifnot( attr(left, "len") == attr(right, "len") );
	stopifnot( offset >= 0L )
	stopifnot( offset < attr(right, "len") )
	.Call("CbitSumAndYinX", left, right[[1 + (offset) %% 32]], offset %/% 32, PACKAGE = "fastCircularPermutations")
}

circularPermutationAnalysis = function(left, right, npermute, margin = 0.05, alsoDoFisher = TRUE, returnPermOverlaps = FALSE) {
	rez = list();
	stopifnot( class(left) == "fcpLeft" );
	stopifnot( class(right) == "fcpRight" );
	stopifnot( attr(left, "len") == attr(right, "len") );
	stopifnot( npermute > 0 );
	sum1 = attr(left, "sum");
	sum2 = attr(left, "sum");
	sum12 = singleCircularPermutation(left, right, 0);
	len =  attr(left, "len");
	if(alsoDoFisher) {
		fisherMat = matrix(c(sum12,sum1-sum12,sum2-sum12,len-sum1-sum2+sum12),2,2);
		fisherTest = fisher.test(fisherMat);
		rez = list(fisherTest = fisherTest,
					  fisherMat = fisherMat)
	}

	overlapsPerm = integer(npermute);
	offsets = sample(len*(1-2*margin), size = npermute, replace = FALSE) + len*margin;
	for( i in seq_len(npermute) ) 
		overlapsPerm[i] = .Call("CbitSumAndYinX", left, right[[1 + (offsets[i]) %% 32]], offsets[i] %/% 32, PACKAGE = "fastCircularPermutations")
	permPVenrich  = max(mean(overlapsPerm >= sum12), 0.5/npermute);
	permPVdeplete = max(mean(overlapsPerm <= sum12), 0.5/npermute);
	permPV = min(permPVenrich, permPVdeplete, 0.5)*2;
	
	meanO = mean(overlapsPerm);
	stdO = sd(overlapsPerm);
	
	permZ = (sum12 - meanO) / stdO;

	rez = c(rez,list(
		nfeatures = len,
		lfeatures = sum1,
		rfeatures = sum2,
		overlap = sum12,
		overlapUnderNull = sum1 / len * sum2,
		enrichment = sum12 / sum1 / sum2 * len,
		permPVenrich = permPVenrich,
		permPVdeplete = permPVdeplete,
		permPV = permPV,
		permZ = permZ
		));
	if(returnPermOverlaps) {
		rez = c(list(overlapsPerm = overlapsPerm), rez);
	}
	return(rez);
}
