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
	n = length(set);
	rez = packBits( c(set, set, rep(zero, (-n)%%32)), type = 'integer');
	class(rez) = "fcpLeft";
	return(rez);
}

prepareBinaryDataRight = function(set0) {
	rez = vector('list', 32);
	for( i in 0:31 ) {
		# cat('i',i,'\n')
		bin2[[i+1]] = packBits( c(rep(zero, i), set2, rep(zero, padding-i)), type = 'integer');
	}
	class(rez) = "fcpRight";
	return(rez);
}
	
