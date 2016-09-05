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
