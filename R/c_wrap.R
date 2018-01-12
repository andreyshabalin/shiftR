bitSum = function(x){
    .Call("CbitSum", x, PACKAGE = "shiftR");
}

bitSumAnd = function(x,y){
    .Call("CbitSumAnd", x, y, PACKAGE = "shiftR");
}

bitSumOr = function(x,y){
    .Call("CbitSumOr", x, y, PACKAGE = "shiftR");
}

bitSumAndShifted = function(x, y, yoffset){
    .Call("CbitSumAndShifted", x, y, yoffset, PACKAGE = "shiftR");
}

bitSumAndYinX = function(x, y, yoffset){
    .Call("CbitSumAndYinX", x, y, yoffset, PACKAGE = "shiftR");
}
