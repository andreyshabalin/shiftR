#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <R_ext/Rdynload.h>

// __builtin_popcount
static inline int popcount( unsigned int i ){
    // bytesum = function(x){ a = rawToBits(as.raw(x)); return(sum(as.integer(a)));}
    // paste(sapply(0:255, bytesum), collapse = ',')
    static unsigned char bitsums[256] = {
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
    return(  
        bitsums[ i      & 0xFF] + 
        bitsums[(i>>8)  & 0xFF] +
        bitsums[(i>>16) & 0xFF] +
        bitsums[(i>>24)       ]);
}

int NumberOfSetBits(unsigned int i){
    // Java: use >>> instead of >>
    // C or C++: use uint32_t
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

SEXP CbitSum(SEXP x){
    x = PROTECT(coerceVector(x, INTSXP));
    int n = length(x);
    int *px = INTEGER(x);
    
    int sum = 0;
    for( int i=0; i<n; i++) {
        sum += NumberOfSetBits(px[i]);
    }
    SEXP rez = PROTECT(allocVector(INTSXP, 1)); 
    INTEGER(rez)[0] = sum;
    UNPROTECT(2);
    return(rez);
}

SEXP CbitSumAnd(SEXP x, SEXP y){
    x = PROTECT(coerceVector(x, INTSXP));
    int n1 = length(x);
    y = PROTECT(coerceVector(y, INTSXP));
    int n2 = length(x);
    int n = n1>n2 ? n2 : n1;
    int *px = INTEGER(x);
    int *py = INTEGER(y);
    
    int sum = 0;
    for( int i=0; i<n; i++) {
        sum += popcount(px[i] & py[i]);
    }
    SEXP rez = PROTECT(allocVector(INTSXP, 1)); 
    INTEGER(rez)[0] = sum;
    UNPROTECT(3);
    return(rez);
}

SEXP CbitSumOr(SEXP x, SEXP y){
    x = PROTECT(coerceVector(x, INTSXP));
    int n1 = length(x);
    y = PROTECT(coerceVector(y, INTSXP));
    int n2 = length(x);
    int n = n1>n2 ? n2 : n1;
    int *px = INTEGER(x);
    int *py = INTEGER(y);
    
    int sum = 0;
    for( int i=0; i<n; i++) {
        sum += popcount(px[i] | py[i]);
    }
    SEXP rez = PROTECT(allocVector(INTSXP, 1)); 
    INTEGER(rez)[0] = sum;
    UNPROTECT(3);
    return(rez);
}

SEXP CbitSumAndShifted(SEXP x, SEXP y, SEXP yoffset) {
    x = PROTECT(coerceVector(x, INTSXP));
    int n1 = length(x);
    y = PROTECT(coerceVector(y, INTSXP));
    int n2 = length(x);
    int n = n1>n2 ? n2 : n1;
    int *px = INTEGER(x);
    int *py = INTEGER(y);
    yoffset = PROTECT(coerceVector(yoffset, INTSXP));
    
    int offset = INTEGER(yoffset)[0];
    
    int sum = 0;
    for( int i=0; i<offset; i++) {
        sum += popcount(px[i-offset+n] & py[i]);
    }
    for( int i=offset; i<n; i++) {
        sum += popcount(px[i-offset  ] & py[i]);
    }
    
    SEXP rez = PROTECT(allocVector(INTSXP, 1)); 
    INTEGER(rez)[0] = sum;
    UNPROTECT(4);
    return(rez);
}

SEXP CbitSumAndYinX(SEXP x, SEXP y, SEXP xstart) {
    x = PROTECT(coerceVector(x, INTSXP));
    int n1 = length(x);
    y = PROTECT(coerceVector(y, INTSXP));
    int n = length(y);
    xstart = PROTECT(coerceVector(xstart, INTSXP));
    int start = INTEGER(xstart)[0];
    
    if( start < 0 ) return R_NilValue;
    if( n1 < start + n ) return R_NilValue;
    
    int *px = INTEGER(x);
    int *py = INTEGER(y);
    
    int sum = 0;
    for( int i=0; i<n; i++) {
    sum += popcount(px[i+start] & py[i]);
    }
    
    SEXP rez = PROTECT(allocVector(INTSXP, 1)); 
    INTEGER(rez)[0] = sum;
    UNPROTECT(4);
    return(rez);
}

static R_CallMethodDef callMethods[] = {
    {"CbitSum",					(DL_FUNC) &CbitSum, 1},
    {"CbitSumAnd",				(DL_FUNC) &CbitSumAnd, 2},
    {"CbitSumOr",				(DL_FUNC) &CbitSumOr, 2},
    {"CbitSumAndShifted",	(DL_FUNC) &CbitSumAndShifted, 3},
    {"CbitSumAndYinX",		(DL_FUNC) &CbitSumAndYinX, 3},
    {NULL, NULL, 0}
};

void R_init_shiftR(DllInfo *info)	{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
