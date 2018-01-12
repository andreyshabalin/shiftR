cramerV = function(sum12, sum1, sum2, len){
    # Cramer's (V) coefficient
    # https://en.wikipedia.org/wiki/Phi_coefficient
    # http://www.aliquote.org/pub/3531_ftp.pdf
    len = as.numeric(len);
    sqrt(len / (sum1 * (len - sum1) * (sum2) * (len - sum2))) * 
        (sum12 * (len + (sum12- sum1 - sum2)) - 
            as.numeric(sum1-sum12) * (sum2-sum12));
}
