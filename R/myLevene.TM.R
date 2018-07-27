myLevene.TM = function (value, memSubjects, trim.alpha = 0.25) 
{
    u.memSubjects = sort(unique(memSubjects))
    if (length(u.memSubjects) != 2) {
        stop("memSubjects must take 2 and only 2 values\n")
    }
    if (!identical(u.memSubjects, c(0, 1))) {
        stop("memSubjects must only take values 0 or 1\n")
    }
    if (length(value) != length(memSubjects)) {
        stop("value must have the same length as memSubjects\n")
    }
    pos1 = which(memSubjects == 1)
    pos0 = which(memSubjects == 0)
    value1 = value[pos1]
    value0 = value[pos0]

    var1=var(value1, na.rm=TRUE)
    var0=var(value0, na.rm=TRUE)

    value1.m = mean(value1, na.rm = TRUE, trim = trim.alpha)
    value0.m = mean(value0, na.rm = TRUE, trim = trim.alpha)
    z1 = abs(value1 - value1.m)
    z0 = abs(value0 - value0.m)
    z = c(z1, z0)
    zbar = mean(z, na.rm = TRUE)
    z1.m = mean(z1, na.rm = TRUE)
    z0.m = mean(z0, na.rm = TRUE)
    n1 = length(pos1)
    n0 = length(pos0)
    numer = n1 * (z1.m - zbar)^2 + n0 * (z0.m - zbar)^2
    denom = sum((z1 - z1.m)^2, na.rm = TRUE)
    denom = denom + sum((z0 - z0.m)^2, na.rm = TRUE)
    N = length(value)
    W = (N - 2) * numer/denom
    pval = 1 - pf(W, df1 = 1, df2 = N - 2)
    #res = list(stat = W, pval = pval)
    if(var1 < var0)
    {
      W = -W
    }
    res = c(pval, W)
    names(res) = c("pvalue", "stat")
    return(res)
}

