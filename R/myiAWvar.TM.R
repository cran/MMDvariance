myiAWvar.TM=function (value, memSubjects, trim.alpha = 0.25) 
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

    var1 = var(value1, na.rm=TRUE)
    var0 = var(value0, na.rm=TRUE)

    m.value1 = mean(value1, na.rm = TRUE, trim = trim.alpha)
    m.value0 = mean(value0, na.rm = TRUE, trim = trim.alpha)
    value1.2 = abs(value1 - m.value1)
    value0.2 = abs(value0 - m.value0)
    z = rep(NA, length(value))
    z[pos1] = value1.2
    z[pos0] = value0.2
    ybar = mean(memSubjects, na.rm = TRUE)
    U2 = sum((memSubjects - ybar) * z, na.rm = TRUE)
    zbar = mean(z, na.rm = TRUE)
    varU2 = ybar * (1 - ybar) * sum((z - zbar)^2, na.rm = TRUE)
    T2 = U2^2/varU2
    pval = 1 - pchisq(T2, df = 1)
    #res = list(U2 = U2, varU2 = varU2, stat = T2, pval = pval, 
    #    z = z, zbar = zbar)
    if(var1 < var0)
    {
      T2 = -T2
    }
    res = c(pval, T2)
    names(res) = c("pvalue", "stat")
    return(res)
}

