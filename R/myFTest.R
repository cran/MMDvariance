myFTest=function (value, memSubjects) 
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
    res.L = stats::var.test(x = value[which(memSubjects == 0)], y = value[which(memSubjects == 
        1)])

    value1 = value[which(memSubjects == 1)]
    value0 = value[which(memSubjects == 0)]

    var1 = var(value1, na.rm=TRUE)
    var0 = var(value0, na.rm=TRUE)

    #res = list(stat = res.L$statistic, pval = res.L$p.value)
    pval=res.L$p.value
    stat = res.L$statistic
    if(var1 < var0)
    {
      stat = -stat
    }

    res=c(pval, stat)
    names(res)=c("pvalue", "stat")
    return(res)
}

