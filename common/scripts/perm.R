upermn <- function(x)
    factorial(length(x)) / prod(factorial(as.numeric(table(x))))

uperm <- function(x) {
    require(combinat)
    u <- sort(unique(x))
    l <- length(u)

    if (l == length(x))
        do.call(rbind, permn(x))
    else if (l == 1)
        x
    else {
        result <- matrix(NA, upermn(x), length(x))
        index <- 1
        for (i in 1 : l) {
            v <- x[-which(x == u[i])[1]]
            newindex <- upermn(v)
            result[index : (index + newindex - 1), ] <-
                cbind(u[i], if (tabe(x)[i] == 1)
                      do.call(rbind, unique(permn(v)))
                else
                    uperm(v))

            index <- index + newindex
        }
        result
    }
}
