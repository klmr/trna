require(combinat)

# Calculate unique permutations. Like `combinat::permn`, but disregarding
# duplicates. I.e. c(1, 1, 2) will generate permutations `c(1, 1, 2)`,
# `c(1, 2, 1)` and `c(2, 1, 1)` only.
upermnsize <- function(x)
    factorial(length(x)) / prod(factorial(as.numeric(table(x))))

upermn <- function(x) {
    u <- sort(unique(x))
    l <- length(u)

    if (l == length(x))
        do.call(rbind, permn(x))
    else if (l == 1)
        x
    else {
        result <- matrix(NA, upermnsize(x), length(x))
        index <- 1
        for (i in 1 : l) {
            v <- x[-which(x == u[i])[1]]
            newindex <- upermnsize(v)
            result[index : (index + newindex - 1), ] <-
                cbind(u[i], if (table(x)[i] == 1)
                      do.call(rbind, unique(permn(v)))
                else
                    uperm(v))

            index <- index + newindex
        }
        result
    }
}
