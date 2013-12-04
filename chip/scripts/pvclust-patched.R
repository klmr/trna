library(pvclust)

# In the following we patch the pvclust functions to use Spearman’s rank
# correlation in place of Pearson’s product–moment correlation.

dist.pvclust <- function (x, method = "euclidean", use.cor = "pairwise.complete.obs") {
    if (!is.na(pmatch(method, "correlation"))) {
        res <- as.dist(1 - cor(x, method = "spearman", use = use.cor))
        attr(res, "method") <- "correlation"
        return(res)
    }
    else if (!is.na(pmatch(method, "abscor"))) {
        res <- as.dist(1 - abs(cor(x, method = "spearman", use = use.cor)))
        attr(res, "method") <- "abscor"
        return(res)
    }
    else if (!is.na(pmatch(method, "uncentered"))) {
        if (sum(is.na(x)) > 0) {
            x <- na.omit(x)
            warning("Rows including NAs were omitted")
        }
        x <- as.matrix(x)
        P <- crossprod(x)
        qq <- matrix(diag(P), ncol = ncol(P))
        Q <- sqrt(crossprod(qq))
        res <- as.dist(1 - P/Q)
        attr(res, "method") <- "uncentered"
        return(res)
    }
    else dist(t(x), method)
}

pvclust <- function (data, method.hclust = "average", method.dist = "correlation", 
                     use.cor = "pairwise.complete.obs", nboot = 1000,
                     r = seq(0.5, 1.4, by = 0.1), store = FALSE, weight = FALSE) {
    n <- nrow(data)
    p <- ncol(data)
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- dist.pvclust(data, method = method.dist, use.cor = use.cor)
    data.hclust <- hclust(distance, method = method.hclust)
    size <- floor(n * r)
    rl <- length(size)
    if (rl == 1) {
        if (r != 1) 
            warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
        r <- list(1)
    }
    else r <- as.list(size/n)
    mboot <- lapply(r, pvclust:::boot.hclust, data = data, object.hclust = data.hclust, 
                    nboot = nboot, method.dist = method.dist, use.cor = use.cor, 
                    method.hclust = method.hclust, store = store, weight = weight)
    result <- pvclust:::pvclust.merge(data = data, object.hclust = data.hclust, 
                                      mboot = mboot)
    return(result)
}
