

initpoints <- function(boundaries, n.points, integer, method = "halton") {

    if (method == "halton") {
        s <- as.matrix(randtoolbox::halton(n = n.points -
            2, dim = length(boundaries)), ncol = length(boundaries))
        pmin <- pmax <- c()
        for (i in seq_len(length(boundaries))) {
            dmin <- boundaries[[i]][1]
            dmax <- boundaries[[i]][2]
            s[, i] <- dmin + s[, i] * (dmax - dmin)
            pmin[i] <- dmin
            pmax[i] <- dmax
        }
        if(integer) s <- round(s)
        points <- rbind(s, as.numeric(pmin), as.numeric(pmax))
    }

    if (method == "sobol") {
        s <- as.matrix(randtoolbox::sobol(n = n.points,
            dim = length(boundaries), scrambling = 0),
            ncol = length(boundaries))

        for (i in seq_len(length(boundaries))) {
            dmin <- boundaries[[i]][1]
            dmax <- boundaries[[i]][2]
            s[, i] <- dmin + s[, i] * (dmax - dmin)
        }
        if(integer) points <- round(s) else points <- s
    }


    if (method == "random") {

        s <- matrix(stats::runif(n.points * length(boundaries)),
            ncol = length(boundaries))

        for (i in seq_len(length(boundaries))) {
            dmin <- boundaries[[i]][1]
            dmax <- boundaries[[i]][2]
            s[, i] <- dmin + s[, i] * (dmax - dmin)
        }
        if(integer) points <- round(s) else points <- s
    }


    colnames(points) <- names(boundaries)

    return(points)
}


