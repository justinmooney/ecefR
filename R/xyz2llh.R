#' Title
#'
#' @param xvec
#'
#' @export
#'
xyz2llh <- function(xvec) {
    # xyz vector  to  lat,lon,height
    #
    # input:
    #    xvec[3]   xyz ECEF location
    # output:
    #
    #    llhvec[3] with components
    #
    # flat      geodetic latitude in deg
    # flon      longitude in deg
    # altkm     altitude in km
    dtr <-  pi / 180.0
    rrnrm <- vector(length = 3)
    llhvec <- vector(length = 3)

    geodGBL()

    esq   <- EARTH_Esq

    x <- xvec[1]
    y <- xvec[2]
    z <- xvec[3]

    x <- as.numeric(x)
    y <- as.numeric(y)
    z <- as.numeric(z)

    rp <- sqrt(x * x + y * y + z * z)

    flatgc <- asin(z / rp) / dtr

    testval <- abs(x) + abs(y)

    if (testval < 1.0e-10) {
        flon <- 0.0
    } else {
        flon <- atan2(y, x) / dtr
    }
    if (flon < 0.0) {
        flon <- flon + 360.0
    }

    p <- sqrt(x * x + y * y)

    # on pole special case

    if (p < 1.0e-10) {
        flat <- 90.0
        if (z < 0.0) {
            flat <- -90.0
        }

        altkm <- rp - rearth(flat)
        llhvec[1]  <- flat
        llhvec[2]  <- flon
        llhvec[3]  <- altkm

        llhvec
    }

    # //        first iteration, use flatgc to get altitude
    # //        and alt needed to convert gc to gd lat.

    rnow  <- rearth(flatgc)
    altkm <- rp - rnow
    flat <- gc2gd(flatgc, altkm)

    rrnrm <- radcur(flat)
    rn    <- rrnrm[2]

    for (kount in 1:5)
    {
        slat  <- sin(dtr * flat)
        tangd <- (z + rn * esq * slat) / p
        flatn <- atan(tangd) / dtr

        dlat <- flatn - flat
        flat <- flatn
        clat <- cos(dtr * flat)

        rrnrm <- radcur(flat)
        rn    <- rrnrm[2]

        altkm <- (p / clat) - rn

        if (abs(dlat) < 1.0e-12) {
            break
        }

    }

    llhvec[1]  <- flat
    llhvec[2]  <- flon
    llhvec[3]  <- altkm

    llhvec

}
