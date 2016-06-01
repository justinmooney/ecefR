#' ECEF Coordinates to Lat, Long, Altitude
#'
#' @param x x-coord
#' @param y y-coord
#' @param z z-coord
#'
#' @return A 3-vec containing lat(deg), long(deg), altitude(km)
#'
#' @details
#' This converter uses N latitude and E longitude.
#' Heights are in meters, and are ellipsoidal heights.
#' There is no geoid model included.
#' The WGS 84 ellipsoid is used.
#'
#' @export
#'
xyz2lla <- function(x, y, z) {
    dtr <-  pi / 180.0

    EARTH <- wgs84()

    esq <- EARTH$Esq

    rp <- sqrt(x ^ 2 + y ^ 2 + z ^ 2)

    flatgc <- asin(z / rp) / dtr

    testval <- abs(x) + abs(y)
    if (testval < 1.0e-10)
        flon <- 0.0
    else
        flon <- atan2(y, x) / dtr

    if (flon < 0.0)
        flon <- flon + 360.0

    p <- sqrt(x ^ 2 + y ^ 2)

    # on pole special case
    if (p < 1.0e-10) {
        flat <- 90.0

        if (z < 0.0)
            flat <- -90.0

        altkm <- rp - rearth(flat)

        return(c(flat, flon, altkm))  # i dont even know how to test this
    }

    # first iteration, use flatgc to get altitude
    # and alt needed to convert gc to gd lat.

    rnow  <- rearth(flatgc)
    altkm <- rp - rnow
    flat  <- gc2gd(flatgc, altkm)

    rrnrm <- radcur(flat)
    rn    <- rrnrm[2]

    #for (i in 1:5) {  # why only five times?
    while (TRUE) {
        slat  <- sin(dtr * flat)
        tangd <- (z + rn * esq * slat) / p
        flatn <- atan(tangd) / dtr

        dlat <- flatn - flat
        flat <- flatn
        clat <- cos(dtr * flat)

        rrnrm <- radcur(flat)
        rn    <- rrnrm[2]

        altkm <- (p / clat) - rn

        if (abs(dlat) < 1.0e-12)
            break
    }

    return(c(flat, flon, altkm))
}
