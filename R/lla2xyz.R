#' Convert LLA to ECEF
#'
#' @param vec lat)deg), long(deg), and altitude(km)
#'
#' @return A 3-vec containing xyz for ECEF (km)
#'
#' @details
#' This converter uses N latitude and E longitude. Heights are in meters, and are ellipsoidal heights. There is no geoid model included.
#' The WGS 84 ellipsoid is used.
#'
#' @export
#'
lla2xyz <- function(vec) {

    flat  <- as.numeric(vec[1])
    flon  <- as.numeric(vec[2])
    altkm <- as.numeric(vec[3])

    dtr <- pi / 180.0

    geodGBL()

    clat <- cos(dtr * flat)
    slat <- sin(dtr * flat)
    clon <- cos(dtr * flon)
    slon <- sin(dtr * flon)

    rrnrm <- radcur(flat)
    rn    <- rrnrm[2]

    ecc <- EARTH_Ecc
    esq <- ecc ** 2

    x <- (rn + altkm) * clat * clon
    y <- (rn + altkm) * clat * slon
    z <- ((1 - esq) * rn + altkm) * slat

    xvec <- c(x,y,z)

    xvec
}

