#' Convert LLA to ECEF
#'
#' @param lat latitude in degrees
#' @param lon longitude in degrees
#' @param alt altitude in km
#'
#' @return A 3-vec containing xyz for ECEF (km)
#'
#' @details
#' This converter uses N latitude and E longitude. Heights are in meters, and are ellipsoidal heights. There is no geoid model included.
#' The WGS 84 ellipsoid is used.
#'
#' @export
#'
lla2xyz <- function(lat, lon, alt) {

    flat  <- lat
    flon  <- lon
    altkm <- alt

    dtr <- pi / 180.0

    EARTH <- wgs84()

    clat <- cos(dtr * flat)
    slat <- sin(dtr * flat)
    clon <- cos(dtr * flon)
    slon <- sin(dtr * flon)

    rrnrm <- radcur(flat)
    rn    <- rrnrm[2]

    ecc <- EARTH$Ecc
    esq <- ecc ** 2

    x <- (rn + altkm) * clat * clon
    y <- (rn + altkm) * clat * slon
    z <- ((1 - esq) * rn + altkm) * slat

    c(x,y,z)
}

