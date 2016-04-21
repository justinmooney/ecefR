#' Title
#'
#' @param vec
#'
#' @details
#' This converter uses N latitude and E longitude. Heights are in meters, and are ellipsoidal heights. There is no geoid model included.
#' The WGS 84 ellipsoid is used.
#' @export
#'
llh2xyz <- function(vec) {
    # lat,lon,height to xyz vector
    #
    # input:
    #    flat      geodetic latitude in deg
    # flon      longitude in deg
    # altkm     altitude in km
    # output:
    #    returns vector x 3 long ECEF in km
    flati <- vec[1]
    floni <- vec[2]
    altkmi <- vec[3]

    dtr <- pi / 180.0
    rrnrm <- vector(length = 3)
    xvec <- vector(length = 3)

    geodGBL()

    flat  <- as.numeric(flati)
    flon  <- as.numeric(floni)
    altkm  <- as.numeric(altkmi)

    clat <- cos(dtr * flat)
    slat <- sin(dtr * flat)
    clon <- cos(dtr * flon)
    slon <- sin(dtr * flon)

    rrnrm  <- radcur(flat)
    rn     <- rrnrm[2]
    re     <- rrnrm[1]

    ecc    <- EARTH_Ecc
    esq    <- ecc * ecc

    x      <-  (rn + altkm) * clat * clon
    y      <-  (rn + altkm) * clat * slon
    z      <-  ((1 - esq) * rn + altkm) * slat

    xvec[1]  <- x
    xvec[2]  <- y
    xvec[3]  <- z

    xvec

}

