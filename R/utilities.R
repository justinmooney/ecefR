#' Sets Earth Constants as global variables.
#'
#' @param a no idea
#' @param b no idea
#'
earthcon <- function(a, b) {

    f     <- 1 - b / a
    eccsq <- 1 - b ** 2 / a ** 2
    ecc   <- sqrt(eccsq)

    list(A = a, B = b, F = f, Ecc = ecc, Esq = eccsq)
}


#' WGS84 Earth Constants
wgs84 <- function() {
    wgs84a <- 6378.137
    wgs84f <- 1.0 / 298.257223563
    wgs84b <- wgs84a * (1.0 - wgs84f)

    earthcon(wgs84a, wgs84b)
}


#' compute the radii at the geodetic latitude lat (in degrees)
#'
#' @param lat geodetic latitude in degrees
#'
#' @return r, rn, rm in km
#'
radcur <- function(lat) {

    dtr <- pi / 180.0

    EARTH <- wgs84()

    a <- EARTH$A
    b <- EARTH$B

    asq   <- a ** 2
    bsq   <- b ** 2
    eccsq <- 1 - bsq / asq

    clat <- cos(dtr * lat)
    slat <- sin(dtr * lat)

    dsq <- 1.0 - eccsq * slat * slat
    d   <- sqrt(dsq)

    rn <- a / d
    rm <- rn * (1.0 - eccsq) / dsq

    rho <- rn * clat
    z   <- (1.0 - eccsq) * rn * slat
    rsq <- rho * rho + z * z
    r   <- sqrt(rsq)

    c(r, rn, rm)
}


#' physical radius of earth from geodetic latitude
#'
#' @param lat latitude at whic to calculate radius
#'
#' @return radius of earth at lat
#'
rearth <- function(lat) {
    radcur(lat)[1]
}


#' geocentric latitude to geodetic latitude
#'
#' @param flatgc geocentric latitude deg
#' @param altkm altitide in km
#'
#' @return geodetic latitude in deg
#'
gc2gd <- function(flatgc, altkm) {

    dtr <- pi / 180.0
    rtd <- 1 / dtr

    EARTH <- wgs84()

    ecc <- EARTH$Ecc
    esq <- ecc * ecc

    # approximation by stages
    # 1st use gc-lat as if is gd, then correct alt dependence

    altnow  <-  altkm

    rrnrm   <-  radcur(flatgc)
    rn      <-  rrnrm[2]

    ratio   <- 1 - esq * rn / (rn + altnow)

    tlat    <- tan(dtr * flatgc) / ratio
    flatgd  <- rtd * atan(tlat)

    #now use this approximation for gd-lat to get rn etc.

    rrnrm   <-  radcur(flatgd)
    rn      <-  rrnrm[2]

    ratio   <-  1  - esq * rn / (rn + altnow)
    tlat    <-  tan(dtr * flatgc) / ratio
    flatgd  <-  rtd * atan(tlat)

    flatgd

}


#' geodetic latitude to geocentric latitude
#'
#' @param flatgd geodetic latitude deg
#' @param altkm altitide in km
#'
#' @return geocentric latitude in deg
#'
gd2gc <- function(flatgd, altkm) {

    dtr <- pi / 180.0
    rtd <- 1 / dtr

    EARTH <- wgs84()

    ecc <- EARTH$Ecc
    esq <- ecc ** 2

    altnow <- altkm

    rrnrm <- radcur(flatgd)
    rn    <- rrnrm[2]

    ratio <- 1 - esq * rn / (rn + altnow)

    tlat   <- tan(dtr * flatgd) * ratio
    flatgc <- rtd * atan(tlat)

    flatgc

}


#' latitude longitude to east,north,up unit vectors
#'
#' @param flat latitude in degees N
#' @param flon longitude in degrees E
#'
#' @return enu[3[3]]  packed 3-unit vectors / each a 3 vector
#'
llenu <- function(flat, flon) {

    dtr <- pi / 180.0

    clat <- cos(dtr * flat)
    slat <- sin(dtr * flat)
    clon <- cos(dtr * flon)
    slon <- sin(dtr * flon)

    ee <- c(-slon, clon, 0.0)
    en <- c(-clon*slat, -slon * slat, clat)
    eu <- c(clon * clat, slon * clat, slat)

    c(ee, en, eu)
}



