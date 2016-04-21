#' A workaround using global namespaces
#'
#' @export
#'
geodGBL <- function() {
    # this should probably not be a thing
    if (!exists("EARTH_A"))
        wgs84()

}


#' Sets Earth Constants as global variables.
#'
#' @param ai no idea
#' @param bi no idea
#'
#' @export
#'
#'
earthcon <- function(ai, bi) {

    a <- as.numeric(ai)
    b <- as.numeric(bi)

    f <- 1 - b / a
    eccsq <- 1 - b * b / (a * a)
    ecc <- sqrt(eccsq)

    EARTH_A  <<-  a
    EARTH_B  <<-  b
    EARTH_F  <<- f
    EARTH_Ecc <<-  ecc
    EARTH_Esq <<- eccsq
}


#' WGS84 Earth Constants
#'
#' @export
#'
wgs84 <- function() {
    wgs84a <- 6378.137
    wgs84f <- 1.0 / 298.257223563
    wgs84b <- wgs84a * (1.0 - wgs84f)

    earthcon(wgs84a, wgs84b)
}


#' compute the radii at the geodetic latitude lat (in degrees)
#'
#' @param lati geodetic latitude in degrees
#'
#' @return r, rn, rm in km
#'
#' @export

radcur <- function(lati) {

    rrnrm = vector(length = 3)

    dtr   = pi / 180.0

    geodGBL()

    a     = EARTH_A
    b     = EARTH_B

    asq   = a * a
    bsq   = b * b
    eccsq  =  1 - bsq / asq
    ecc = sqrt(eccsq)

    lat   =  as.numeric(lati)

    clat  =  cos(dtr * lat)
    slat  =  sin(dtr * lat)

    dsq   =  1.0 - eccsq * slat * slat
    d     =  sqrt(dsq)

    rn    =  a / d
    rm    =  rn * (1.0 - eccsq) / dsq

    rho   =  rn * clat
    z     =  (1.0 - eccsq) * rn * slat
    rsq   =  rho * rho + z * z
    r     =  sqrt(rsq)

    rrnrm[1]  =  r
    rrnrm[2]  =  rn
    rrnrm[3]  =  rm

    rrnrm
}


#' physical radius of earth from geodetic latitude
#'
#' @param lati
#'
#' @return radius of earth at lat
#'
#' @export
#'
rearth <- function(lati) {

    lat <- as.numeric(lati)

    rrnrm <- radcur(lat)
    r     <-  rrnrm[1]

    r
}


#' geocentric latitude to geodetic latitude
#'
#' @param flatgci geocentric latitude deg
#' @param altkmi altitide in km
#'
#' @return geodetic latitude in deg
#'
#' @export
#'
gc2gd <- function(flatgci, altkmi) {

    dtr <- pi / 180.0
    rtd <- 1 / dtr

    rrnrm <- vector(length = 3)

    geodGBL()

    flatgc <- as.numeric(flatgci)
    altkm  <- as.numeric(altkmi)

    ecc <- EARTH_Ecc
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
#' @param flatgdi geodetic latitude deg
#' @param altkmi altitide in km
#'
#' @return geocentric latitude in deg
#'
#' @export
#'
gd2gc <- function(flatgdi, altkmi) {

    dtr <- pi / 180.0
    rtd <- 1 / dtr

    rrnrm <- vector(3)

    geodGBL()

    flatgd <- as.numeric(flatgdi)
    altkm  <- as.numeric(altkmi)

    ecc <- EARTH_Ecc
    esq <- ecc * ecc

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
#' @param flati latitude in degees N
#' @param floni longitude in degrees E
#'
#' @return enu[3[3]]  packed 3-unit vectors / each a 3 vector
#'
#' @export

llenu <- function(flati, floni) {

    ee <- vector(length = 3)
    en <- vector(length = 3)
    eu <- vector(length = 3)

    enu <- vector(length = 3)

    dtr <- pi / 180.0

    flat <- as.numeric(flati)
    flon <- as.numeric(floni)

    clat <- cos(dtr * flat)
    slat <- sin(dtr * flat)
    clon <- cos(dtr * flon)
    slon <- sin(dtr * flon)

    ee[1] <- -slon
    ee[2] <- clon
    ee[3] <- 0.0

    en[1] <- -clon * slat
    en[2] <- -slon * slat
    en[3] <- clat

    eu[1] <- clon * clat
    eu[2] <- slon * clat
    eu[3] <- slat

    enu[1] <- ee
    enu[2] <- en
    enu[3] <- eu

    enu
}



