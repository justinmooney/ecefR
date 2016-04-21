

# geodGBL -----------------------------------------------------------------

#' Title
#'
#' @export
#'
geodGBL <- function() {

   if (!exists("EARTH_A"))
      wgs84()

}


# earthcon ----------------------------------------------------------------

#' Title
#'
#' @param ai
#' @param bi
#'
#' @export
#'
earthcon <- function(ai, bi) {
   # Sets Earth Constants as globals
   # --  input a,b
   # --  Leaves Globals
   # EARTH_A      EARTH_B   EARTH_F  EARTH_Ecc
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



# wgs84 -------------------------------------------------------------------

#' Title
#'
#' @export
#'
wgs84 <- function() {
   # WGS84 Earth Constants
   # --  returns a,b,f,e  --
   #    --  Leaves Globals
   # EARTH_A      EARTH_B   EARTH_F  EARTH_Ecc
   #
   #


   wgs84a <- 6378.137
   wgs84f <- 1.0 / 298.257223563
   wgs84b <- wgs84a * (1.0 - wgs84f)

   earthcon(wgs84a, wgs84b)

}


# radcur ------------------------------------------------------------------

#' Title
#'
#' @param lati
#'
#' @export

radcur <- function(lati) {
   #    compute the radii at the geodetic latitude lat (in degrees)
   #
   # input:
   #    lat       geodetic latitude in degrees
   # output:
   #    rrnrm     an array 3 long
   # r,  rn,  rm   in km

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


# rearth ------------------------------------------------------------------


#' Title
#'
#' @param lati
#'
#' @export
#'
rearth <- function(lati) {
   # physical radius of earth from geodetic latitude

   lat <- as.numeric(lati)

   rrnrm <- radcur(lat)
   r     <-  rrnrm[1]

   r

}



# gc2gd -------------------------------------------------------------------

#' Title
#'
#' @param flatgci
#' @param altkmi
#'
#' @export
#'
gc2gd <- function(flatgci, altkmi) {
   # geocentric latitude to geodetic latitude

   # Input:
   #    flatgc    geocentric latitude deg.
   # altkm     altitide in km
   # ouput:
   #    flatgd    geodetic latitude in deg
   dtr   = pi / 180.0
   rtd   = 1 / dtr

   rrnrm <- vector(length = 3)

   geodGBL()

   flatgc <- as.numeric(flatgci)
   altkm <- as.numeric(altkmi)

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

# gd2gc -------------------------------------------------------------------

#' Title
#'
#' @param flatgdi
#' @param altkmi
#'
#' @export
#'
gd2gc <- function(flatgdi, altkmi) {
   # geodetic latitude to geocentric latitude
   #
   # Input:
   #    flatgd    geodetic latitude deg.
   # altkm     altitide in km
   # ouput:
   #    flatgc    geocentric latitude in deg

   dtr <- pi / 180.0
   rtd <- 1 / dtr

   rrnrm <- vector(3)

   geodGBL()

   flatgd <- as.numeric(flatgdi)
   altkm <- as.numeric(altkmi)

   ecc   <-  EARTH_Ecc
   esq   <-  ecc * ecc

   altnow  <-  altkm

   rrnrm   <-  radcur(flatgd)
   rn      <-  rrnrm[2]

   ratio   <- 1 - esq * rn / (rn + altnow)

   tlat    <- tan(dtr * flatgd) * ratio
   flatgc  = rtd * atan(tlat)

   flatgc

}

# llenu -------------------------------------------------------------------

#' Title
#'
#' @param flati
#' @param floni
#'
#' @export

llenu <- function(flati, floni) {
   # latitude longitude to east,north,up unit vectors
   #
   # input:
   #    flat      latitude in degees N
   # [ gc -> gc enu,  gd usual enu ]
   # flon      longitude in degrees E
   # output:
   #    enu[3[3]]  packed 3-unit vectors / each a 3 vector

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

   ee[1]  <-  -slon
   ee[2]  <-   clon
   ee[3]  <-   0.0

   en[1]  <-  -clon * slat
   en[2]  <-  -slon * slat
   en[3]  <-          clat

   eu[1]  <-   clon * clat
   eu[2]  <-   slon * clat
   eu[3]  <-          slat

   enu[1] <-  ee
   enu[2] <-  en
   enu[3] <-  eu

   enu

}



