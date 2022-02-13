from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10
#from orbital_elements import *
from utils import *
from conversions import cartesian_to_spherical, spherical_to_cartesian, ecliptic_to_equatorial, elements_to_ecliptic
from corrections import moon_perts, jupiter_lon_perts, saturn_lon_perts, saturn_lat_perts, uranus_lon_perts


import numpy as np

sun_oe = {
    'N' : (20.552775873722602, 0.2136832880120676, -3.580978347278808, 174.97389301955806),
    'i' : (0.001101954782290859, 0.21369711394557486, -5.275796592097575, 0.00015489408592731231, 3.621325113811721e-07),
    'w' : (-20.630936449162093, 0.21368386463476155, -3.5858347681354386, 288.0401174521317),
    'a' : (0.0008991876561645148, 0.21276978132643024, 6.440250190274306, 1.000002237215457),
    'e' : (0.000624492722824836, 0.2299825556661901, 6.272769979941623, 0.01670036421395417),
    'M' : [ 9.85645883e-01, -6.84436431e+03]
    }

mercury_oe = {
    'N' : [2.5943137649076465e-33, -1.6458405468796637e-28, 4.1886116153189024e-24, -4.986022118098174e-20, 1.425664896833767e-16, 3.579902179195779e-12, -5.08222878798609e-08, 0.0003084130928174047, -0.9418197727279626, 1229.9985299396456],
    'i' : [-5.377291695654706e-36, 4.985297704367898e-31, -1.7171150986776112e-26, 2.7243783637920533e-22, -1.432580698949224e-18, -1.7317540599834052e-14, 3.4635949171992256e-10, -2.5036145236034457e-06, 0.00875602156565588, -5.3234684354500335],
    'w' : [7.266185018994313e-06, 29.130043647166488],
    'a' : (-7.302881780107657e-07, 0.015618668472226156, -0.16626523478656402, 0.3870983517945791),
    'e' : [-9.795222209711759e-79, 4.948344530013818e-74, -7.165646271347012e-70, -1.1487378401636948e-66, 6.912912897843176e-62, 4.067964777551289e-58, -5.014514038358126e-54, -8.163511441235564e-50, -8.754412696402857e-47, 8.607607039568674e-42, 8.339252306647572e-38, -3.0418308726218236e-34, -1.2193128129962863e-29, -5.338302719089814e-26, 1.1024320681178665e-21, 1.1731761292564323e-17, -9.765034347403027e-14, -1.4031494864886622e-09, 1.9017897264846455e-05, -0.08182346678242457, 125.15061900135669],
    'M' : [ 4.09233458e+00, -2.97113448e+04]
    }

class sun:
    """
    L   : Sun's mean longitude
    lon : Sun's true longitude
    """
    def __init__(self, d, obs_loc=None):
        self.name = 'sun'
        ecl = obl_ecl(d)
        self.d = d
        #self.N, self.i, self.w, self.a, self.e, self.M = sun_oe(d)

        A, w, p, c = sun_oe['N'] ; f_N = lambda t: A * np.sin(w*t + p) + c
        A, w, p, c, m = sun_oe['i'] ; f_i = lambda t: A * np.sin(w*t + p) + c + m*t
        A, w, p, c = sun_oe['w'] ; f_w = lambda t: A * np.sin(w*t + p) + c
        A, w, p, c = sun_oe['a'] ; f_a = lambda t: A * np.sin(w*t + p) + c
        A, w, p, c = sun_oe['e'] ; f_e = lambda t: A * np.sin(w*t + p) + c
        f_M = np.poly1d(sun_oe['M'])

        self.N = f_N(d)
        self.i = f_i(d)
        self.w = f_w(d)
        self.a = f_a(d)
        self.e = f_e(d)
        self.M = rev(f_M(d))

        
        self.L = rev(self.w + self.M)
        self.x_ecl, self.y_ecl, self.z_ecl, self.lon = elements_to_ecliptic('sun', self.N, self.i, self.w, self.a, self.e, self.M)
        self.x_equ, self.y_equ, self.z_equ = ecliptic_to_equatorial(self.x_ecl, self.y_ecl, self.z_ecl, d)
        self.ra, self.dec, self.r = cartesian_to_spherical(self.x_equ, self.y_equ, self.z_equ)

        if obs_loc is None:
            self.az, self.alt = None, None
        else:
            self.az, self.alt = self.altaz(obs_loc)
        

    def altaz(self, obs_loc):
        LON, LAT = obs_loc
        UT = getUT(self.d)
        
        GMST0 = rev(self.L + 180) / 15
        SIDTIME = GMST0 + UT + LON/15 # in hours
        RA_hour = self.ra/15 # convert RA from degrees to hours
        HA_hour = SIDTIME - RA_hour
        HA = HA_hour *15 # converts HA from hours to degrees

        x = cos(HA *rd) * cos(self.dec*rd)
        y = sin(HA *rd) * cos(self.dec*rd)
        z = sin(self.dec*rd)

        xhor = x * sin(LAT*rd) - z * cos(LAT*rd)
        yhor = y
        zhor = x * cos(LAT*rd) + z * sin(LAT*rd)

        az  = arctan2(yhor, xhor)*(180/pi) + 180
        alt = arcsin(zhor)*(180/pi)
        return az, alt
        

class moon:
    """
    Moon positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        N       : longitude of the ascending node
        i       : inclination to the ecliptic
        w       : argument of perihelion
        a       : semi-major axis, or mean distance from Sun
        e       : eccentricity
        M       : mean anomaly
        v       : true anomaly
        E       : eccentric anomaly
        L       : mean longitude
        ra      : Right Ascension (GCRS or topocentric if obs_loc is provided)
        dec     : Declination (GCRS or topocentric if obs_loc is provided)
        r       : distance to Earth
        ecl_lon : ecliptic longitude (GCRS)
        ecl_lat : ecliptic latitude (GCRS)
        x_ecl   : ecliptic x (GCRS)
        z_ecl   : ecliptic y (GCRS)
        y_ecl   : ecliptic z (GCRS)
        x_equ   : equatorial x (GCRS)
        y_equ   : equatorial y (GCRS)
        z_equ   : equatorial z (GCRS)
        elongation : elongation
        FV         : phase angle
    """
        
    def __init__(self, d, obs_loc=None):
        self.name = 'moon'
        ecl = obl_ecl(d)
        #self.obs_loc = obs_loc
        self._sun = sun(d)
        
        N, self.i, w, self.a, self.e, M = moon_oe(d)
        self.N, self.w, self.M = rev(N), rev(w), rev(M)
        
        x_ecl, y_ecl, z_ecl = elements_to_ecliptic('moon', self.N, self.i, self.w, self.a, self.e, self.M)
        ecl_lon, ecl_lat, ecl_r = cartesian_to_spherical(x_ecl, y_ecl, z_ecl)

        # CONSIDERING Perturbations
        # =========================

        self.L = rev(self.N + self.w + self.M) # Moon's mean longitude

        Ls = self._sun.L
        Ms = self._sun.M
        
        Lm = self.L
        Mm = self.M
        Nm = self.N

        D = Lm - Ls # Moon's mean elongation
        F = Lm - Nm # Moon's argument of latitude

        pert_lon, pert_lat, pert_r = moon_perts(Ls, Ms, Lm, Mm, Nm, D, F)

        # Add this to the ecliptic positions we earlier computed:
        self.ecl_lon = ecl_lon + pert_lon
        self.ecl_lat = ecl_lat + pert_lat
        self.r       = ecl_r   + pert_r

        r_ = 1

        x_ecl = r_ * cos(self.ecl_lon*rd) * cos(self.ecl_lat*rd)
        y_ecl = r_ * sin(self.ecl_lon*rd) * cos(self.ecl_lat*rd)
        z_ecl = r_ * sin(self.ecl_lat*rd)

        x_equ, y_equ, z_equ = ecliptic_to_equatorial(x_ecl, y_ecl, z_ecl, d)

        self.ra, self.dec, _ = cartesian_to_spherical(x_equ, y_equ, z_equ)

        # kh: in qesmato khodam anjam dadam (motmaen nistam)
        # badan check shavad
        self.x_ecl = x_ecl * self.r
        self.y_ecl = y_ecl * self.r
        self.z_ecl = z_ecl * self.r

        self.x_equ, self.y_equ, self.z_equ = ecliptic_to_equatorial(self.x_ecl, self.y_ecl, self.z_ecl, d)


        if obs_loc is not None:
            self.ra, self.dec = self.topocentric_correction(obs_loc)

        self.elongation = arccos( cos((self._sun.lon-self.ecl_lon)*rd) * cos(self.ecl_lat*rd) )*(180/pi)
        self.FV = 180 - self.elongation

    def topocentric_correction(self, obs_loc):
        LON, LAT = obs_loc
        mpar = arcsin(1/self.r)*(180/pi) # moon parallax
        gclat = LAT - 0.1924 * sin(2*LAT*rd)
        rho   = 0.99833 + 0.00167 * cos(2*LAT*rd)

        UT = getUT(d)

        GMST0 = rev(self._sun.L + 180) / 15
        LST = GMST0 + UT + LON/15 # in hours
        LST_deg = LST * 15
        HA = rev(LST_deg - self.ra)

        # auxiliary angle
        g = arctan( tan(gclat*rd) / cos(HA*rd) )*(180/pi)

        topRA  = self.ra  - mpar * rho * cos(gclat*rd) * sin(HA*rd) / cos(self.dec*rd)
        topDEC = self.dec - mpar * rho * sin(gclat*rd) * sin((g-self.dec)*rd) / sin(g*rd)
        return topRA, topDEC
        

class planet:
    """
    Planets positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        name (str) : name of the planet
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        
        lon_h_ecl  : Heliocentric ecliptic longitude
        lat_h_ecl  : Heliocentric ecliptic latitude
        r_h_ecl    : Heliocentric ecliptic distance
        xh_ecl     : Heliocentric ecliptic x
        yh_ecl     : Heliocentric ecliptic y
        zh_ecl     : Heliocentric ecliptic z
        lon_g_ecl  : Geocentric ecliptic longitude
        lat_g_ecl  : Geocentric ecliptic latitude
        r_g_ecl    : Geocentric ecliptic distance
        xg_ecl     : Geocentric ecliptic x
        yg_ecl     : Geocentric ecliptic y
        zg_ecl     : Geocentric ecliptic z
        ra         : Right Ascension (GCRS)
        dec        : Declination (GCRS)
        r          : Distance to Earth
        x          : equatorial x (GCRS)
        y          : equatorial y (GCRS)
        z          : equatorial z (GCRS)
        elongation : elongation
        FV         : phase angle
        mag        : Apparent magnitude
        diameter   : Apparent diameter
    """
    def __init__(self, name, d, obs_loc=None, epoch=None):
        self.name = name.lower()
        #self.obs_loc = obs_loc
        ecl = obl_ecl(d)
        self._sun = sun(d)
        
        if self.name=='mercury':
            N,i,w,a,e,M = mercury_oe(d)
        elif self.name=='venus':
            N,i,w,a,e,M = venus_oe(d)
        elif self.name=='mars':
            N,i,w,a,e,M = mars_oe(d)
        elif self.name=='jupiter':
            N,i,w,a,e,M = jupiter_oe(d)
        elif self.name=='saturn':
            N,i,w,a,e,M = saturn_oe(d)
        elif self.name=='uranus':
            N,i,w,a,e,M = uranus_oe(d)
        elif self.name=='neptune':
            N,i,w,a,e,M = neptune_oe(d)
        else:
            raise Exception('Planet name not valid!')

        #self.N, self.i, self.w, self.a, self.e, self.M = N,i,w,a,e,M

        #self.E = getE(e, M, dp=5)
        E = getE(e, M, dp=5)
        
        xv = a * (cos(E*rd) - e)
        yv = a * (sqrt(1 - e**2) * sin(E*rd))

        r = sqrt(xv**2 + yv**2) # in AU
        v = rev(arctan2(yv, xv) *(180/pi)) # true anomaly

        # rectangular heliocentric ecliptic
        xh_ecl = r * ( cos(N*rd) * cos((v+w)*rd) - sin(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        yh_ecl = r * ( sin(N*rd) * cos((v+w)*rd) + cos(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        zh_ecl = r * ( sin((v+w)*rd) * sin(i*rd) )

        # spherical heliocentric ecliptic
        self.lon_h_ecl, self.lat_h_ecl, self.r_h_ecl = \
                        cartesian_to_spherical(xh_ecl, yh_ecl, zh_ecl)

        # Correcting perturbations of Jupiter, Saturn and Uranus
        if self.name in ['jupiter', 'saturn', 'uranus']:
            Mj = jupiter_oe(d)[-1]
            Ms = saturn_oe(d)[-1]
            Mu = uranus_oe(d)[-1]
            if self.name=='jupiter':
                self.lon_h_ecl = self.lon_h_ecl + jupiter_lon_perts(Mj, Ms, Mu)                
            elif self.name=='saturn':
                self.lon_h_ecl = self.lon_h_ecl + saturn_lon_perts(Mj, Ms, Mu)
                self.lat_h_ecl = self.lat_h_ecl + saturn_lat_perts(Mj, Ms, Mu)
            elif self.name=='uranus':
                self.lon_h_ecl = self.lon_h_ecl + uranus_lon_perts(Mj, Ms, Mu)
                
        # Precession
        if epoch is not None:
            self.lon_h_ecl = self.lon_h_ecl + 3.82394E-5 * (365.2422 * (epoch-2000) - d)

        # kh: hypatie
        self.xh_ecl, self.yh_ecl, self.zh_ecl = \
                     spherical_to_cartesian(self.lon_h_ecl, self.lat_h_ecl, self.r_h_ecl)

        self.xg_ecl = self._sun.x_ecl + self.xh_ecl
        self.yg_ecl = self._sun.y_ecl + self.yh_ecl
        self.zg_ecl = self._sun.z_ecl + self.zh_ecl

        self.lon_g_ecl, self.lat_g_ecl, self.r_g_ecl = \
                        cartesian_to_spherical(self.xg_ecl, self.yg_ecl, self.zg_ecl)

        # equatorial rectangular geocentric
        self.x, self.y, self.z = ecliptic_to_equatorial(self.xg_ecl, self.yg_ecl, self.zg_ecl, d)

        # converet to spherical coordinates (radec)
        self.ra, self.dec, self.r = cartesian_to_spherical(self.x, self.y, self.z)

        # Phase angle and the elongation
        R = self.r
        r = self.r_h_ecl
        s = self._sun.r

        self.elongation = arccos((s**2 + R**2 - r**2)/(2*s*R))*(180/pi)
        FV    = arccos((r**2 + R**2 - s**2)/(2*r*R))*(180/pi)
        self.FV = FV
        #self.phase =  (1 + cos(self.FV*rd))/2

        # Magnitude
        if self.name=='mercury':
            d0 = 6.74
            mag = -0.36 + 5*log10(r*R) + 0.027 * FV + 2.2E-13 * FV**6
        elif self.name=='venus':
            d0 = 16.92
            mag = -4.34 + 5*log10(r*R) + 0.013 * FV + 4.2E-7  * FV**3
        elif self.name=='mars':
            d0 = 9.32
            mag = -1.51 + 5*log10(r*R) + 0.016 * FV
        elif self.name=='jupiter':
            d0 = 191.01
            mag = -9.25 + 5*log10(r*R) + 0.014 * FV
        elif self.name=='saturn':
            d0 = 158.2
            ir = 28.06 # tilt rings to ecliptic
            Nr = 169.51 + 3.82E-5 * d # ascending node of plane of rings
            los = self.lon_g_ecl # Saturn's geocentric ecliptic longitude
            las = self.lat_g_ecl # Saturn's geocentric ecliptic latitude
            # B : tilt of Saturn's rings
            B = arcsin(sin(las*rd) * cos(ir*rd) - cos(las*rd) * sin(ir*rd) * sin((los-Nr)*rd))*(180/pi)
            ring_magn = -2.6 * sin(abs(B)*rd) + 1.2 * (sin(B*rd))**2
            mag = -9.0  + 5*log10(r*R) + 0.044 * FV + ring_magn
        elif self.name=='uranus':
            d0 = 63.95
            mag = -7.15 + 5*log10(r*R) + 0.001 * FV
        elif self.name=='neptune':
            d0 = 61.55
            mag = -6.90 + 5*log10(r*R) + 0.001 * FV
        else:
            mag = None
        self.mag = mag
        self.diameter = d0 / self.r

        

#d = day(1990, 4, 19, 0)
d = day(2022, 2, 11, 10)

s = sun(d)

from astropy.coordinates import get_sun
from astropy.time import Time

t = Time('2022-02-11 10:00:00')
c = get_sun(t)
