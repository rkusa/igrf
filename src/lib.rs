#![allow(unused)]

// Sources:
// - https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
// - https://github.com/nsdecicco/libgeomag

use std::borrow::Cow;
use std::cmp::Ordering;

#[derive(Debug, PartialEq)]
pub struct Field {
    /// Declination of the field from the geographic north (deg).
    pub d: f64,
    /// Inclination of the field (deg).
    pub i: f64,
    pub h: f64,
    /// Total field intensity.
    pub f: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Annual rate of change of declination. (arc-min/yr).
    pub ddot: f64,
    pub fdot: f64,
    pub hdot: f64,
    /// Annual rate of change of inclination. (arc-min/yr).
    pub idot: f64,
    pub xdot: f64,
    pub ydot: f64,
    pub zdot: f64,
}

#[derive(Debug)]
struct FieldModel {
    name: &'static str,
    /// Epoch of model.
    epoch: f64,
    /// Min year of model.
    yrmin: f64,
    /// Max year of model.
    yrmax: f64,
    /// Minimum height of each model.
    altmin: f64,
    /// Maximum height of each model.
    altmax: f64,
    /// Main field coefficient.
    max1: usize,
    /// Secular variation coefficient.
    max2: usize,
    /// Acceleration coefficient.
    max3: usize,
    /// Schmidt quasi-normal internal spherical harmonic coeff.
    gh1: &'static [f64],
    /// Schmidt quasi-normal internal spherical harmonic coeff.
    gh2: &'static [f64],
}

include!(concat!(env!("OUT_DIR"), "/data.rs"));

pub fn declination(lat: f64, lon: f64, alt: u32, date: time::Date) -> Result<Field, Error> {
    if !(-90.0..=90.0).contains(&lat) {
        return Err(Error::LatitudeOutOfBounds);
    }

    if !(-360.0..=360.0).contains(&lon) {
        return Err(Error::LongitudeOutOfBounds);
    }

    let alt = f64::from(alt) / 1_000.0; // form meters to km
    let date = decimal_day_of_year(date);

    let min_year = MODELS[0].yrmin;
    let max_year = MODELS[0].yrmin;
    let is_date_in_range = (min_year..max_year).contains(&date);

    let i = MODELS
        .iter()
        .position(|m| (m.yrmin..m.yrmax).contains(&date))
        .unwrap_or_else(|| if date < min_year { 0 } else { MODELS.len() - 1 });

    let model = &MODELS[i];
    let next_model = MODELS.get(i + 1);

    eprintln!("{:#?}", model);

    let mut gha = [0f64; 13 * (13 + 2) + 1];
    let mut ghb = [0f64; 13 * (13 + 2) + 1];

    let nmax = if let Some(next_model) = next_model {
        interpsh(
            date,
            model.yrmin,
            model.max1,
            next_model.yrmin,
            next_model.max1,
            model.gh1,
            model.gh2,
            &mut gha,
        );
        interpsh(
            date + 1.0,
            model.yrmin,
            model.max1,
            next_model.yrmin,
            next_model.max1,
            model.gh1,
            model.gh2,
            &mut ghb,
        )
    } else {
        extrapsh(
            date,
            model.epoch,
            model.max1,
            model.max2,
            model.gh1,
            model.gh2,
            &mut gha,
        );
        extrapsh(
            date + 1.0,
            model.epoch,
            model.max1,
            model.max2,
            model.gh1,
            model.gh2,
            &mut ghb,
        )
    };

    const IEXT: usize = 0;
    const EXT_COEFF1: f64 = 0.0;
    const EXT_COEFF2: f64 = 0.0;
    const EXT_COEFF3: f64 = 0.0;

    // Do the first calculations
    let (mut x, mut y, mut z) = shval3(
        lat,
        lon,
        alt,
        nmax,
        IEXT,
        EXT_COEFF1,
        EXT_COEFF2,
        EXT_COEFF3,
        &gha[..],
    );

    let (d, i, h, f) = dihf(x, y, z);
    let (xtemp, ytemp, ztemp) = shval3(
        lat,
        lon,
        alt,
        nmax,
        IEXT,
        EXT_COEFF1,
        EXT_COEFF2,
        EXT_COEFF3,
        &ghb[..],
    );
    let (dtemp, itemp, htemp, ftemp) = dihf(xtemp, ytemp, ztemp);

    let mut ddot = (dtemp - d).to_degrees();
    if ddot > 180.0 {
        ddot -= 360.0
    };
    if ddot <= -180.0 {
        ddot += 360.0
    };
    ddot *= 60.0;

    let idot = (itemp - i).to_degrees() * 60.0;
    let mut d = d.to_degrees();
    let i = i.to_degrees();
    let hdot = htemp - h;
    let mut xdot = xtemp - x;
    let mut ydot = ytemp - y;
    let zdot = ztemp - z;
    let fdot = ftemp - f;

    /* Deal with geographic and magnetic poles */
    if h < 100.0
    /* at magnetic poles */
    {
        d = f64::NAN;
        ddot = f64::NAN;
        /* while rest is ok */
    }

    /* at geographic poles */
    if (90.0 - lat.abs() <= 0.001) {
        x = f64::NAN;
        y = f64::NAN;
        d = f64::NAN;
        xdot = f64::NAN;
        ydot = f64::NAN;
        ddot = f64::NAN;
        /* while rest is ok */
    }

    Ok(Field {
        d,
        i,
        h,
        f,
        x,
        y,
        z,
        ddot,
        fdot,
        hdot,
        idot,
        xdot,
        ydot,
        zdot,
    })
}

fn decimal_day_of_year(date: time::Date) -> f64 {
    let day = date.day() as i32;
    let month = date.month() as i32;
    let year = date.year();

    const DAYS: [i32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let leap_year = if ((year % 4) == 0) && (((year % 100) != 0) || ((year % 400) == 0)) {
        1
    } else {
        0
    };
    let mut day_in_year = DAYS[(month - 1) as usize] + day;
    if month > 2 {
        day_in_year += leap_year;
    }
    f64::from(year) + (f64::from(day_in_year) / (365.0 + f64::from(leap_year)))
}

/// Extrapolates linearly a spherical harmonic model with a rate-of-change model.
fn extrapsh(
    date: f64,
    dte1: f64,
    nmax1: usize,
    nmax2: usize,
    gh1: &[f64],
    gh2: &[f64],
    gh: &mut [f64],
) -> usize {
    let factor = date - dte1;
    let (k, nmax) = match nmax1.cmp(&nmax2) {
        Ordering::Equal => {
            let k = nmax1 * (nmax1 + 2);
            let nmax = nmax1;
            (k, nmax)
        }
        Ordering::Greater => {
            let k = nmax2 * (nmax2 + 2);
            let l = nmax1 * (nmax1 + 2);
            gh[(k + 1)..l].clone_from_slice(&gh1[(k + 1)..l]);
            let nmax = nmax1;
            (k, nmax)
        }
        Ordering::Less => {
            let k = nmax1 * (nmax1 + 2);
            let l = nmax2 * (nmax2 + 2);
            for i in (k + 1)..l {
                gh[i] = factor * gh2[i];
            }
            let nmax = nmax2;
            (k, nmax)
        }
    };

    for i in 1..k {
        gh[i] = gh1[i] + factor * gh2[i];
    }

    nmax
}

/// Interpolates linearly, in time, between two spherical harmonic models.
fn interpsh(
    date: f64,
    dte1: f64,
    nmax1: usize,
    dte2: f64,
    nmax2: usize,
    gh1: &[f64],
    gh2: &[f64],
    gh: &mut [f64],
) -> usize {
    let factor = (date - dte1) / (dte2 - dte1);
    let (k, nmax) = match nmax1.cmp(&nmax2) {
        Ordering::Equal => {
            let k = nmax1 * (nmax1 + 2);
            let nmax = nmax1;
            (k, nmax)
        }
        Ordering::Greater => {
            let k = nmax2 * (nmax2 + 2);
            let l = nmax1 * (nmax1 + 2);
            for i in (k + 1)..l {
                gh[i] = gh1[i] + factor * -gh1[i];
            }
            let nmax = nmax1;
            (k, nmax)
        }
        Ordering::Less => {
            let k = nmax1 * (nmax1 + 2);
            let l = nmax2 * (nmax2 + 2);
            for i in (k + 1)..l {
                gh[i] = factor * gh2[i];
            }
            let nmax = nmax2;
            (k, nmax)
        }
    };

    for i in 1..k {
        gh[i] = gh1[i] + factor * (gh2[i] - gh1[i]);
    }

    nmax
}

/// Calculates field components from spherical harmonic (sh) models.
///
/// Based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin, report no. 71/1, institute
/// of geological sciences, U.K.
fn shval3(
    flat: f64,
    flon: f64,
    elev: f64,
    nmax: usize,
    iext: usize,
    ext1: f64,
    ext2: f64,
    ext3: f64,
    gh: &[f64],
) -> (f64, f64, f64) {
    let earth_radius = 6371.2;
    let dtr = 0.01745329;
    // Square of the semi-major axes of the WGS84 reference sphereoid used for transforming between
    // geodetic and geocentric coordinates.
    let a2 = 40680631.59;
    // Square of the semi-minor axes of the WGS84 reference sphereoid used for transforming between
    // geodetic and geocentric coordinates.
    let b2 = 40408299.98;

    let r = elev;
    let argument = flat * dtr;
    let slat = argument.sin();

    let aa = if (90.0 - flat) < 0.001 {
        // 300 ft. from North pole
        89.999
    } else if (90.0 + flat) < 0.001 {
        // 300 ft. from South pole
        -89.999
    } else {
        flat
    };

    let clat = (aa * dtr).cos();

    let mut sl = [0f64; 14];
    let mut cl = [0f64; 14];

    let argument = flon * dtr;
    sl[1] = argument.sin();
    cl[1] = argument.cos();

    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;

    let sd = 0.0;
    let cd = 1.0;
    let mut l = 1;
    let mut n: usize = 0;
    let mut m: usize = 1;
    let npq = (nmax * (nmax + 3)) / 2;

    // to geocentric
    let aa = a2 * clat * clat;
    let mut bb = b2 * slat * slat;
    let mut cc = aa + bb;
    let argument = cc;
    let dd = argument.sqrt();
    let argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
    let r = argument.sqrt();
    let cd = (elev + dd) / r;
    let sd = (a2 - b2) / dd * slat * clat / r;
    let aa = slat;
    let slat = slat * cd - clat * sd;
    let clat = clat * cd + aa * sd;

    let ratio = earth_radius / r;
    let mut aa = (3.0f64).sqrt();

    let mut p = [0f64; 119];
    let mut q = [0f64; 119];

    p[1] = 2.0 * slat;
    p[2] = 2.0 * clat;
    p[3] = 4.5 * slat * slat - 1.5;
    p[4] = 3.0 * aa * clat * slat;
    q[1] = -clat;
    q[2] = slat;
    q[3] = -3.0 * clat * slat;
    q[4] = aa * (slat * slat - clat * clat);

    let mut fn_ = 0.0;
    let mut fm = 0.0;

    for k in 1..=npq {
        let mut rr = 0.0;

        if n < m {
            m = 0;
            n = n + 1;
            rr = ratio.powi(n as i32 + 2);
            fn_ = n as f64;
        }

        fm = m as f64;

        if k >= 5 {
            if m == n {
                let argument = (1.0 - 0.5 / fm);
                aa = argument.sqrt();
                let j = k - n - 1;
                p[k] = (1.0 + 1.0 / fm) * aa * clat * p[j];
                q[k] = aa * (clat * q[j] + slat / fm * p[j]);
                sl[m] = sl[m - 1] * cl[1] + cl[m - 1] * sl[1];
                cl[m] = cl[m - 1] * cl[1] - sl[m - 1] * sl[1];
            } else {
                let argument = fn_ * fn_ - fm * fm;
                aa = argument.sqrt();
                let argument = ((fn_ - 1.0) * (fn_ - 1.0)) - (fm * fm);
                bb = argument.sqrt() / aa;
                cc = (2.0 * fn_ - 1.0) / aa;
                let ii = k - n;
                let j = k - 2 * n + 1;
                p[k] = (fn_ + 1.0) * (cc * slat / fn_ * p[ii] - bb / (fn_ - 1.0) * p[j]);
                q[k] = cc * (slat * q[ii] - clat / fn_ * p[ii]) - bb * q[j];
            }
        }

        aa = rr * gh[l];

        if (m == 0) {
            x += aa * q[k];
            z -= aa * p[k];
            l = l + 1;
        } else {
            bb = rr * gh[l + 1];
            cc = aa * cl[m] + bb * sl[m];
            x = x + cc * q[k];
            z = z - cc * p[k];
            if clat > 0.0 {
                y += (aa * sl[m] - bb * cl[m]) * fm * p[k] / ((fn_ + 1.0) * clat);
            } else {
                y += (aa * sl[m] - bb * cl[m]) * q[k] * slat;
            }
            l = l + 2;
        }
        m = m + 1;
    }

    if (iext != 0) {
        aa = ext2 * cl[1] + ext3 * sl[1];
        x = x - ext1 * clat + aa * slat;
        y = y + ext2 * sl[1] - ext3 * cl[1];
        z = z + ext1 * slat + aa * clat;
    }

    aa = x;
    x = x * cd + z * sd;
    z = z * cd - aa * sd;

    (x, y, z)
}

/// Computes the geomagnetic declination (d), inclination (i), horizontal field strength (h), and
/// total field strength (f) from x, y, and z.
fn dihf(x: f64, y: f64, z: f64) -> (f64, f64, f64, f64) {
    let mut sn = 0.0;
    let mut hpx = 0.0;

    sn = 0.0001;

    let mut d = 0.0;
    let mut i = 0.0;
    let mut h = (x * x + y * y).sqrt(); /* calculate horizontal intensity */
    let mut f = (x * x + y * y + z * z).sqrt(); /* calculate total intensity */
    if (f < sn) {
        // If d and i cannot be determined, set them equal to NaN.
        d = f64::NAN;
        i = f64::NAN;
    } else {
        i = z.atan2(h);
        if (h < sn) {
            d = f64::NAN;
        } else {
            hpx = h + x;
            if (hpx < sn) {
                d = std::f64::consts::PI;
            } else {
                d = 2.0 * y.atan2(hpx);
            }
        }
    }

    (d, i, h, f)
}

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum Error {
    #[error("Year must be between 1900 and 2030")]
    InvalidYear,
    #[error("Latitude is out of bounds")]
    LatitudeOutOfBounds,
    #[error("Longitude is out of bounds")]
    LongitudeOutOfBounds,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_declination() {
        // TODO: compare with
        // geomag70.exe IGRF13.COF 1999.0 D M100 41.60911536877931 41.6009718546448
        let field = declination(
            41.60911536877931,
            41.6009718546448,
            100,
            time::Date::try_from_yo(1999, 1).unwrap(),
        );
        assert_eq!(field, Err(Error::InvalidYear));
    }
}
