
""" 
Prepare HARPS data for Wobble. 
"""

import logging
import os
import argparse
import astropy.units as u
import numpy as np
import tarfile
from astropy.coordinates import (SkyCoord, EarthLocation)
from astropy.io import fits
from astropy.time import Time
from astroquery.gaia import Gaia

from harps.client import Harps

parser = argparse.ArgumentParser(description="""
Search for, retrieve, and prepare HARPS data for use with Wobble.

Given a provided source position, this script will perform a cone search against
Gaia DR2 (Brown et al. 2018) and use Gaia observables to select likely spectra
of this source from the HARPS archive. Quality constraints as to what should be
considered a good spectrum can be provided by the user. This script will then
retrieve all selected HARPS spectra, and process them into a HDF5 file that can
be used with Wobble (Bedell et al. 2019).
""")

def ValidFloatActionFactory(lower, upper):
    """ A factory for actions that will validate float inputs. """
    class ValidateAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not (upper >= values >= lower):
                raise ValueError(f"{self.dest} must be between [{lower}, {upper}]")
            setattr(namespace, self.dest, values)
    return ValidateAction

class ValidateRadius(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            values = values.to(u.arcsecond)

        except:
            logger.warning("Assuming radius is in arcseconds!")
            values = values * u.arcsecond

        setattr(namespace, self.dest, values)


parser.add_argument("ra", metavar="ra",
                    type=float, action=ValidFloatActionFactory(0, 360),
                    help="right ascension of the source [degrees]")
parser.add_argument("dec", metavar="dec",
                    type=float, action=ValidFloatActionFactory(-90, 90),
                    help="declination of the source [degrees]")
parser.add_argument("--name", metavar="source_name", type=str,
                    help="the name of the source (used for local naming only)")
parser.add_argument("--radius", metavar="radius", type=float,
                    action=ValidFloatActionFactory(0, np.inf),
                    default=5 * u.arcsecond,
                    help="cone search radius [arcseconds] to use for external "\
                         "services (e.g., Gaia). Default is 5 arcseconds.")
parser.add_argument("--gaia-adql-constraints", type=str, default="",
                    metavar="gaia_adql_constraints",
                    help="ADQL constraints to provide in ESA/Gaia query (e.g., ")
parser.add_argument("--limit", metavar="limit", type=int,
                    default=10000,
                    help="limit the number of exposures to retrieve")
parser.add_argument("--eso-credentials-path", metavar="eso_credentials_path",
                    default="eso_credentials.yml",
                    help="local path containing ESO credentials")
parser.add_argument("--working-directory", metavar="working_directory",
                    help="local working directory for files")
parser.add_argument("--verbose", "-v", action="store_true",
                    help="verbose logging")
args = parser.parse_args()

# Prepare logging.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG if args.verbose else logging.INFO)
ch.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(ch)

# Step 1: Gaia search.
input_coord = SkyCoord(ra=args.ra, dec=args.dec, unit=(u.degree, u.degree),
                       frame="icrs") # todo: allow user-specified frame
gaia_adql_constraints = f" AND {args.gaia_adql_constraints}" if args.gaia_adql_constraints else ""

logger.info(f"""Querying Gaia at {input_coord}
with radius {args.radius}
and constraints: {gaia_adql_constraints}""")

#j = Gaia.cone_search(input_coord, args.radius)

j = Gaia.launch_job(f"""SELECT DISTANCE(POINT('ICRS', ra, dec),
                                        POINT('ICRS', {args.ra}, {args.dec})) AS dist, *
                       FROM gaiadr2.gaia_source
                       WHERE CONTAINS(POINT('ICRS', ra, dec),
                                      CIRCLE('ICRS', {args.ra}, {args.dec}, {args.radius.to(u.deg).value})) = 1
                       {gaia_adql_constraints}
                       ORDER BY dist ASC""")

gaia_results = j.get_results()
G = len(gaia_results)

if G < 1:
    raise NotImplementedError("no sources found in Gaia query")

if G > 1:
    raise NotImplementedError("multiple sources found in Gaia query")

#python wobble_prepare.py 344.36658333333327 20.768833333333333
gaia_result = gaia_results[0]
if args.name is None:
    args.name = gaia_result["designation"].decode()

logger.info(f"Source identified as {gaia_result['designation'].decode()}")

gaia_ra, gaia_dec = (gaia_result["ra"], gaia_result["dec"])

# Step 3: Search HARPS.
logger.info("Logging in to ESO archive..")
harps = Harps(args.eso_credentials_path)

# todo: check that 'box' should actually be in degrees
harps_results = harps.query_position(gaia_ra, gaia_dec,
                                     box=f"{args.radius.to(u.deg)}")

logger.info(f"HARPS archive returned {len(harps_results)} rows:")
logger.info(harps_results)

if len(harps_results) < 1:
    raise NotImplementedError

# Step 4: Download files.
logger.info(f"Limiting to first {args.limit} rows")
rid, paths = harps.get_dataset_identifiers(harps_results["dataset_identifier"][:args.limit])
N = len(paths)

if args.working_directory is None:
    args.working_directory = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        args.name)

os.makedirs(args.working_directory, exist_ok=True)
logger.info(f"Downloading {N} files to {args.working_directory}")

headers = dict()
header_keys = (
    "ORIGFILE",
    "ASSON1",
    "HIERARCH ESO DRS CAL TH FILE", 
    "HIERARCH ESO DRS BERV", 
    "HIERARCH ESO DRS BJD",
    "HIERARCH ESO DRS BERVMX"
)

for i, remote_path in enumerate(paths, start=1):
    local_path = os.path.join(args.working_directory, os.path.basename(remote_path))
    logger.info(f"  ({i}/{N}): {remote_path} -> {local_path}")
    if os.path.exists(local_path):
        logger.info("  -> Skipping")

    else:
        harps.get_remote_path(remote_path, local_path)

    # If it's a FITS file, open it and get the wavelength calibration file.
    if local_path.lower().endswith(".fits"):
        with fits.open(local_path) as image:
            headers[local_path] = \
                dict([(k, image[0].header.get(k, None)) for k in header_keys])

# Download the unique calibration files.
calibration_remote_paths = []
for local_path, h in headers.items():
    for k, v in h.items():
        if k == "HIERARCH ESO DRS CAL TH FILE":
            calibration_remote_paths.append(f"SAF+{v[:-12]}.tar")
calibration_remote_paths = list(set(calibration_remote_paths))

cal_rid, cal_paths = harps.get_dataset_identifiers(calibration_remote_paths)

logger.info(f"Downloading {len(cal_paths)} calibration files")

for i, remote_path in enumerate(cal_paths, start=1):
    local_path = os.path.join(args.working_directory, os.path.basename(remote_path))
    logger.info(f"  ({i}/{len(cal_paths)}): {remote_path} -> {local_path}")
    if os.path.exists(local_path):
        logger.info("  -> Skipping")
        continue

    harps.get_remote_path(remote_path, local_path)

# Step 5: Extract all downloaded TAR files.
rv_headers = dict()

all_paths = list(paths) + list(cal_paths)
for i, path in enumerate(all_paths):
    local_path = os.path.join(args.working_directory, os.path.basename(path))
    if not local_path.endswith(".tar"): continue

    # Un tar it.
    with tarfile.open(local_path) as tar:

        for member in tar.getmembers():
            tar.extract(member)
            member_local_path = os.path.join(args.working_directory,
                                             os.path.basename(member.name))

            os.system(f"mv {member.name} {member_local_path}")



default_values = {"HIERARCH ESO DRS CAL TH FILE": ""}
logger.info("Bookkeeping...")

trim_eso_header = lambda _: _[13:] if _.startswith("HIERARCH ESO ") else _

for header_key in header_keys:

    default_value = default_values.get(header_key, np.nan)

    values = []
    for dataset_identifier in harps_results["dataset_identifier"]:

        basename = "+".join(dataset_identifier.split("+")[2:])
        local_path = os.path.join(args.working_directory, f"{basename}.fits")

        value = headers.get(local_path, dict(header_key=default_value))\
                       .get(header_key, default_value)

        values.append(value)


    harps_results[trim_eso_header(header_key)] = values

# Values to extract from the calibration file.
cal_header_keys = ("HIERARCH ESO DRS CCF RVC",
                   "HIERARCH ESO DRS DVRMS")

for h in cal_header_keys:
    harps_results[trim_eso_header(h)] = np.nan * np.ones(len(harps_results))

for i, (dataset_identifier, associated_file) \
in enumerate(zip(harps_results["dataset_identifier"], harps_results["ASSON1"])):
    
    # Get the calibration file.
    bis_g2_path = os.path.join(args.working_directory,
                               associated_file.replace("DRS_HARPS_3.5.tar", "bis_G2_A.fits"))

    if os.path.exists(bis_g2_path):
        with fits.open(bis_g2_path) as image:
            for h in cal_header_keys:
                harps_results[trim_eso_header(h)][i] = image[0].header.get(h, np.nan)
                print(h, image[0].header.get(h, np.nan))

    else:
        logger.warning(f"Could not find calibration file ({bis_g2_path}) for {dataset_identifier}")


# Step 6: Add an astropy-calculated barycentric velocity correction to each
#         header.

logger.info("Calculating barycentric velocity corrections")

observatory = EarthLocation.of_site("lasilla")
gaia_coord = SkyCoord(ra=gaia_result["ra"] * u.degree,
                      dec=gaia_result["dec"] * u.degree,
                      pm_ra_cosdec=gaia_result["pmra"] * u.mas/u.yr,
                      pm_dec=gaia_result["pmdec"] * u.mas/u.yr,
                      obstime=Time(2015.5, format="decimalyear"),
                      frame="icrs")

def apply_space_motion(coord, time):
    return SkyCoord(ra=coord.ra + coord.pm_ra_cosdec/np.cos(coord.dec.radian) * (time - coord.obstime),
                    dec=coord.dec + coord.pm_dec * (time - coord.obstime),
                    obstime=time, frame="icrs")


berv_key = "ASTROPY BERV"
harps_results[berv_key] = np.nan * np.ones(len(harps_results))

diffs = []
for i, dataset_identifier in enumerate(harps_results["dataset_identifier"]):

    basename = "+".join(dataset_identifier.split("+")[2:])
    local_path = os.path.join(args.working_directory, f"{basename}.fits")

    if os.path.exists(local_path):

        image = fits.open(local_path)

        time = Time(image[0].header["HIERARCH ESO DRS BJD"], format="jd")
        coord = apply_space_motion(gaia_coord, time)
        berv = coord.radial_velocity_correction(location=observatory)\
                    .to(u.km/u.s).value

        harps_berv = image[0].header["HIERARCH ESO DRS BERV"]
        harps_results[berv_key][i] = berv
        image[0].header[berv_key] = berv

        logger.info(f"{i}: {os.path.basename(local_path)} {time} ({coord.ra}, {coord.dec}) {harps_berv} {berv} {berv - harps_berv}")
        image.writeto(local_path, overwrite=True)
        diffs.append(berv - harps_berv)

    else:
        logger.warning(f"Skipping missing file {local_path}")

logger.info(f"Of {len(diffs)} measurements, the mean and standard deviation of "\
            f"(HARPS BERV - ASTROPY BERV) is {np.mean(diffs):.3f} km/s, {np.std(diffs):.3f} km/s")

# Write HARPS summary results to disk.
summary_path = os.path.join(args.working_directory, "summary.csv")
harps_results.write(summary_path, overwrite=True)

logger.info("Fin")