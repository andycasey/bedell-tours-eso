import os
import re
import time
import numpy as np
from glob import glob
from bs4 import BeautifulSoup

from astroquery.eso import Eso as ESO 

cwd = os.path.dirname(os.path.realpath(__file__))

ESO_USERNAME = "andycasey"

SKIP = 0 # Skip how many batches at the start? (for if you are re-running this..)
BATCH = 2000 # How many datasets should we get per ESO request?
WAIT_TIME = 60 # Seconds between asking ESO if they have prepared our request
DATA_DIR = os.path.join(cwd, "data")
MISSING_FILES_PATH = os.path.join(cwd, "missing_files.txt")

if not os.path.exists(DATA_DIR):
    os.mkdir(DATA_DIR)

with open(MISSING_FILES_PATH, "r") as fp:
    missing_files = [line.strip() for line in fp.readlines()]

# Check to see if we have this filename already.
records = ["SAF+{}.tar".format(each[:-12]) \
    for each in missing_files \
    if not os.path.exists(os.path.join(DATA_DIR, each[:-12] + ".tar"))]

N = len(records)
I = int(np.ceil(N/float(BATCH)))

remote_paths = []

print("In total we will request {} records in {} requests".format(N, I))

for i in range(I):

    if i < SKIP:
        print("Skipping {}".format(i + 1))
        continue

    print("Starting with batch number {}/{}".format(i + 1, I))

    data = [("dataset", dataset) for dataset in records[i*BATCH:(i + 1)*BATCH]]

    # Login to ESO.
    eso = ESO()
    eso.login(ESO_USERNAME)

    prepare_response = eso._session.request("POST",
        "http://dataportal.eso.org/rh/confirmation", data=data)
    assert prepare_response.ok

    # Additional payload items required for confirmation.
    data += [
        ("requestDescription", ""),
        ("deliveryMediaType", "WEB"), # OR USB_DISK --> Holy shit what the fuck!
        ("requestCommand", "SELECTIVE_HOTFLY"),
        ("submit", "Submit")
    ]

    confirmation_response = eso._session.request("POST", 
        "http://dataportal.eso.org/rh/requests/{}/submission".format(ESO_USERNAME),
        data=data)
    assert confirmation_response.ok

    # Parse the request number so that we can get a download script from ESO later
    _ = re.findall("Request #[0-9]+\w", confirmation_response.text)[0].split()[-1]
    request_number = int(_.lstrip("#"))

    print("Retrieving remote paths for request number {}/{}: {}".format(
        i + 1, I, request_number))

    # Check if ESO is ready for us.
    while True:

        url = "https://dataportal.eso.org/rh/requests/{}".format(ESO_USERNAME)
        check_state = eso._request("GET", url, cache=False)
        root = BeautifulSoup(check_state.text, "html5lib")

        link = root.find(href="/rh/requests/{}/{}".format(
            ESO_USERNAME, request_number))

        image = link.find_next("img")
        state = image.attrs["alt"]

        print("Current state {} on request {} ({}/{})".format(
            state, request_number, i + 1, I))

        if state != "COMPLETE":

            # Remove anything from the astroquery cache.
            for cached_file in glob(os.path.join(eso.cache_location, "*")):
                os.remove(cached_file)

            print("Sleeping for {} seconds..".format(WAIT_TIME))
            time.sleep(WAIT_TIME)

        else:
            break

    response = eso._request("GET", "{}/{}/script".format(url, request_number))
    
    paths = response.text.split("__EOF__")[-2].split("\n")[1:-2]
    print("Found {} remote paths for request_number {}".format(
        len(paths), request_number))
    remote_paths.extend(paths)

    # Remove anything from the astroquery cache.
    for cached_file in glob(os.path.join(eso.cache_location, "*")):
        os.remove(cached_file)
    
    # We have all the remote paths for this request. At ESO's advice, let's
    # wait another 60 seconds before starting our new request.
    if I > i + 1:
        time.sleep(60)

# Prepare the script for downloading.
template_path = os.path.join(cwd, "download.sh.template")
with open(template_path, "r") as fp:
    contents = fp.read()

script_path = os.path.join(DATA_DIR, "download.sh")
with open(script_path, "w") as fp:
    fp.write(contents.replace("$$REMOTE_PATHS$$", "\n".join(remote_paths))\
                     .replace("$$ESO_USERNAME$$", ESO_USERNAME))

print("Created script {0}".format(script_path))
print("Now run `cd {}; sh {}` and enter your ESO password when requested."\
      .format(DATA_DIR, os.path.basename(script_path)))

