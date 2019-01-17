
import re
import os
import requests
import tarfile
import tempfile
import yaml
from requests.auth import HTTPBasicAuth

from astropy.io import fits
from astropy.table import Table
from collections import OrderedDict
from time import sleep

from bs4 import BeautifulSoup

from time import mktime
from datetime import datetime

def convert_obs_time(obs_time, fmt="%Y-%m-%dT%H:%M:%S.%f"):
    return mktime(datetime.strptime(obs_time, fmt).timetuple())


class Harps(object):

    row_limit = 100000

    def __init__(self, credentials_path):
        self._credentials_path = credentials_path

        with open(credentials_path, "r") as fp:
            credentials = yaml.load(fp)["eso"]
            self._eso_credentials = (credentials["username"], credentials["password"])
            
        self.login(*self._eso_credentials)

        return None



    def login(self, username, password):

        self.session = requests.Session()
        prepare = self.session.get("https://www.eso.org/sso")

        root = BeautifulSoup(prepare.content, 'html5lib')
        login_input = root.find(name='input', attrs={'name': 'execution'})

        if login_input is None:
            raise ValueError("ESO page did not have the correct attribute")

        login = self.session.post(prepare.url, 
                                  data=dict(
                                     username=username,
                                     password=password,
                                     execution=login_input.get("value"),
                                     _eventId="submit",
                                     geolocation=""
                                 ))
        login.raise_for_status()

        root = BeautifulSoup(login.content, "html5lib")
        authenticated = (root.find("h4").text == "Login successful")

        assert authenticated



    def query_position(self, ra, dec, **kwargs):
        return self._query(coord1=("", f"{ra}"), coord2=("", f"{dec}"))

    def query_target(self, target_name, **kwargs):
        return self._query(target=("", f"{target_name}"))
        

    def _query(self, **params):
        payload = OrderedDict([
            ("wdbo", ("", "html/display")),
            ("max_rows_returned", ("", "{:.0f}".format(self.row_limit))),
            ("target", ("", "")),
            ("resolver", ("", "simbad")),
            ("wdb_input_file", ("", "", "application/octet-stream")),
            ("coord_sys", ("", "eq")),
            ("coord1", ("", "")),
            ("coord2", ("", "")),
            ("box", ("", "02 09 00")),
            ("tab_ra", ("", "on")),
            ("tab_dec", ("", "on")),
            ("tab_filter", ("", "on")),
            ("filter", ("", "Any")),
            ("tab_wavelength", ("", "on")),
            ("wavelength", ("", "Any")),
            ("tab_dataproduct_type", ("", "on")),
            ("dataproduct_type", ("", "Any")),
            ("tel_id", ("", "Any")),
            ("tab_ins_id", ("", "on")),
            ("ins_id", ("", "HARPS")),
            ("obstech", ("", "Any")),
            ("tab_date_obs", ("", "on")),
            ("date_obs", ("", "")),
            ("mjd_obs", ("", "")),
            ("tab_exptime", ("", "on")),
            ("exptime", ("", "")),
            ("multi_ob", ("", "%")),
            ("tab_collection_name", ("", "on")),
            ("tab_prog_id", ("", "on")),
            ("prog_id", ("", "")),
            ("username", ("", "")),
            ("p3orig", ("", "%")),
            ("tab_origfile", ("", "on")),
            ("origfile", ("", "")),
            ("tab_dp_id", ("", "on")),
            ("dp_id", ("", "")),
            ("rel_date", ("", "")),
            ("tab_referenc", ("", "on")),
            ("referenc", ("", "")),
            ("batch_id", ("", "")),
            ("publication_date", ("", "")),
            ("wdb_input_file_raw", ("", "", "application/octet-stream")),
            ("order_main", ("", "dummy"))
        ])
        payload.update(params)

        response = self.session.post(
            "http://archive.eso.org/wdb/wdb/adp/phase3_main/query", 
            files=payload)

        content = str(response.content)

        rows = "\n".join([r for r in content.split("\n") if "PHASE3+" in r])

        names = ("Mark", "More", "ARCFILE", "HDR", "Object", "RA", "DEC", "Filter",
                "ABMAGLIM", "Wavelength", "SNR", "Resolution", "Product category",
                "Instrument", "Date Obs", "Exptime", "Collection", "Product version",
                "Release Description", "Run/Program ID", "ORIGFILE", "REFERENCE Catalog",
                "Interface")
        delete_names = ("Mark", "More", "HDR", "Filter", "ABMAGLIM",
            "Product category", "Collection", "Product version", "Release Description",
            "Interface", "REFERENCE Catalog")

        _ = "<TR"
        rows = _ + _.join(rows.replace("[doc&nbsp;id:", "[doc:").split(_)[1:])
        if rows == _:
            return Table(names=[n for n in names if n not in delete_names])

        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(bytes(f"<table>{rows}</table>", encoding="utf-8"))
        f.close()


        table = Table.read(f.name, format="ascii.html", names=names)

        # Delete unnecessary columns.
        for column_name in delete_names:
            del table[column_name]

        # Parse the PHASE3 identifiers.
        table["dataset_identifier"] = re.findall(
            "PHASE3\+[0-9]+\+ADP\.[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}", 
            content)

        # Add obs_seconds?
        table["date_obs_seconds"] = [convert_obs_time(ea) for ea in table["Date Obs"]]

        return table




    def _prepare_dataset_request(self, dataset_identifier):

        if isinstance(dataset_identifier, str):
            dataset_identifier = [dataset_identifier]

        data = [("dataset", dataset) for dataset in dataset_identifier]

        prepare = self.session.post(
            "http://dataportal.eso.org/rh/confirmation", data=data)

        # Additional payload items required for confirmation.
        data += [
            ("requestDescription", ""),
            ("deliveryMediaType", "WEB"), 
            ("requestCommand", "SELECTIVE_HOTFLY"),
            ("submit", "Submit")
        ]

        username = self._eso_credentials[0]
        confirm = self.session.post(self._data_portal_api_end_point(
                f"rh/requests/{username}/submission"),
                data=data)
        confirm.raise_for_status()

        _ = re.findall("Request #[0-9]+\w", confirm.text)[0].split()[-1]
        request_number = int(_.lstrip("#"))

        return request_number



    def _get_dataset_state(self, request_number):

        username = self._eso_credentials[0]
        status = self.session.get(self._data_portal_api_end_point(
            f"rh/requests/{username}/recentRequests"))

        idx = status.text.index(f"/rh/requests/{username}/{request_number}")
        state = status.text[idx:].split("<img")[1].split("alt=")[1].split('"')[1]

        return state


    def _get_dataset_download_script(self, request_number):
        username = self._eso_credentials[0]
        url = self._data_portal_api_end_point(
                f"rh/requests/{username}/{request_number}/script")
        r = self.session.get(url)
        return r.text
        

    def _get_dataset_remote_paths(self, request_number):
        content = self._get_dataset_download_script(request_number)
        paths = content.split("__EOF__")[-2].split("\n")[1:-2]
        return tuple([path.strip('"') for path in paths])
    

    def _data_portal_api_end_point(self, end_point):
        return f"https://dataportal.eso.org/{end_point}"


    def get_dataset_identifiers(self, dataset_identifiers, check_interval=1):

        if len(dataset_identifiers) < 1:
            return (None, [])
        
        request_number = self._prepare_dataset_request(dataset_identifiers)

        while self._get_dataset_state(request_number) != "COMPLETE":
            sleep(check_interval)

        return (request_number, self._get_dataset_remote_paths(request_number))




    def get_spectrum(self, observation):

        # Instead of doing this we should request it from harps.guru API end
        # point, and harps.guru should find it if it doesn't already have it.


        # This should be magical.
        request = self._prepare_dataset_request(observation["dataset_identifier"])

        while self._get_dataset_state(request) != "COMPLETE":
            sleep(1)

        remote_paths = self._get_dataset_remote_paths(request)

        for remote_path in remote_paths:
            if remote_path.endswith(".tar"):
                local_path = self._get_dataset(remote_path)
                
                with tarfile.open(local_path) as tar:
                    for member in tar.getmembers():
                        if "s1d_A" in member.name:
                            tar.extract(member)
                            image = fits.open(member.name)
        return image



    def _get_dataset(self, remote_path):

        r = self.session.get(remote_path, auth=HTTPBasicAuth(*self._eso_credentials))

        local_path = os.path.basename(remote_path)

        with open(local_path, "wb") as fp:
            fp.write(r.content)

        return local_path


    def get_remote_path(self, remote_path, local_path):

        r = self.session.get(remote_path, auth=HTTPBasicAuth(*self._eso_credentials))
        if not r.ok:
            r.raise_for_status()

        with open(local_path, "wb") as fp:
            fp.write(r.content)

        return True


# API end points.
# api.harps.guru/identifier/{date_obs}
# api.harps.guru/metadata/{identifier}
# api.harps.guru/spectrum/{identifier}

