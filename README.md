Bedell tours ESO!
-----------------
Superstar Megan Bedell is touring the ESO archive! Get your tickets online!


Touring plan
------------
1. If you don't already have them, install `astroquery` and `BeautifulSoup`.

2. Update `prepare_download.py` to use your ESO username so I don't get blamed
   when the ESO archive crashes (again).

3. Run `python prepare_download.py`

4. Do as instructed by the script: `cd data; sh download.sh`, and enter your
   ESO password when asked.

5. When ESO asks, deny all responsibility.


Retrieving data on a single star
--------------------------------

Use the `query_object.py` script:

`python query_object.py <ra> <dec> --name my_star`

Use the `-h` flag to get more informatioon on the options available.

`python query_object.py -h`

You will need to provide a file with your ESO credentials, which by default is set to `eso_credentials.yml`.
The contents of the ESO credentials file should be something like:

````
eso:
  username: bedell
  password: pupper
````
