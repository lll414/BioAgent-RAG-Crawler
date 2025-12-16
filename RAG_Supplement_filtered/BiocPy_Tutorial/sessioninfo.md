Source URL: https://biocpy.github.io/tutorial/chapters/sessioninfo.html
Date Scraped: 2025-12-15

---

The code base for this repository is available at <https://github.com/BiocPy/tutorial>.

This book is automatically built to identify issues caused by changes in the dependencies. We use [quarto](https://quarto.org/docs/books) for reproducing the scripts provided in this book.

```
import sys
import subprocess

from rich import print
```

Lets make sure we have all packages we need for this section

```
sys.version_info
```

```
sys.version_info(major=3, minor=9, micro=18, releaselevel='final', serial=0)
```

```
import math

import biocframe
import biocutils
import genomicranges
import summarizedexperiment
import singlecellexperiment
import multiassayexperiment
import rds2py
import session_info


session_info.show()
```

Click to view session information

```
-----
biocframe                   0.5.9
biocutils                   0.1.5
genomicranges               0.4.12
multiassayexperiment        0.4.2
rds2py                      0.4.2
rich                        NA
session_info                1.0.0
singlecellexperiment        0.4.4
summarizedexperiment        0.4.3
-----
```

Click to view modules imported as dependencies

```
PIL                         10.2.0
anyio                       NA
arrow                       1.3.0
asttokens                   NA
attr                        23.2.0
attrs                       23.2.0
babel                       2.14.0
certifi                     2024.02.02
cffi                        1.16.0
charset_normalizer          3.3.2
comm                        0.2.2
cycler                      0.12.1
cython_runtime              NA
dateutil                    2.9.0.post0
debugpy                     1.8.1
decorator                   5.1.1
defusedxml                  0.7.1
exceptiongroup              1.2.0
executing                   2.0.1
fastjsonschema              NA
fqdn                        NA
idna                        3.6
importlib_metadata          NA
importlib_resources         NA
ipykernel                   6.29.3
iranges                     0.2.3
isoduration                 NA
jedi                        0.19.1
jinja2                      3.1.3
json5                       0.9.22
jsonpointer                 2.4
jsonschema                  4.21.1
jsonschema_specifications   NA
jupyter_events              0.9.1
jupyter_server              2.13.0
jupyterlab_server           2.25.4
kiwisolver                  1.4.5
markupsafe                  2.1.5
matplotlib                  3.8.3
matplotlib_inline           0.1.6
mpl_toolkits                NA
nbformat                    5.10.3
numpy                       1.26.4
overrides                   NA
packaging                   24.0
pandas                      2.2.1
parso                       0.8.3
pexpect                     4.9.0
platformdirs                4.2.0
prometheus_client           NA
prompt_toolkit              3.0.43
psutil                      5.9.8
ptyprocess                  0.7.0
pure_eval                   0.2.2
pydev_ipython               NA
pydevconsole                NA
pydevd                      2.9.5
pydevd_file_utils           NA
pydevd_plugins              NA
pydevd_tracing              NA
pygments                    2.17.2
pyparsing                   3.1.2
pythonjsonlogger            NA
pytz                        2024.1
referencing                 NA
requests                    2.31.0
rfc3339_validator           0.1.4
rfc3986_validator           0.1.1
rpds                        NA
scipy                       1.12.0
send2trash                  NA
six                         1.16.0
sniffio                     1.3.1
stack_data                  0.6.3
tornado                     6.4
traitlets                   5.14.2
typing_extensions           NA
uri_template                NA
urllib3                     2.2.1
wcwidth                     0.2.13
webcolors                   1.13
websocket                   1.7.0
yaml                        6.0.1
zipp                        NA
zmq                         25.1.2
zoneinfo                    NA
```

 

```
-----
IPython             8.18.1
jupyter_client      8.6.1
jupyter_core        5.7.2
jupyterlab          4.1.5
notebook            7.1.2
-----
Python 3.9.18 (main, Aug 28 2023, 08:38:32) [GCC 11.4.0]
Linux-6.5.0-1016-azure-x86_64-with-glibc2.35
-----
Session information updated at 2024-03-15 22:54
```