from yaml import load
import sys
import pkg_resources
import os

__version__ = '0.1'

def get_exec_condition():
    if sys.platform == 'linux':
        import socket
        hostname = socket.gethostname()
        if hostname.startswith("ebi") or hostname.startswith("gpu-"):
            return "ebi_testing"
        else:
            return "production"
    else:
        return "osx_testing"

def get_local_paths():
    with open(pkg_resources.resource_filename(__name__, "../local_paths.yaml"), "r") as ofh:
        config = load(ofh)
    return config[get_exec_condition()]

PATHS = get_local_paths()
host_interfix = ""
if 'INTERFIX' in os.environ and os.environ['INTERFIX'] != "":
    host_interfix = "/"+os.environ['INTERFIX'].strip("/")

if 'MONGO_HOST' in os.environ and os.environ['MONGO_HOST'] != "":
    PATHS["mongo_host"] = os.environ['MONGO_HOST']

max_req_per_sec = 20
max_req_per_min = 200

if 'REQ_PER_MIN' in os.environ and os.environ['REQ_PER_MIN'] != "":
    max_req_per_min = int(os.environ['REQ_PER_MIN'])

if 'REQ_PER_SEC' in os.environ and os.environ['REQ_PER_SEC'] != "":
    max_req_per_sec = int(os.environ['REQ_PER_SEC'])
