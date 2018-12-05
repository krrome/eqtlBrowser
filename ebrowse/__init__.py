from yaml import load
import sys
import pkg_resources

__version__ = '0.1'


def get_local_paths():
    with open(pkg_resources.resource_filename(__name__, "../local_paths.yaml"), "r") as ofh:
        config = load(ofh)
    if sys.platform == 'linux':
        import socket
        if socket.gethostname().startswith("ebi"):
            return config["ebi_testing"]
        else:
            return config["production"]
    else:
        return config["osx_testing"]

PATHS = get_local_paths()