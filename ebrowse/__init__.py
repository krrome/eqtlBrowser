from yaml import load
import sys

def get_local_paths():
    with open("./local_paths.yaml", "r") as ofh:
        config = load(ofh)
    if sys.platform == 'linux':
        return config["ebi_testing"]
    else:
        return config["osx_testing"]

PATHS = get_local_paths()