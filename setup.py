from distutils.core import setup


#primary_deps mongoengine, zlib, pyyaml, flask, setuptools
#development_deps yaml, h5py, numpy, pandas, tqdm, cyvcf2, sqlalchemy, mongoengine, json, zlib, flask, shutil, getpass

requirements = [
    "future",
    "mongoengine==0.16.1",
    "zlib",
    "pyyaml",
    "flask==1.0.2",
    "setuptools",
]


setup(
    name='eqtlBrowser',
    version='0.1',
    packages=['ebrowse'],
    install_requires=requirements,
    url='',
    license='MIT',
    author='Roman Kreuzhuber',
    author_email='',
    description=''
)
