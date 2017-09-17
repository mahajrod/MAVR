# MAVR v0.1

## How to install on Manjaro Linux

```
sudo pacman -S python-pip # installs pip, a python package manager

sudo pip install virtualenv # installs virtualenv, a tool for creating isolated directories for installing project-specific python dependencies

virtualenv ~/.virtualenvs/MAVR --python=python2.7 # creates a virtual environment for MAVR, where all its package dependencies will be installed. We tell it to use python 2.7

. ~/.virtualenvs/MAVR/bin/activate # activates the MAVR virtual environment. now every "pip install" command will install packages into this virtual environment (folder), and every python execution will look up packages from this environment

pip install -e git+https://github.com/mahajrod/MAVR#egg=MAVR # installs MAVR as a pip package
```

## Troubleshooting
The following packages are required by some of the pip dependencies
```
sudo pacman -S gcc # installs gcc which is required for one of the dependencies of MAVR
sudo pacman -S tk # required by matplotlib
```