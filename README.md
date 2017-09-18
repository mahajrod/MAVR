# MAVR v0.1

## How to install
1. Install `pip` (python package manager)
```bash
# Manjaro Linux
sudo pacman -S python-pip 

# Ubuntu
sudo apt install python-pip
```
2. Set up `virtualenv` (a tool for creating isolated directories for installing project-specific python dependencies)
```bash
sudo pip install virtualenv  

# create a virtual environment for MAVR, where all its package dependencies will be installed, we tell it to use python 2.7
virtualenv ~/.virtualenvs/MAVR --python=python2.7 

# activates the MAVR virtual environment, now every "pip install" command will install packages into this virtual 
# environment (folder), and every python execution will look up packages from this environment
. ~/.virtualenvs/MAVR/bin/activate 
```
3. Use `pip` to install MAVR directly from GitHub
```bash
pip install -e git+https://github.com/mahajrod/MAVR#egg=MAVR 
```

### Troubleshooting
The following packages are required by some of the pip dependencies: `gcc` and `tk`
```bash
# Manjaro Linux
sudo pacman -S gcc # installs gcc which is required for one of the dependencies of MAVR
sudo pacman -S tk # required by matplotlib

# Ubuntu
sudo apt install gcc
sudo apt install tk
```
