# Description
**MAVR** is a package of various scripts related to the multiple areas of the genomics: genome assembly, 
phylogeny, comparative genomics and many others.
It is under a continuous development and includes scripts of various level: from basic manipulation with sequences to complex pipelines.
Also it is an incubator for new independent packages. 
For example, MACE, RouToolPa, ChromoDoter in the beginning were just single scripts within the MAVR. 

# Installation
MAVR is available via conda.

For installation just run command below. It should pull all the dependencies necessary for scripts itself, but you might need to install tools if they are required for scripts.

```shell
mamba install -c mahajrod routoolpa mavr
```

or clone this and RouToolPa's repositories, add them to the PYTHONPATH environment variable and install all the dependencies   
