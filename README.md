# repack
Re-pack and Compress Line-transition Data for Ratiative-tranfer Calculations

This code identifies the strong lines that dominate the spectrum from
the large-majority of weaker lines.  The code returns a binary
line-by-line (LBL) file with the strong lines info (wavenumber, Elow,
gf, and isotope ID), and an ascii file with the combined contribution
of the weaker lines compressed into a continuum extinction coefficient
(in cm-1 amagat-1) as function of wavenumber and temperature.

Currently available databases:
* ExoMol (http://www.exomol.com/)
* HITRAN (https://www.cfa.harvard.edu/hitran/)
* Kurucz's TiO (http://kurucz.harvard.edu/molecules/tio)

### Table of Contents
* [Team Members](#team-members)
* [Install and Compile](#install-and-compile)
* [Getting Started](#getting-started)
* [ Notes](#notes)
* [Be Kind](#be-kind)
* [License](#license)

### Team Members
* [Patricio Cubillos](https://github.com/pcubillos/) (IWF) <patricio.cubillos@oeaw.ac.at>

### Install and Compile
``repack`` is compatible with both Python2 and Python3, and runs (at least) in both Linux and OSX.  
To obtain the ``repack`` code, clone this repository to your local machine with the following terminal commands:  
```shell
# Clone the repository to your working directory:  
git clone https://github.com/pcubillos/repack/

# Compile the C-extensions of the program:
cd repack
make  
```


### Getting Started

The following script processes the Exomol HCN line-transition data.  Returning a LBL and continuum files with compressed info.

```shell
# Go to the repositoy demo folder:
cd demo/
# Download and uncompress Exomol HCN data:
wget -i wget_exomol_hcn.txt
bzip2 -d *.bz2

# Call the repack command-line executable for the HCN demo config file:
../repack.py repack_HCN.cfg
```

### Notes
- The input configuration file specifies the wavenumber and temperature sampling arrays, the input files, and the line-intensity threshold.  
- The binary LBL output file contains the wavenumber (in cm-1), energy of the lower state (in cm-1), the gf value, and the isotope index, encoded as three doubles and one integer, storing the info for each line after the other.  
- The ascii continuum output file contains the extinction-coefficient data (in cm-1 amagat-1) as function of wavenumber and temperature.  


### Be Kind

Please, be kind and acknowledge the effort of the authors by citing the article asociated to this project:  

  [Cubillos (2017): An Algorithm to Compress Line-transition Data for Radiative-transfer Calculations](http://adsabs.harvard.edu/abs/2017ApJ...850...32C), ApJ 850, 32.  


### License

Copyright (c) 2017 Patricio Cubillos and contributors.
``repack`` is open-source software under the MIT license (see LICENSE).

