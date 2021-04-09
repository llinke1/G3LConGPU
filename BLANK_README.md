<!--
***
***
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** github_username, repo_name, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />

  <h3 align="center">G3LConGPU</h3>

  <p align="center">
    Pipeline that measures the Galaxy-Galaxy-Galaxy-Lensing (G3L) Correlation function ON GPUs, provided galaxy catalogs
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This code can measure the Galaxy-Galaxy-Galaxy lensing correlation function and converts it into the corresponding aperture statistics. Given a catalog of lens and source galaxies, it can calculate the correlation function either with kd-trees on a standard CPU (parallelized over all available Cores), or brute-force on a GPU. It was used for the analysis in <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..13L/abstract">Linke et al.(2020a)</a> and <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...640A..59L/abstract"> Linke et al.(2020b)</a>. To learn more about Galaxy-Galaxy-Galaxy-Lensing, you can also check out <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...622A.104S/abstract"> Simon et al.(2019)</a>. 

**To avoid retyping too much info. Do a search and replace with your text editor for the following:**
`github_username`, `repo_name`, `twitter_handle`, `email`, `project_title`, `project_description`


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites
To use this code solely with the CPU version, these simple requirements are all that is needed:
* **g++** (Tested for version 9.3.0). 
Under Ubuntu this can be installed with
  ```sh
  sudo apt install build-essential
  ```
* **bash**. Should be available under Linux distributions, for use under Windows, consult how to create a Windows Subsytem for Linux (e.g. here [https://docs.microsoft.com/de-de/windows/wsl/install-win10].
* **openMP** (Tested for version 4.5). Under Ubuntu this can be installed with
```sh
sudo apt-get install libomp-dev
```
* **GNU Scientific Library** (Tested for version 2.6). Check here for how to install it: [https://www.gnu.org/software/gsl/]
* **python** (Tested for version 3.8.5, at least 3 is needed!). On Ubuntu this can be installed with
```sh
sudo apt-get install python3.8
```

To use the GPU accelerated version, additionally the following is needed:

* **NVIDIA graphics card with CUDA capability of at least 2**. Check here to see, if your card works: [https://en.wikipedia.org/wiki/CUDA#GPUs_supported].
* **CUDA SDK Toolkit** (Tested for version 10.1, at least 7 needed!)
Can be downloaded here [https://developer.nvidia.com/accelerated-computing-toolkit]

In general, some knowledge on CUDA and how GPUs work is useful to understand the code!

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/github_username/repo_name.git
   ```
2. Install missing prerequisites
If you only want to use the CPU version, go to step 5
3. Go into source directory and open Makefile
```sh
cd src
xdg-open Makefile
```
4. Adapt the `-arch=` parameter to the architecture of your graphics card.
5. run `make`, if you only install the CPU part **or** run `make all` to install all
6. Now, check if the folder `bin` was created and contains the necessary executables.


<!-- USAGE EXAMPLES -->
## Usage

### Input
#### Galaxy Catalogs
To use the code, you need at least three galaxy catalogs, one (or two) for the lens galaxies, one for the source galaxies, and one (or two) for randoms, that mimic the selection function of the lens galaxies, but are uncorrelated. These catalogs should be ASCII-files, where each line corresponds to a different galaxy, and the columns denote:
1. x-position in arcmin
2. y-position in arcmin
3. First component of galaxy ellipticity
4. Second component of galaxy ellipticity
5. Galaxy redshift
6. Weight of galaxy ellipticity (most likely the output of the shear estimation code like lensfit)

For source galaxies the redshift is ignored, for random and lens galaxies the ellipticities and the weight is ignored. If no redshift weighting of lens pairs is used,
the lens redshift is also ignored. 

#### Source redshift distribution
To calculate the "physical" aperture statistics weighted by the critical surface mass density, the redshift distribution of source galaxies is needed. This needs to be provided as an ASCII file, where the first column lists z and the second column n(z). 

### Catalog --> Aperturestatistics
To determine the aperture statistics from the provided catalogs, a `run` script needs to be started. An example is provided in `run`. This script does the following:
1. It estimates the angular two-point correlation function of the lens galaxies using the Landy-Szalay Estimator for auto- and the Szalay-Szapudi Estimator for cross-correlations.
2. It calculates the average critical surface mass density for the given source redshift distribution for a range of lens redshifts
3. It calculates the angular diameter distance for a range of lens redshifts
4. It estimates the G3L correlation function Gtilde using the estimator in [Linke et al(2020a)](https://ui.adsabs.harvard.edu/abs/2020A%26A...640A..59L/abstract)
5. It converts the Gtilde into aperture statistics, using the formula in [Schneider & Watts (2005)](https://ui.adsabs.harvard.edu/abs/2005A%26A...432..783S/abstract)

To do all of this, `run_example` needs to be adapted to the specific purpose. In particular, folders containing the data, the names of the data tiles, whether or not redshift weighting should be used, and which cosmology should be used for the calculation of the physical aperture statistics.


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* The kdTree code parts are based on code by P. Simon (psimon1@uni-bonn.de)
* This ReadMe is based on [https://github.com/othneildrew/Best-README-Template]




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username
