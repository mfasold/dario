
# DARIO - a web server for the analysis of short RNAs from high throughput sequencing data

High throughput sequencing offers the opportunity to investigate the entire ncRNAome in an essentially unbiased way. Here we present a web service called 'DARIO' that helps researchers to study short read sequencing libraries in a platform independent way and predict new ncRNAs. Small noncoding RNAs, including microRNAs, snoRNAs and tRNAs, represent a diverse collection of molecules with several important biological functions. Exploring small RNA biology or characterizing differential expression profiles by sequencing and comparing small RNA transcriptomes is also an exciting possibility to get more and deeper information about the world of non-coding RNAs (ncRNAs).

## Features and Usage

Please find a description of DARIOs features as well as usage instructons
at the dedicated help page http://dario.bioinf.uni-leipzig.de/help



## Live Site / Demo

DARIO can be accessed on http://dario.bioinf.uni-leipzig.de/

![App Screenshot](src/static/tour/DARIO-Introduction_files/DARIO-Introduction.007-001.jpg)


## Tech Stack

**Worker / Data Analysis:** Perl, Python, R, Bash, WEKA

**Web-Server:** Python, Flask, Mako

## Basic Architecture

DARIO is designed to run on two computers: 

* *WEBSERVER* serves the web pages, receives the jobs/data and sends them to the WORKER.
* *WORKER* processes the jobs, i.e. runs sophisticated bioinformatic data analysis, and then copies the results back to the WEBSERVER.

The interaction of WEBSERVER and WORKER:

1. WEBSERVER/Apache provides start page dario.bioinf.uni-leipzig
2. User sets analysis parameters, server received input (incl. files) via POST
3. WEBSERVER: Uploaded files are stored in a new directory $WEBSERVER_JOBS_PATH/$JOB_HASH and 
   store $JOB_HASH in a file $NEW_PROCESS_FILENAME.
4. WEBSERVER: a running daemon process (daemon_webserver.py) copies
   this directory via ssh to WORKER and adds $JOB_HASH to the
   job queue stored in a file $JOB_QUEUE_FILENAME there (via ssh "echo") 
5. WORKER: a running daemon process (daemon_worker.py) processes jobs $JOB_QUEUE_FILENAME threaded
6. WORKER: On completion, the result is moved via ssh to $WEBSERVER_RESULTS_PATH
7. WEBSERVER: the user watches /wait/$JOB_HASH until directory $WEBSERVER_RESULTS_PATH/$JOB_HASH 
   exists, then is forwarded to the results 
8. WEBSERVER: results are deleted after a while using a cronjob


All essential configuration items / paths are in config.py.

All web endpoints on the server (e.g. dario.*.de/help) are routed / processed by dario_app.py.


## Installation, Configuration and Deployment

### Webserver 

The following packages are required:
* python3
* Python3 packages:

      pip3 install Flask mako
* Linux tools cron, ssh, ... 

Required steps:
* Clone this repository
* Edit config.py and enter the correct paths for WEBSERVER variables 
  as WEBSERVER_SOURCE_PATH. Make sure all paths defined there exist
  and are writiable by the user (and/or the web user).
* Setup SSH access without passphrase from WEBSERVER to WORKER, and back
* Setup webserver with WSGI and point the WSGI app to dario_app ([Example for mod_wsgi](https://modwsgi.readthedocs.io/en/master/user-guides/configuration-guidelines.html))
* Start the daemon, the easiest way is to run `python3 daemon_webserver.py` in a screen session

### Analysis server (WORKER)

The following packages are required:
* perl
* Java and WEKA4
* python3
* Python3 packages:

      pip3 install mako
* R (Recommendation: setup local R Library directory: `mkdir ~/R_Libs` and add it to ~/.Renviron : `R_LIBS="~/R_Libs"`)
* R packages:

      install.packages(c("ggplot2")) # >= 2.9.8
      install.packages("maptools")
      install.packages("calibrate")
      source("http://www.bioconductor.org/biocLite.R")
      biocLite("genomeIntervals")    
      install.packages(c("RColorBrewer", "plyr", "reshape", "sp"))

You can test whether the tools are properly setup by running a  
standalone data analysis (see below).
Required steps:
* Clone this repository
* Edit config.py and enter the correct paths for WORKER variables 
  as WORKER_SOURCE_PATH. Make sure all paths defined there exist
  and are writiable by the user (and/or the web user).
* Start the daemon, the easiest way is to run `python3 daemon_worker.py` in a screen session


### Standalone analysis
If you just want to run the analysis, without the webserver, only the steps from the
WORKER are needed, without starting the daemon. However, this analysis also 
requires the input be in a specific folder and in a specific format. Essentially, you
need a folder with two files, "job_params.txt" and "mapping_loci.upload". An example 
is provided in the repository.

via `python3 analysis.py -i JOB_PATH`.

## License

[![AGPL License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0)

## Project Status 

The project is on low-maintenance mode: the web service is running, but no 
additions / improvements are made. Please feel free to use and adapt this 
code to analyse your own species. Unfortunately, no support can be provided.

## Authors

- Mario Fasold (focus on python & R & bash code)
- David Langenberger (focus on perl code)


## Acknowledgements

- Felix Gärtner & Alex Donath for the Mito Webserver, where our webdesign & architecture is based on
- Steve Hoffman & Peter F. Stadler
