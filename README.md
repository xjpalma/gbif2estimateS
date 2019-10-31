# gbif2EstimateS


[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

*gbif2EstimateS* is a python script to convert [GBIF] output to [EstimateS] program.

  - *GBIF* - the Global Biodiversity Information Facility—is an international network and research infrastructure funded by the world’s governments and aimed at providing anyone, anywhere, open access to data about all types of life on Earth.
  - *EstimateS* is a free software application for Windows and Macintosh operating systems that computes a variety of biodiversity statistics, estimators, and indices based on biotic sampling data. Some features require species relative abundance data, others only species presence/absence data

# Features
  - use geohash with variable precision. Geohash is a public domain geocode system invented in 2008 by Gustavo Niemeyer, which encodes a geographic location into a short string of letters and digits
  - use daterange filter
  - use geographic filter by shapefile, geojson or bbox
  - use species filter

You can also:
  - Import and export files to a especific location directory
  - use verbose mode to see what's happening


# Tech

*gbif2EstimateS* is written in Python version 3 and of course *gbif2EstimateS* itself is open source with a [public repository][gbif2EstimateS] on GitHub.

## Installation

gbif2EstimateS requires [Python](https://www.python.org) v3+ to run.

Just install the dependencies and run it.
Open your favorite Terminal and run these commands. Replace <dir_to_install> to you instalation directory
```sh
$ cd <dir_to_install>
$ git clone https://github.com/xjpalma/gbif2estimateS.git
$ cd gbif2EstimateS

$ pip install -r dist/requirements.txt
```
We strongly advise to use a [conda] environment. Install the dependencies with conda and run gbif2EstimateS in it to not mess up your system
```sh
$ conda install --file dist/requirements.txt
```
Run test to see if everything is installed correctly:
```
$ cd test
$ ./test.sh

 1. build grid
 2. init big data
 3. read each gbif observation (be patient...)
 ============================================================ 100.0% ...
   observations outside grid, time or selected species: 530
   observations wrong format (no date): 7
   observations repeated: 1142
   observations rejected: 1679
   observations accepted: 6890
   observations total: 8569

 4. process big data and output results

 == Metacommunities with more then 2 individuals:
   774 communities
   25 species
   5542 individuals

 Time elapsed: 0:00:25.425246
```
# Help, Bugs, Feedback
If you need help with gbif2EstimateS, you can hang out by mail: <xjpalma@gmail.com>. To report bugs, open a github issue and explain what is the problem.
Want to contribute? Great! Contact me and I can help you understand and update the script. Also, you are free to make changes and instantanously to see your updates!

# License
MIT
**Free Software, Hell Yeah!**

[GBIF]: <https://www.gbif.org>
[EstimateS]: <http://viceroy.eeb.uconn.edu/estimates>
[gbif2EstimateS]: <https://github.com/xjpalma/gbif2estimateS>
[git-repo-url]: <https://github.com/xjpalma/gbif2estimateS.git>
[conda]: <https://docs.conda.io/en/latest/miniconda.html>
