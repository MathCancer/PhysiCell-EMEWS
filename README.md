# PhysiCell-EMEWS

Software to accompany our paper submitted to [BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/):

_High-throughput cancer hypothesis testing with an integrated PhysiCell-EMEWS workflow._
Jonathan Ozik, Nicholson Collier, Justin M. Wozniak, Charles Macal, Chase Cockrell, Samuel H. Friedman, Ahmadreza Ghaffarizadeh, Randy Heiland, Gary An and Paul Macklin

The primary study of interest is the _cancer-immune_ one. To compile [PhysiCell](http://physicell.mathcancer.org/) for this model:
```
$ cd studies/cancer-immune/PhysiCell-src
$ make -f Makefile-immune

# your compiler needs to support OpenMP, so you may need to specify it if the previous fails:
$ make -f Makefile-immune CC=g++-7
```
Assuming it builds without errors, you can then test that it runs, e.g.:
```
$ mkdir run
$ cd run
$ ../cancer-immune-EMEWS    # this will provide sample parameters for a run.
```

To run [EMEWS](http://emews.org/), see the `EMEWS-scripts` directory for instructions. To perform final analysis, see the `analysis-scripts` directory.
