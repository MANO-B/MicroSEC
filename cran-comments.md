## Test environments
* local: aarch64-apple-darwin20 (macOS Sonoma 14.5), R version 4.3.2
* r-hub: linux, macos, macos-arm64, windows on Github Actions  
  https://github.com/MANO-B/MicroSEC/actions/runs/10106932404  
  
## R CMD check results
Duration: 4m 1.9s  
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

It also takes time to process test data because of the large genomic data.  
─  checking examples ... [100s/100s] OK (1m 40.1s)  
   Examples with CPU (user + system) or elapsed time > 5s  
                    user system elapsed  
   fun_homology   78.915  0.644  79.687  
   fun_read_check 15.382  0.479  15.867  

* New submission
Package was archived on CRAN  
CRAN repository db overrides:  
  X-CRAN-Comment: Archived on 2022-03-31 as check problems were not  
    corrected in time.  
  Updated: version 1.1.3 to version 2.1.3.  

## revdepcheck results
There are currently no downstream dependencies for this package.

