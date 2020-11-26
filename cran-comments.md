## Resubmission
Thank you for the comment:

You have examples for unexported functions.
  fun_hairpin_check() in:
     fun_hairpin_check.Rd
  fun_hairpin_trimming() in:
     fun_hairpin_trimming.Rd
  fun_repeat_check() in:
     fun_repeat_check.Rd
  fun_setting() in:
     fun_setting.Rd
  fun_support() in:
     fun_support.Rd
  fun_zero() in:
     fun_zero.Rd
Please either omit these examples or export the functions.

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...)
by the user. That's why wrapping examples in \dontrun{} adds the comment 
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec,
or replace \dontrun{} with \donttest{}.

This is a resubmission. In this version I have:  
* Version changed: 1.1.1 -> 1.1.2
* Example data are added in data.
* Example files are added in inst/extdata.
* Examples for unexported functions are removed.
* Examples for exported functions are unwrapped.
* Examples for fun_homology and fun_read_check are wrapped with \dontrun{}.

## Test environments
* macOS 10.15.7, R 3.5.3
* ubuntu 18.04, R 3.6.3
* Windows 10, R 4.0.3
Environments with the following packages already installed:
- BSgenome.Hsapiens.UCSC.hg38  
- BSgenome.Hsapiens.UCSC.hg19  
- BSgenome.Mmusculus.UCSC.mm10  

## R CMD check results
Duration: 4m 32.3s
0 errors ✓ | 0 warnings ✓ | 0 notes ✓
R CMD check succeeded  

Examples with CPU or elapsed time > 5s
                user system elapsed
fun_homology   66.401  2.349  68.770
fun_read_check 20.256  0.792  21.049
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
