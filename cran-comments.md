## Resubmission
This is a resubmission. In this version I have:  
* Version changed: 1.1.3 -> 2.0.0
* Many functions were removed and the structure was simplified.
- fun_convert
- fun_hairpin_trimming
- fun_insert_length
- fun_load_id
- fun_load_mutation_gz
- fun_save_gz

## Test environments
* macOS 14.5, R 4.2.3
* ubuntu 18.04, R 3.6.3
* Windows 10, R 4.0.3
Environments with the following packages already installed:
- BSgenome.Hsapiens.UCSC.hg38  
- BSgenome.Hsapiens.UCSC.hg19  
- BSgenome.Mmusculus.UCSC.mm10  

## R CMD check results
Duration: 3m 0.7s
0 errors ✓ | 0 warnings ✓ | 0 notes ✓
R CMD check succeeded  

Examples with CPU or elapsed time > 5s
                user system elapsed
fun_homology   66.401  2.349  68.770
fun_read_check 20.256  0.792  21.049
fun_convert     6.587  0.703   7.306
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
