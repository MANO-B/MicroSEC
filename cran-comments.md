## Resubmission
This is a resubmission. In this version I have:  
* Version changed: 1.1.3 -> 2.0.0  
* Many functions were removed and the structure was simplified.  
* Bugs fixed.
- fun_convert  
- fun_hairpin_trimming  
- fun_insert_length  
- fun_load_id  
- fun_load_mutation_gz  
- fun_save_gz  
  
## Test environments  
Local environments  
* macOS-arm64 14.5, R 4.2.3

R-hub environments  
* gcc13 R-devel (2024-06-28 r86850), Fedora Linux 38 (Container Image)  
* linux R-4.2.3, ubuntu-latest on GitHub  
* macos R-4.2.3, macos-13 on GitHub  
  
## R CMD check results  
Duration: 4m 9.5s  
0 errors ✔ | 0 warnings ✔ | 0 notes ✔  
R CMD check succeeded  
  
Examples with CPU (user + system) or elapsed time > 5s  
                user system elapsed  
fun_homology   76.889  0.611  77.783  
fun_read_check 15.179  0.461  15.847  
  
## Downstream dependencies  
There are currently no downstream dependencies for this package.  
