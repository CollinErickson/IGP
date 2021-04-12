This was archived on CRAN on 1/28/21 since it used the package PythonInR.
I switched so it uses reticulate to call Python instead.

It took me longer than expected to get this submitted again because I had
to update one of my other packages and I ran into some Ubuntu issues.

## Test environments
* local Windows 10 install, R 4.0.5
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

local Windows:

    0 errors | 0 warnings | 0 notes

R-hub only has harmless notes:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    New submission
    Package was archived on CRAN
    
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2021-01-28 as requires archived package
    
        'PythonInR'.

Winbuilder (win_devel and win_release both) just has 3 NOTES. The same as R-hub, plus some for slow examples:

    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    New submission
    
    Package was archived on CRAN
    
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2021-01-28 as requires archived package
        'PythonInR'.
        
    
    * checking examples ...
    ** running examples for arch 'i386' ... [25s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                             user system elapsed
    IGP_LOOEC_GauPro_kernel 10.92   0.01   11.75
    ** running examples for arch 'x64' ... [27s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                            user system elapsed
    IGP_LOOEC_GauPro_kernel 11.7   0.02   13.26
  

## Reverse dependencies

There are no reverse dependencies.
