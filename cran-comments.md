This packaged was archived by CRAN after I failed to update it on 10/3/22.
I have updated everything so now it should be good to go back on CRAN.

## Test environments
* local Windows 10 install, R 4.0.5
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)
* R-hub

## R CMD check results

local Windows 11 R 4.2.2 (12/29/22):

    0 errors | 0 warnings | 0 notes

R-hub Ubuntu 20.04 and R-hub Fedora Linux (12/29/22) have notes for a slow example,
being archived from CRAN, and potentially misspelled words that are fine.
None of these are problems.

R-hub Windows Server (12/29/22) has OK.

Winbuilder (win_devel and win_release both) just has 1 NOTE. 

    * checking CRAN incoming feasibility ... [9s] NOTE
    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    New submission
    
    Package was archived on CRAN
    
    Possibly misspelled words in DESCRIPTION:
      al (15:30)
      et (15:27)
    
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2022-10-03 as check issues were not
        corrected despite reminders.
  

## Reverse dependencies

There are no reverse dependencies.
