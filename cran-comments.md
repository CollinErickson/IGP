This packaged was archived by CRAN after I failed to update it on 10/3/22.
I have updated everything so now it should be good to go back on CRAN.

## Test environments
* local Windows 10 install, R 4.0.5
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)
* R-hub
* macOS builder

## R CMD check results

local Windows 11, R 4.2.2 (1/13/23):

  0 errors ✔ | 0 warnings ✔ | 0 notes ✔


Winbuilder (win_devel and win_release both, 1/13/23, 1/15/23) just has 1 NOTE. 

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

macOS builder (1/13/23), R 4.2.1:
  OK

R-Hub Windows Server (1/14/23):
  OK

R-Hub Ubuntu Linuxand Fedora Linux (1/14/23):
  NOTEs for new submission, archived on CRAN, misspelled "et" and "al",
  and some slow examples.

## Reverse dependencies

There are no reverse dependencies.
