I submitted this to CRAN on 4/13/21. It was rejected for using dontrun in
examples instead of donttest. I removed some instances of dontrun where
the examples could actually be run. I changed the remaining to donttest.
These are examples that take too long or call Python through reticulate.

In the email from Gregor Sayer, he also asked said I should add any
references to the description file. But there are not any papers to
reference, this is a standalone software project.

I have rerun on all tests on my laptop, R-hub, Travis, and win-devel,
with no new errors, warnings, or notes.

(Below is comment from 4/13)

This was archived on CRAN on 1/28/21 since it used the package PythonInR.
I switched so it uses reticulate to call Python instead.

It took me longer than expected to get this submitted again because I had
to update one of my other packages and I ran into some Ubuntu issues.

I submitted this on 4/12/21, but it was rejected by the auto test
because there was a NOTE for an example that took more than 10 seconds.
I have changed the example to avoid this issue. Now there is only a note
saying that it was archived on CRAN.


## Test environments
* local Windows 10 install, R 4.0.5
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)
* R-hub

## R CMD check results

local Windows:

    0 errors | 0 warnings | 0 notes

R-hub only has one NOTE for all three:
  
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Collin Erickson <collinberickson@gmail.com>’
    
    New submission
    
    Package was archived on CRAN
    
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2021-01-28 as requires archived package
        'PythonInR'.

Winbuilder (win_devel and win_release both) just has 1 NOTE. 

    * package encoding: UTF-8
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    New submission
    
    Package was archived on CRAN
    
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2021-01-28 as requires archived package
        'PythonInR'.
  

## Reverse dependencies

There are no reverse dependencies.
