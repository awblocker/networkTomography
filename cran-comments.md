## Release summary

This is a minor release to move ipfp from an included function to a dependency.
This change removes all C code from networkTomography.

## Test environments

* local x86_64-pc-linux-gnu R 4.2.0
* Ubuntu 16.04.7 LTS (on travis-ci), R 4.0.2 and dev 2022-02-09 r81690
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

## revdepcheck results

We checked 10 reverse dependencies, comparing R CMD check results
across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

