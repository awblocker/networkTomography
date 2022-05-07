## Release summary

This is a minor release to move ipfp from an included function to a dependency.
This change removes all C code from networkTomography.

## Test environments

* local x86_64-pc-linux-gnu R 4.2.0
* Windows Server 2022, R-devel, 64 bit (via R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (via R-hub)
* Fedora Linux, R-devel, clang, gfortran (via R-hub)
* mac-builder release
* win-builder devel

## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and
dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

