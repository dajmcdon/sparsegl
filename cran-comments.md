# Revision to address CRAN check results

* Simplified Makevars to better match "Writing R Extensions"

## Test environments

R-lib:

- {os: macOS-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}

Rhub:


|   |name           |r_version                                          |os_name                           |
|:--|:--------------|:--------------------------------------------------|:---------------------------------|
|24 |linux          |*                                                  |NA                                |
|22 |macos          |*                                                  |NA                                |
|23 |macos-arm64    |*                                                  |NA                                |
|21 |windows        |*                                                  |NA                                |
|6  |c23            |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|7  |clang-asan     |R Under development (unstable) (2024-06-24 r86823) |Ubuntu 22.04.4 LTS                |
|8  |clang16        |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|9  |clang17        |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|10 |clang18        |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|11 |clang19        |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|18 |donttest       |R Under development (unstable) (2024-06-22 r86814) |Ubuntu 22.04.4 LTS                |
|12 |gcc13          |R Under development (unstable) (2024-06-24 r86823) |Fedora Linux 38 (Container Image) |
|13 |gcc14          |R Under development (unstable) (2024-06-24 r86823) |Fedora Linux 40 (Container Image) |
|16 |mkl            |R Under development (unstable) (2024-06-24 r86823) |Fedora Linux 38 (Container Image) |
|14 |nold           |R Under development (unstable) (2024-06-24 r86823) |Ubuntu 22.04.4 LTS                |
|20 |rchk           |R Under development (unstable) (2024-06-24 r86823) |Ubuntu 22.04.4 LTS                |
|1  |ubuntu-clang   |R Under development (unstable) (2024-06-24 r86823) |Ubuntu 22.04.4 LTS                |
|2  |ubuntu-gcc12   |R Under development (unstable) (2024-06-24 r86823) |Ubuntu 22.04.4 LTS                |
|3  |ubuntu-next    |R version 4.4.1 Patched (2024-06-20 r86823)        |Ubuntu 22.04.4 LTS                |
|4  |ubuntu-release |R version 4.4.1 (2024-06-14)                       |Ubuntu 22.04.4 LTS                |


## R CMD check results

0 errors | 0 warnings | 1 notes

❯ checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

Response: This note does not seem to be related to our R package and can likely be ignored. 
