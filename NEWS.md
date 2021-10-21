test
# coMethDMR 0.99.9
We have been working on bug fixes and formatting/documentation changes as requested during the Bioconductor review. See these requests here: <https://github.com/Bioconductor/Contributions/issues/2064>

## Breaking Changes

- In the `CpGsInfoAllRegions()`, `CpGsInfoOneRegion()`, and `GetCpGsInRegion()` functions, note that we have added a new argument `region_gr` (in the second position) which allows for the input of a `GRanges` object to these functions. This will cause errors in existing code if that code uses positional argument matching. Any existing code that matches arguments by name will be unaffected.


# coMethDMR 0.99.2

## Major Changes
Bolstered documentation across all package help / manual files

## Bug Fixes


# coMethDMR 0.99.1
Migrated all clean files to this repository.
All development history in https://github.com/TransBioInfoLab/coMethDMR_old



# coMethDMR 0.0.0.9001

## Major Changes
Included EPIC arrays annotation (#1)
Added a check that: (#5)
The rownames of the beta matrix are probe IDs, and
There are only numeric columns in the beta matrix

## Bug Fixes
