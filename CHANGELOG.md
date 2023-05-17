# Changelog

## 0.8.1 - 2023-05-17

This is a major release.

### Added

- Integrated congruency of SNVs within paired-end reads into contamination calculations (a1744a4).
- Base cutoff values are now dynamically calculated based upon gene-specific quality score, length, and depth of coverage, with a starting cutoff of 3 which can be changed using `--base_cutoff` (880445d, 84b7d91, c06d438).

### Changed

- Refactored code by moving methods to `methods.py` (a7e9af6).
- Improved `README.md` and MkDocs documentation with increased accuracy and readability.
- Pytest tests now use downsampled samples from the originally published ConFindr benchmarking dataset, and instructions for running these tests have been added to the MkDocs (cc27c94).
- Enforced Phred33 encoding for `bbduk.sh` calls to support future development with Nanopore reads (#39) (ac3b976).

### Removed

- Percentage contamination reporting; this was found to be unreliable and sometimes misleading (ec3ae7a).
- `--cross_details` flag; analysis is now always continued after cross-genus contamination has been detected (ec3ae7a).

### Fixed

- TypeError that occurred when using an older version of BioPython (#27, #30, #38, #41) (19d0d1d, 96e1c7d).
- Error in `install.md` which suggested that rMLST databases are freely available to all users (a1ce7dc).