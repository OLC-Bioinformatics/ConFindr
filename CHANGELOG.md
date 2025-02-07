# Changelog

## 0.8.2 - 2025-02-07

This is a minor release. The most important changes include updates to URL schemes, improvements in handling temporary directories, and enhancements in read name processing.

### Added

- Added checks to ensure temporary directories exist before attempting to remove them in `confindr_src/confindr.py` (#65).
- Initialized dictionaries for summarizing base types with specific keys in `confindr_src/methods.py` to prevent key errors (#65).
- Added Pytest unit tests and test FASTQ files for ensuring proper parsing of paired-end Illumina FASTQ headers (#70).

### Changed

- Improved logic for using temporary directories during database creation in `confindr_src/methods.py` (#65). Fixes #57.
- Updated the installation guide to reflect changes in the PubMLST API key generation process in `docs/install.md` (#62, #71).
- Corrected the URL for the example dataset in `docs/usage.md` to point to the latest version.

### Fixed

- Updated URLs from HTTP to HTTPS in the `__init__` method of `confindr_src/database_setup.py` to ensure secure connections and allow for successful downloads of databases (#69). Fixes #67.
- Improved addition of read direction suffixes to read names in `confindr_src/methods.py` to handle cases where the suffix already exists (#70). Fixes #63 and resolves #54.

## 0.8.1 - 2023-05-19

This is a major release.

### Added

- Integrated congruency of SNVs within paired-end reads into contamination calculations (a1744a4).
- Base cutoff values are now dynamically calculated based upon gene-specific quality score, length, and depth of coverage, with a starting cutoff of 3 which can be changed using `--base_cutoff` (880445d, 84b7d91, c06d438).
- Option to download rMLST databases using `-u/--unverified` within the `confindr_database_setup` command, for downloading databases behind a firewall and/or have a self-signed certificate.

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
