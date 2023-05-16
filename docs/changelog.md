
# 0.7.0

- Removed `consensus_masking_threshold` option.

- Removed `keep_all_intermediate_files` option.

- Removed `masked_nextclade` option.

- Removed `masked_ivar` option.

- Removed `testbiin` option

- Added unittests for modules `bin.best_reference`, `bin.detect_variants`.

- `runtest` now outputs to the default temporary directory.

- `runtest` now has `compare_TestResults` to validate the observed results with known results.

- Added `gap_open_penalty` option.

- Added `gap_extension_penalty` option.


# 0.6.0

- Added picard to remove duplicates at the alignment level

- Fixed bug where files used for processing in best_reference.py were not removed after script was finished.

- Removed sampleReports directory

- Options using masked consensus sequences are availble but will still used unmasked consensensus sequence. This option may be updated or removed in the future.

# 0.5.0

- Major rewrite to remove R dependency.
