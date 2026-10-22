# Compact result-file smoke fixtures

These small examples are derived from SPR-KKR result files in ../../examples,
not intended to replace the full-resolution numerical reference data.

- SnTe BSF / BSF-SPOL: the original EK-REL data had NE=10, NK=15,
  NQ_eff=4. Selected energy indices [0,3,6,9] and k indices [0,3,6,9,14],
  updating both NE headers, both NK declarations, INDKDIR, the K grid,
  and the interleaved spin/channel RAW_DATA blocks in their native ordering.
  The source EMIN/EMAX and nonuniform interior samples are retained.
- Mo2Se4 energy SPC: seven selected rows from original 300 energies,
  with header NE=7, NT=NP=1 and all eight native columns retained.
- Mo2Se4 kx_ky SPC: synthetic orthogonal 5x7 mesh derived from the
  Mo2Se4 intensity columns, using seven native columns and NP=5/NT=7.

The full files remain available for integration, high-resolution visualization,
quantitative regression and physically faithful spectra. These smoke tests only
check parser shape contracts and finite values.

## Regenerating smaller fixtures

From the repository checkout, run the utility on the original examples:

```sh
python src/ase2sprkkr/result_files/test/compact_spectral_data.py \
  src/ase2sprkkr/result_files/examples/FePt.new_BSFEK_BLOCHSF.bsf \
  src/ase2sprkkr/result_files/test/data/FePt_BSFEK_small.bsf \
  --energies 8 --kpoints 12

python src/ase2sprkkr/result_files/test/compact_spectral_data.py \
  src/ase2sprkkr/result_files/examples/WSe2_ARPES.spc \
  src/ase2sprkkr/result_files/test/data/WSe2_ARPES_small.spc \
  --energies 8 --kpoints 12
```

The utility validates each native data block before selecting its subgrid.
Generated files should be parsed by `ResultFile.from_file`, their declared
dimensions compared to the actual data shape, and the plotting smoke test run
before promoting them as replacements. Keep original large datasets for
numerical comparisons requiring exact original dimensions, energies or pixels.
