Result files, plots and NeXus export
====================================

Structured SPR-KKR result files
-------------------------------

ASE2SPRKKR includes parsers for several structured scientific
SPR-KKR outputs (such as DOS, BSF, SPC, JXC, RAT and SFN).
Plotting requires the optional ``plot`` dependencies; install them
if they are not already available::

   python -m pip install 'ase2sprkkr[plot]'

Then plot a recognized file::

   # Plot the result values.
   ase2sprkkr plot calculation.dos

   # Inspect which result-file formats are supported.
   ase2sprkkr info --output-file

   # Show the names of the values available in this file.
   ase2sprkkr plot calculation.dos --show-values

Available values depend on the format and installed package version.

Load a recognized structured file using ``ResultFile``:

.. code-block:: python

   from ase2sprkkr import ResultFile

   data = ResultFile.from_file("calculation.dos")

``ResultFile`` is the current public name; ``OutputFile`` is a
deprecated compatibility name. Result-file parsing is independent of
running a new calculation.

Grace .agr projects
-------------------

Grace ``.agr`` projects are supported separately from the structured
SPR-KKR file registry:

.. code-block:: python

   from ase2sprkkr.result_files import read_agr, read_plottable_file

   grace = read_agr("bands.agr")
   print(grace.datasets)
   data = read_plottable_file("bands.agr")

The latter helper dispatches by filename and handles both Grace
projects and recognized SPR-KKR result files.

To plot a Grace file or save a figure::

   ase2sprkkr plot bands.agr
   ase2sprkkr plot bands.agr -o bands.png

To list the values in a Grace file, use::

   ase2sprkkr plot bands.agr --show-values

The ``--value`` option selects a value to plot and can be repeated. SFN shape-function files additionally support
``--mesh``, ``--cell``, ``--periodic``, ``--shape``,
``--slice`` and ``--radial`` plotting modes.

ARPES SPC to NeXus/HDF5
-----------------------

The 4.0 series can export supported **ARPES .spc spectra** in NeXus
format. Install the separate optional dependencies::

   python -m pip install 'ase2sprkkr[nexus]'

Export from the command line::

   ase2sprkkr nexus arpes.spc
   ase2sprkkr nexus arpes.spc -o arpes.nxs -i arpes.inp

By default, the output has the same basename as the source and a
``.nxs`` suffix. The optional ``-i`` argument supplies the native
SPR-KKR input file for provenance; an unambiguous accompanying input
file can also be found automatically. The exporter retains available
axes, spectra and source settings.

From Python, use the same writer on a compatible SPC result:

.. code-block:: python

   from ase2sprkkr import ResultFile

   spectrum = ResultFile.from_file("arpes.spc", try_only="spc")
   spectrum.to_nexus("arpes.nxs", input_file="arpes.inp")

Missing experimental geometry or calibration cannot be reconstructed
from the SPC file. Intensities may therefore be labeled as arbitrary
units; this does **not** mean they have been normalized or calibrated.

For API help and the CLI syntax, see :doc:`cli`.
