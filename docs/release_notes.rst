Releases Notes
==============

ASYNCH release notes provide information on the features and improvements in each release. This page includes release notes for major releases and minor (bugfix) releases. If you are upgrading from an earlier version of ASYNCH, you will find essential information in the Breaking Changes associated with the relevant release notes.

Version 1.3
-----------

This release is the result of a loong run profiling and optimization work.

Memory footprint improvements, better data structure reduce the memory usage typically by a factor two. If you have memory available you may want to increase the number of buffers, see :ref:`Buffer Sizes`.

Performance improvements, ASYNCH runs about 30% faster.

Version 1.2
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

The global file structure has changed. The time of the simulation is now given in absolute time, see :ref:`Simulation period`.

The snapshots in HDF5 format has changed, see :ref:`Ini HDF5 Files`.

New Features
~~~~~~~~~~~~

Time series outputs can be written in HDF5 format, see :ref:`Time Series Location`.
