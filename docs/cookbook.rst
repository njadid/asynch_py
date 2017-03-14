Cookbook
========

This is not a section about French cuisine, although you will find pretty good recipes here and contribution from our users.

Preprocessing
--------------



Postprocessing
--------------

Reading the HDF5 outputs with python with h5py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

  import numpy as np
  import h5py

  # open file and iterate over each row
  with h5py.File(h5_file_path, "r") as hdf_file:
      hdf_file_content = np.array(hdf_file.get("outputs"))
      for i in range(len(hdf_file_content)):
          # read current row content
          cur_row = hdf_file_content[i]
          # get the values of each columns
          cur_linkid = cur_row[0]
          cur_time = cur_row[1]
          cur_state0 = float(cur_row[2])
