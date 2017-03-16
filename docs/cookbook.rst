Cookbook
========

This is not a section about French cuisine, although you will find pretty good recipes here and contribution from our users.

Preprocessing
--------------

Todo.

Postprocessing
--------------

Reading the HDF5 outputs with Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may have to install the ``tables`` package.

.. code-block:: sh

  pip install tables

Reading `PyTables <http://www.pytables.org>`_ is pretty straighforward:

.. code-block:: python

  import tables

  # Open the file in read-only mode
  h5file = tables.open_file("outputs.h5", "r")

  # Print some info about the file
  print(h5file)

  for row in h5file.root.outputs:

      # Read current row content
      # The name of the comuns depends on your outputs settings
      link_id = row['LinkID']
      time = row['Time']
      state_0 = row['State0']

      # Do something with the values
      print "link_id:", link_id, "time:", time, "state_0:", state_0

    # Close the file
    h5file.close()

Here is an other way using the the ``numpy`` and ``h5py`` package:

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
