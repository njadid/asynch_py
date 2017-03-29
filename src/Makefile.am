noinst_LIBRARIES = libasynch.a
libasynch_a_SOURCES = \
  advance.c \
  asynch_interface.c \
  blas.c \
  comm.c \
  compression.c \
  config_gbl.c \
  data_types.c \
  date_manip.c \
  db.c \
  forcings.c forcings_io.c \
  io.c \
  misc.c \
  outputs.c \
  partition.c \
  processdata.c \
  riversys.c \
  rksteppers.c \
  sort.c \
  system.c \
  models/check_consistency.c \
  models/check_state.c \
  models/definitions.c \
  models/equations.c \
  solvers/dopri5_dense.c \
  solvers/lagrange.c \
  solvers/radau.c \
  solvers/rk3_2_dense.c \
  solvers/rk4_3_dense.c \
  steppers/explicit.c \
  steppers/explicit_index1.c \
  steppers/explicit_index1_dam.c \
  steppers/forced.c

bin_PROGRAMS = asynch
asynch_SOURCES = optparse.c asynchdist.c
asynch_LDADD = libasynch.a $(HDF5_LIBS) $(POSTGRESQL_LIBS)
asynch_LDFLAGS = $(HDF5_LDFLAGS) $(POSTGRESQL_LDFLAGS)

include_HEADERS = structs.h structs_fwd.h asynch_interface.h

AM_CFLAGS = $(HDF5_CPPFLAGS) $(POSTGRESQL_CPPFLAGS)
