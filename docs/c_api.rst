C API
=====

ASYNCH comes with a collection of routines to access specific details of the underlying solver structure. This allows a user to construct his or her own programs, while making calls to the solver in particular ways. This is useful for modifying how the solvers behave, and for producing more specialized outputs. A user can also specify his or her own models in a separate module. Currently, interface routines exist for the C/C++ and Python programming languages.

Typical Interface Usage
-----------------------

In this section, we provide an overview of how to use the API routines for running simulations. Regardless of the language, the general order of calling API routines is the same. Further routines may be added to increase flexibility.

A basic program that uses the ASYNCH solvers to perform calculations and takes advantage of all the features of a global file is the following (written in C):

.. code-block:: c

  //Make sure we finalize MPI
  void asynch_onexit(void)
  {
      int flag;
      MPI_Finalized(&flag);
      if (!flag)
          MPI_Finalize();
  }

  int main(int argc,char* argv[])
  {
    int res;
    char *global_filename = argv[1];

    //Initialize MPI stuff
    res = MPI_Init(&argc, &argv);
    if (res == MPI_SUCCESS)
        atexit(asynch_onexit);
    else
    {
        print_err("Failed to initialize MPI");
        exit(EXIT_FAILURE);
    }

    //Init asynch object and the river network
    AsynchSolver asynch;
    Asynch_Init(&asynch, MPI_COMM_WORLD);
    Asynch_Parse_GBL(&asynch, global_filename);
    Asynch_Load_Network(&asynch);
    Asynch_Partition_Network(&asynch);
    Asynch_Load_Network_Parameters(&asynch, 0);
    Asynch_Load_Dams(&asynch);
    Asynch_Load_Numerical_Error_Data(&asynch);
    Asynch_Initialize_Model(&asynch);
    Asynch_Load_Initial_Conditions(&asynch);
    Asynch_Load_Forcings(&asynch);
    Asynch_Load_Save_Lists(&asynch);
    Asynch_Finalize_Network(&asynch);
    Asynch_Calculate_Step_Sizes(&asynch);

    //Prepare output files
    Asynch_Prepare_Temp_Files(&asynch);
    Asynch_Write_Current_Step(&asynch);
    Asynch_Prepare_Peakflow_Output(&asynch);
    Asynch_Prepare_Output(&asynch);

    //Perform the calculations
    Asynch_Advance(&asynch, 1);

    //Create output files
    Asynch_Take_System_Snapshot(&asynch, NULL);
    Asynch_Create_Output(&asynch, NULL);
    Asynch_Create_Peakflows_Output(&asynch);

    //Clean up
    Asynch_Delete_Temporary_Files(&asynch);
    Asynch_Free(&asynch);

    return EXIT_SUCCESS;
  }

Details of each function call can be found in Section [sec: user interface routines].

The program begins by initializing the asynchsolver object. The MPI communicator consisting of all processes (:code:`MPI_COMM_WORLD`) is used for the calculations. Next, the global file specified as a command line argument to the program is parsed. Based upon the information specified by the global file, the different components of the network are constructed.

Next, the outputs for the program are initialized. The files for holding calculation results while the program runs is initialized (i.e. the temporary files). The initial values of the states are written to disk. Any initializations needed for the output peakflow and output time series sources are done.

With a call to :code:`Asynch_Advance`, the calculations are performed. Any time series results are temporarily stored in the temporary files.

Once the calculations are complete, any necessary outputs are created (snapshot, output time series, output peakflow data). This could include writing to files or database tables, depending upon the options selected in the global file.

Lastly, clean up routines are called. The temporary files are deleted. The asynchsolver object is also deleted from memory with a call to . Note that this routine does not need to be called if using the interface functions from a language that supports automatic garbage collection.

Of course, outputting more information might be useful (timing results, command line parameter checking, results printed to screen, etc), and additional features may be needed (outputting data to multiple sources, creating peakflow data over intervals of time, etc). However, this is the basic structure needed to perform simulations. The source files *asynchdist.c* and *asynchdist.py* are essentially duplicates of the above program, but with information printed to screen.

User Interface Routines
-----------------------

In this section, routines for operating the solver are described. These routines can be used to create an instance of an ASYNCH solver, and manipulate properties such as total simulation time, when data output occurs, etc. Creation of custom outputs is discussed in Section [sec: custom outputs] and creation of custom models is discussed in Section [sec: custom models].

Solver Initialization
~~~~~~~~~~~~~~~~~~~~~

C Interface:

.. code-block:: c

  void Asynch_Init(AsynchSolver* asynch, MPI_Comm comm)

-  Inputs

-  MPI_Comm comm: The MPI communicator to use with this solver object.

-  int\* argc: Pointer to the integer holding the number of command line arguments.

-  char\*\* argv[]: Pointer to the command line arguments.

-  Return Value

-  A pointer to a newly created asynchsolver object. If an error occurred in creating the object, the return value is *NULL*.

This routine creates an instance of an ASYNCH solver. Multiple solvers could be created if multiple problems are to be solved. An instance of an *asynchsolver* contains all the relevant data structures and information to solve the equations for an underlying model. An *asynchsolver* object should be destroyed with a call to *Asynch_Free*. The inputs *argc* and *argv* are only passed on to *MPI_Init*, and can be *NULL*.

::

Python Interface:
asynchsolver()

-  Return Value

-  A reference to a newly created *asynchsolver* object.

This routine creates and initializes an ASYNCH solver, and is similar to the corresponding C initializer. The interface functions described below are the only members of the asynchsolver object accessible. An asynchobject is automatically destroyed when leaving scope.

Asynch_Free
~~~~~~~~~~~~

::

C Interface:
void Asynch_Free(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: A pointer to the *asynchsolver* object to free.

This routine deallocates the memory occupied by an *asynchsolver* object created with a call to *Asynch_Init*. The routine is exclusive to the C interface.

Asynch_Parse_GBL
~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Parse_GBL(asynchsolver* asynch,
char* gbl_filename)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* gbl_filename: Filename of a global file.

::

Python Interface:
Parse_GBL(gbl_filename)

-  Inputs

-  gbl_filename: String with the filename of a global file.

This routine opens and processes a global file. It reads all specified database connection files, but does not process any other input file. An error in this routine is considered fatal, and results in a call to the routine *MPI_Abort* on the communicator used to create asynch.

Asynch_Load_Network
~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Network(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Network()

This routine processes topology inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This initializes each Link object and sets their parent and child information. Generally, this is the first initialization routine to call after parsing a global.

Asynch_Partition_Network
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Partition_Network(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Partition_Network()

This routine assigns the Links of the *asynchsolver* object to the MPI processes. This routine must be called after *Asynch_Load_Network*.

Asynch_Load_Network_Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Network_Parameters(
asynchsolver* asynch,short int load_all)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  short int load_all: Flag to load parameters at every Link.

::

Python Interface:
Load_Network_Parameters(load_all)

This routine processes parameter inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. Setting *load_all* to true causes every MPI process to store the parameters at every Link. Setting *load_all* to false causes the MPI processes to only store parameters for Links assigned to them. This routine can be called before *Asynch_Partition_Network* only if *load_all* is set to true.

Asynch_Load_Dams
~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Dams(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Dams()

This routine processes the dam inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Network* and *Asynch_Load_Network_Parameters* have been called.

Asynch_Load_Numerical_Error_Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Numerical_Error_Data(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Numerical_Error_Data()

This routine processes numerical solver inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Network* has been called.

Asynch_Initialize_Model
~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Initialize_Model(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Initialize_Model()

This routine sets the model specific routines for each link for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Network* and *Asynch_Load_Network_Parameters* have been called.

Asynch_Load_Initial_Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Initial_Conditions(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Initial_Conditions()

This routine processes the initial condition inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Network* and *Asynch_Initialize_Model* have been called.

Asynch_Load_Forcings
~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Forcings(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Forcings()

This routine processes the forcing inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Parameters* has been called.

Asynch_Load_Save_Lists
~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Load_Save_Lists(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Load_Save_Lists()

This routine processes save list inputs for the *asynchsolver* object as set in the global file read by *Asynch_Parse_GBL*. This routine should be called after *Asynch_Partition_Network* has been called.

Asynch_Finalize_Newtork
~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Finalize_Network(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Finalize_Network()

This routine checks that all inputs are loaded for the *asynchsolver* object. Some small final initializations are also performed. This routine should be called as the last initialization routine.

Asynch_Calculate_Step_Sizes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Calculate_Step_Sizes(asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Calculate_Step_Sizes()

This routine processes calculates appropriate step sizes for the integrators at each link in the *asynchsolver* object. This routine must be called before a call to *Asynch_Advance*, and after all initializations are performed.

Asynch_Advance
~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Advance(asynchsolver* asynch,
short int print_flag)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  short int print_flag: If 0, no time series information is produced. Otherwise, time series information is produced.

::

Python Interface:
Advance(print_flag)

-  Inputs

-  print_flag: If 0, no time series information is produced. Otherwise, time series information is produced.

This routine advances the numerical solver up to the time set in *maxtime*. See Section [sec: model type and maxtime]. Calculations to solve the model differential and algebraic equations are performed, using forcing data as needed. If *print_flag* is set, requested output time series are produced internally, but not written to a final output source.

Asynch_Take_System_Snapshot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Take_System_Snapshot(
asynchsolver* asynch,char* name)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* name: String to append to the snapshot output filename.

-  Return value

-  int: An error code. Returns 0 if a snapshot was made, 1 if an error was encountered, -1 if no snapshot is made.

::

Python Interface:
Take_System_Snapshot(name)

-  Inputs

-  name: String to append to the snapshot output filename.

-  Return value

-  Returns 0 if a snapshot was made, 1 if an error was encountered, -1 if no snapshot is made.

This routine creates snapshot output data to either a recovery file or a database table. The current value of every state at every link is outputted. If a recovery file was the format specified in the global file, then the string *name* is appended to the end of the recovery filename. This appending does not occur if *name* is *NULL* or if a database table is the selected format.

For the Python routine *Take_System_Snapshot*, a value of *None* for *name* causes no appending to the filename.

If the return value is -1, then no snapshot output has been selected (i.e. the snapshot flag is 0. See Section [sec: snapshot info].).

Asynch_Set_Database_Connection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_Database_Connection(
asynchsolver* asynch,char* database_info,
unsigned int conn_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* database_info: Information to connect to a database.

-  unsigned int conn_idx: The database connection to set.

::

Python Interface:
Set_Database_Connection(database_info,conn_idx)

-  Inputs

-  database_info: String of information to connect to a database.

-  conn_idx: The database connection to set.

This routine sets a new database for an input or output. If information for a database has already been set, it is released, and the new connection information is set. Database information includes hostname, username, password, etc. This is the same information that is available in the header of database connection files (see Section [sec: database connection files]). Several database connections exist for every *asynchsolver* object. *conn_idx* can take the values:

-  ASYNCH_DB_LOC_TOPO

-  ASYNCH_DB_LOC_PARAMS

-  ASYNCH_DB_LOC_INIT

-  ASYNCH_DB_LOC_RSV

-  ASYNCH_DB_LOC_HYDROSAVE

-  ASYNCH_DB_LOC_PEAKSAVE

-  ASYNCH_DB_LOC_HYDRO_OUTPUT

-  ASYNCH_DB_LOC_PEAK_OUTPUT

-  ASYNCH_DB_LOC_SNAPSHOT_OUTPUT

-  ASYNCH_DB_LOC_FORCING_START

The last value for *conn_idx* is the database connection for the first forcing specified in the global file. Database connections for other forcings can be access sequentially. For example, to access the forcing with index 2 (the third forcing in a global file), set *conn_idx* as

ASYNCH_DB_LOC_FORCING_START + 2

Asynch_Get_Total_Simulation_Time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
double Asynch_Get_Total_Simulation_Time(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  double: The current end value of the simulation time of *asynch*.

::

Python Interface:
Get_Total_Simulation_Time()

-  Return value

-  The current end value of the simulation time of the *asynchsolver* object.

This routine returns the value of *maxtime*, as defined in Section [sec: model type and maxtime].

Asynch_Set_Total_Simulation_Time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_Total_Simulation_Time(
asynchsolver* asynch,double new_time)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double new_time: The value to set the maximum simulation time.

::

Python Interface:
Set_Total_Simulation_Time(new_time)

-  Inputs

-  new_time: The value to set the maximum simulation time.

This routine sets the value of *maxtime*, as defined in Section [sec: model type and maxtime], to the value *new_time*.

Asynch_Get_Last_Rainfall_Timestamp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
unsigned int Asynch_Get_Last_Rainfall_Timestamp(
asynchsolver* asynch,unsigned int forcing_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int forcing_idx: An index of a forcing.

-  Return Value

-  unsigned int: The timestamp of the last timestamp for the forcing with index *forcing_idx*.

::

Python Interface:
Get_Last_Rainfall_Timestamp(forcing_idx)

-  Inputs

-  forcing_idx: An index of a forcing.

-  Return Value

-  The timestamp of the last timestamp for the forcing with index *forcing_idx*.

This routine returns the last timestamp for a forcing. This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files, or database table. See Section [sec: forcing inputs] for a description of these formats.

Asynch_Set_Last_Rainfall_Timestamp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_Last_Rainfall_Timestamp(
asynchsolver* asynch,
unsigned int epoch_timestamp,
unsigned int forcing_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int epoch_timestamp: The value to set the last rainfall timestamp.

-  unsigned int forcing_idx: An index of a forcing.

::

Python Interface:
Set_Last_Rainfall_Timestamp(
epoch_timestamp,forcing_idx)

-  Inputs

-  epoch_timestamp: The value to set the last rainfall timestamp.

-  forcing_idx: An index of a forcing.

This routine sets the last timestamp for a forcing. This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files, or database table. See Section [sec: forcing inputs] for a description of these formats.

Asynch_Get_First_Rainfall_Timestamp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
unsigned int Asynch_Get_First_Rainfall_Timestamp(
asynchsolver* asynch,unsigned int forcing_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int forcing_idx: An index of a forcing.

-  Return Value

-  unsigned int: The timestamp of the first timestamp for the forcing with index *forcing_idx*.

::

Python Interface:
Get_First_Rainfall_Timestamp(forcing_idx)

-  Inputs

-  forcing_idx: An index of a forcing.

-  Return Value

-  The timestamp of the first timestamp for the forcing with index *forcing_idx*.

This routine returns the first timestamp for a forcing. This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files, or database table. See Section [sec: forcing inputs] for a description of these formats.

Asynch_Set_First_Rainfall_Timestamp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_First_Rainfall_Timestamp(
asynchsolver* asynch,
unsigned int epoch_timestamp,
unsigned int forcing_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int epoch_timestamp: The value to set first rainfall timestamp.

-  unsigned int forcing_idx: An index of a forcing.

::

Python Interface:
Set_First_Rainfall_Timestamp(
epoch_timestamp,forcing_idx)

-  Inputs

-  epoch_timestamp: The value to set first rainfall timestamp.

-  forcing_idx: An index of a forcing.

This routine sets the first timestamp for a forcing. This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files, or database table. See Section [sec: forcing inputs] for a description of these formats.

Asynch_Set_RainDB_Starttime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_RainDB_Starttime(
asynchsolver* asynch,
unsigned int epoch_timestamp,
unsigned int forcing_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int epoch_timestamp: The value to set the start time for a forcing.

-  unsigned int forcing_idx: An index of a forcing.

::

Python Interface:
Set_RainDB_Starttime(epoch_timestamp,forcing_idx)

-  Inputs

-  epoch_timestamp: The value to set the start time for a forcing.

-  forcing_idx: An index of a forcing.

This routine sets the start time used for a forcing. This value is used for converting between timestamps in a database table and the local time of the solvers. This can only be used if the forcing with index *forcing_idx* is using a format of database table. See Section [sec: forcing inputs] for a description of this format.

Asynch_Set_Init_File
~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_Init_File(
asynchsolver* asynch,
char* filename)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* filename: Filename to use for initial value data.

::

Python Interface:
Set_Init_File(filename)

-  Inputs

-  filename: Filename to use for initial value data.

This routine sets a file for reading initial value data. The *init_flag* is set based upon the extension of *filename*. The initial data is **NOT** read while executing this routine. A call to *Asynch_Load_System* is needed to set the filename. This routine cannot be used to set the format to a database connection.

Asynch_Prepare_Output
~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Prepare_Output(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Prepare_Output()

This routine prepares the output sources for time series data. Preparation includes creating files or database tables. This routine must be called before any time series data can be produced.

Asynch_Prepare_Peakflow_Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Prepare_Peakflow_Output(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Prepare_Peakflow_Output()

This routine prepares the output sources for the peakflow data. Preparation includes creating files or database tables. This routine must be called before any peakflow data can be produced.

Asynch_Prepare_Temp_Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Prepare_Temp_Files(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Prepare_Temp_Files()

This routine prepares the temporary files for time series data. Preparation includes creating files or database tables. This routine must be called before any time series data can be calculated. A call to *Asynch_Advance* with the *print_flag* set before preparing temporary files will create an error.

Asynch_Write_Current_Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Write_Current_Step(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  int: Returns 0 if the step is written, 1 if no temporary file is available.

::

Python Interface:
Write_Current_Step()

-  Return value

-  Returns 0 if the step is written, 1 if no temporary file is available.

This routine writes the current state of each link to temporary files for every link where a time series output has been specified. Normally, calling *Asynch_Advance* with the *print_flag* set is enough to write output time series. However, advancing an *asynchsolver* object without the *print_flag* set or calls that modify how data is outputted to the temporary files will cause data to not be written for some times. This routine can be called to commit missing data.

Asynch_Create_Output
~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Create_Output(
asynchsolver* asynch,
char* additional_out)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* additional_out: String appended to the end of any output filename.

-  Return value

-  int: Returns 0 if output is written, -1 means there is no data to output.

::

Python Interface:
Create_Output(additional_out)

-  Inputs

-  additional_out: String appended to the end of any output filename.

-  Return value

-  Returns 0 if output is written, -1 means there is no data to output.

This routine takes all data written to temporary files and moves them to a final output destination for time series data. If the output format is data file or CSV file, then the string *additional_out* is appended to the filename. If *additional_out* is *NULL*, then no appending to the filename occurs.

For the Python routine *Create_Output*, no appending to the filename occurs if *additional_out* is *None*.

Asynch_Create_Peakflows_Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Create_Peakflows_Output(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  int: Returns 0 if output is written, -1 means there is no data to output.

::

Python Interface:
Create_Peakflows_Output()

-  Return value

-  Returns 0 if output is written, -1 means there is no data to output.

This routine takes all calculated peakflow data and writes them to a final output destination.

Asynch_Get_Number_Links
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
unsigned int Asynch_Get_Number_Links(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  unsigned int: Number of links in the network.

::

Python Interface:
Get_Number_Links()

-  Return value

-  Number of links in the network.

This routine returns the total number of links in the network of *asynch*.

Asynch_Get_Local_Number_Links
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
unsigned int Asynch_Get_Local_Number_Links(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  unsigned int: Number of links assigned to the current MPI process.

::

Python Interface:
Get_Local_Number_Links()

-  Return value

-  Number of links assigned to the current MPI process.

This routine returns the total number of links assigned to the current MPI process in the network of *asynch*. This number represents the total number of links whose differential and algebraic equations are solved by the current process.

Asynch_Set_Temp_Files
~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Temp_Files(
asynchsolver* asynch,
double set_time,
void* set_value,
unsigned int output_idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double set_time: The local time corresponding the step where the temp files will be set.

-  void\* set_value: The value to which the tempfiles will be set.

-  unsigned int output_idx: The index of the series output of which *set_value* corresponds.

-  Return value

-  int: Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.

::

Python Interface:
Set_Temp_Files(set_time,set_value,output_idx)

-  Inputs

-  set_time: The local time corresponding the step where the temp files will be set.

-  set_value: The value to which the tempfiles will be set.

-  output_idx: The index of the series output of which *set_value* corresponds.

-  Return value

-  Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.

This routine moves the temporary file pointer to a previous value. All future values are deleted. The point where the pointer is moved is determined by *set_value*, which is in the output series determined by *output_idx*. The local time where the next data will be written for each link is *set_time*.

A warning occurs if *set_value* is not found for a link, and that link’s tempfile pointer is not changed.

Asynch_Reset_Temp_Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Reset_Temp_Files(
asynchsolver* asynch,
double set_time)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double set_time: The local time corresponding the step where the temp files will be set.

-  Return value

-  int: Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.

::

Python Interface:
Reset_Temp_Files(set_time)

-  Inputs

-  set_time: The local time corresponding the step where the temp files will be set.

-  Return value

-  Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.

This routine moves the tempfile pointer to the beginning of the file for each link. All future values are deleted. The local time where the next data will be written for each link is *set_time*.

A warning occurs if *set_value* is not found for a link, and that link’s tempfile pointer is not changed.

Asynch_Check_Output
~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Check_Output(
asynchsolver* asynch,
char* name)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* name: The name of the time series output to check.

-  Return value

-  int: Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.

::

Python Interface:
Check_Output(name)

-  Inputs

-  name: The name of the time series output to check.

-  Return value

-  Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.

This routine determines whether the time series with name *name* is being used for the current simulations. If *name* is requested in the global file and has been set, the routine returns 1. If the *name* is requested in the global file, but has not yet been set, this routine returns 0. If *name* was not requested in the global file, the return value is -1. An output can be set with a call to *Asynch_Set_Output*. See Section [sec: asynch\ :sub:`s`\ et\ :sub:`o`\ utput] for details about setting outputs.

Asynch_Check_Peakflow_Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Check_Peakflow_Output(
asynchsolver* asynch,
char* name)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* name: The name of the peakflow output to check.

-  Return value

-  int: Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.

::

Python Interface:
Check_Peakflow_Output(name)

-  Inputs

-  name: The name of the peakflow output to check.

-  Return value

-  Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.

This routine determines whether the peakflow output with name *name* is being used for the current simulations. If *name* is requested in the global file and has been set, the routine returns 1. If the *name* is requested in the global file, but has not yet been set, this routine returns 0. If *name* was not requested in the global file, the return value is -1. A peakflow output can be set with a call to *Asynch_Set_Peakflow_Output*. See Section [sec: asynch\ :sub:`s`\ et\ :sub:`p`\ eakflow\ :sub:`o`\ utput] for details about setting peakflow outputs.

Asynch_Delete_Temporary_Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Delete_Temporary_Files(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Return value

-  int: Returns 0 if the tempfiles were deleted, 1 if no tempfiles exist, and 2 if an error occurred.

::

Python Interface:
Delete_Temporary_Files()

-  Return value

-  Returns 0 if the tempfiles were deleted, 1 if no tempfiles exist, and 2 if an error occurred.

This routine deletes the tempfiles for *asynch*. This is useful for cleaning up temporary files at the end of a simulation, or for deleting the files if they must be reconstructed.

Asynch_Activate_Forcing
~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Activate_Forcing(
asynchsolver* asynch,
unsigned int idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int idx: The index of the forcing to activate.

-  Return value

-  int: Returns 0 if the forcing was activated, 1 if an error occurred.

::

Python Interface:
Activate_Forcing(idx)

-  Inputs

-  idx: The index of the forcing to activate.

-  Return value

-  Returns 0 if the forcing was activated, 1 if an error occurred.

This routine activates a forcing for use. This means when a call to the routine *Asynch_Advance* is made, the forcing with index *idx* will be applied for calculations. By default, all forcings set in a global file are initially active. This routine has the opposite effect of *Asynch_Deactivate_Forcing*.

Asynch_Deactivate_Forcing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Deactivate_Forcing(
asynchsolver* asynch,
unsigned int idx)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int idx: The index of the forcing to deactivate.

-  Return value

-  int: Returns 0 if the forcing was deactivated, 1 if an error occurred.

::

Python Interface:
Deactivate_Forcing(idx)

-  Inputs

-  idx: The index of the forcing to deactivate.

-  Return value

-  Returns 0 if the forcing was deactivated, 1 if an error occurred.

This routine deactivates a forcing for use. This means when a call to the routine *Asynch_Advance* is made, all values for the forcing with index *idx* will be taken as 0. By default, all forcings set in a global file are initially active. This routine has the opposite effect of *Asynch_Activate_Forcing*.

Asynch_Set_Forcing_State
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Forcing_State(
asynchsolver* asynch,
unsigned int idx,
double maxtime,
unsigned int first_file,
unsigned int last_file)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  unsigned int idx: The index of the forcing to modify.

-  double maxtime: The new value of *maxtime* for the forcing.

-  unsigned int first_file: The new value of *first_file* for the forcing.

-  unsigned int last_file: The new value of *last_file* for the forcing.

-  Return value

-  int: Returns 0 if the forcing was deactivated, 1 if an error occurred.

::

Python Interface:
Set_Forcing_State(idx,t_0,first_file,last_file)

-  Inputs

-  int idx: The index of the forcing to modify.

-  maxtime: The new value of *maxtime* for the forcing.

-  first_file: The new value of *first_file* for the forcing.

-  last_file: The new value of *last_file* for the forcing.

-  Return value

-  Returns 0 if the forcing was deactivated, 1 if an error occurred.

This routine sets information for a forcing with index *idx*. This routine can only set information if the forcing is of format binary files, gz binary files, or database table. The database start time is set to *first_file*.

Asynch_Reset_Peakflow_Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Reset_Peakflow_Data(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

::

Python Interface:
Reset_Peakflow_Data()

This routine clears all peakflow data for each link. The time to peak is set to the current local time of *asynch*, and the value of the peakflow is set to the current state vector of the link.

Asynch_Set_System_State
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
void Asynch_Set_System_State(
asynchsolver* asynch,
double t_0,
VEC** values)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double t_0: The local time to set at each link.

-  VEC\*\* values: The vectors to apply to the current state of each link.

::

Python Interface:
Set_System_State(t_0,values)

-  Inputs

-  t_0: The local time to set at each link.

-  values: The vectors to apply to the current state of each link.

This routine sets the current time of each link to *t_0* and the current state to the vector in the array *values*. The vectors are assumed to be in the same order as the links provided by the topology source specified in the global file (see Section [sec: topology]).

Asynch_Set_Peakflow_Output_Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Peakflow_Output_Name(
asynchsolver* asynch,
char* peakflowname)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* peakflowname: Name to set the peakflow output filename.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

::

Python Interface:
Set_Peakflow_Output_Name(peakflowname)

-  Inputs

-  peakflowname: Name to set the peakflow output filename.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

This routine sets the filename of the output peakflow file. If a database connection is used, then no changes are made.

Asynch_Get_Peakflow_Output_Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Get_Peakflow_Output_Name(
asynchsolver* asynch,
char* peakflowname)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* peakflowname: Name of the peakflow output filename returned.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

::

Python Interface:
Get_Peakflow_Output_Name()

-  Outputs

-  Returns a list with two entries. First is the string with the peakflow filename. Second is an integer with value 0 if the filename was set successfully. 1 otherwise.

This routine gets the filename of the output peakflow file. If a database connection is used, then the contents of *peakflowname* is not modified. The Python interface routine returns the string with the filename instead of taking it as an argument.

Asynch_Set_Snapshot_Output_Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Snapshot_Output_Name(
asynchsolver* asynch,
char* snapshotname)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* snapshotname: Name to set the snapshot output filename.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

::

Python Interface:
Set_Snapshot_Output_Name(snapshotname)

-  Inputs

-  snapshotname: Name to set the snapshot output filename.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

This routine sets the filename of the output snapshot file. If a database connection is used, then no changes are made.

Asynch_Get_Snapshot_Output_Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Get_Snapshot_Output_Name(
asynchsolver* asynch,
char* snapshotname)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* snapshotname: Name of the snapshot output filename returned.

-  Outputs

-  int: 0 if the filename was set successfully. 1 otherwise.

::

Python Interface:
Get_Snapshot_Output_Name()

-  Outputs

-  Returns a list with two entries. First is the string with the snapshot filename. Second is an integer with value 0 if the filename was set successfully. 1 otherwise.

This routine gets the filename of the output snapshot file. If a database connection is used, then the contents of *snapshotname* is not modified. The Python interface routine returns the string with the filename instead of taking it as an argument.

Asynch_Get_Size_Global_Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
unsigned int Asynch_Get_Size_Global_Parameters(
asynchsolver* asynch)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  Outputs

-  unsigned int: The number of global parameters, typically specified in a global file.

::

Python Interface:
Get_Size_Global_Parameters()

-  Outputs

-  Returns the number of global parameters, typically specified in a global file.

This routine returns the number of parameters for the differential equations that are available to every link. This number is either specified in a global file, or with a call to the routine *Asynch_Set_Global_Parameters* (see Section [sec: asynch\ :sub:`s`\ et\ :sub:`g`\ lobal\ :sub:`p`\ arameters]).

Asynch_Get_Global_Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Get_Global_Parameters(
asynchsolver* asynch,double* gparams)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double\* gparams: Array to store the global parameters.

-  Outputs

-  int: Returns 1 if an error occurred, 0 otherwise.

::

Python Interface:
Get_Global_Parameters()

-  Outputs

-  Returns a list with two entries: A list containing the global parameters, and an integer with value 1 if an error occurred and 0 otherwise.

This routine gets the values of the global parameters. For the C function, the array *gparams* should contain enough memory to hold the global parameters. The number of global parameters can be determined with a call to *Asynch_Get_Size_Global_Parameters* (see Section [sec: asynch\ :sub:`g`\ et\ :sub:`s`\ ize\ :sub:`g`\ lobal\ :sub:`p`\ arameters]). The obtained values are a copy of Asynch’s internal values of the parameters, so modifying them will have no effect on the solver results.

Asynch_Set_Global_Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Global_Parameters(
asynchsolver* asynch,
double* gparams,
unsigned int n)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  double\* gparams: Array of global parameter values to set.

-  unsigned int n: The number of global parameters in gparams.

-  Outputs

-  int: Returns 1 if an error occurred, 0 otherwise.

::

Python Interface:
Set_Global_Parameters(gparams)

-  Inputs

-  gparams: List of global parameter values to set.

-  Outputs

-  Returns 1 if an error occurred, 0 otherwise.

This routine sets the values of the global parameters. The number of global parameters may be different from the number previously stored by Asynch. However, if the new global parameters should be consistent with the model equations at each link (for example, if the model at each hillslope expects 4 global parameters, using this routine to change the number of global parameters to 2 WILL cause an error). The values are copied from the array (or list) into Asynch’s internal memory. So changing their values after this function call will NOT change the solver results.

Custom Outputs
--------------

Time series outputs can be customized so as to produce results specific to a particular model. For instance, fluxes that are used internally to the model, or values that are not needed at all for computing state solutions, may be outputted.

To create a custom output, specify the name of your new output in the global file (see Section [sec: time series output]). Next, in your program after reading the global file but before performing simulations, call the routine *Asynch_Set_Output* to set your output. Then perform any calculations or further modifications as usual.

A custom routine must be provided by the user for determining how the output is calculated. The specifications for this function are provided below.

Similarly, custom peakflow outputs can be created. These are set with a call to the routine *Asynch_Set_Peakflow_Output*.

The routines for setting custom outputs are described below.

Asynch_Set_Output
~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Output(
asynchsolver* asynch,
char* name,
short int data_type,
void (*func)(double,VEC*,VEC*,
VEC*,IVEC*,int,void*),
int* used_states,
int num_states)

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* name: Name of the custom time series.

-  short int data_type: The data type for the time series.

-  void (\*func): Pointer to the routine to call when writing data for the custom output.

-  int\* used_states: Array of all indices in the state vector used by *func*.

-  int num_states: The number of states from the state vector used by *func*.

-  Return value

-  int: 0 if the output is set, 1 if an error occurred.

::

Python Interface:
Set_Output(name,data_type,func,used_states_list)

-  Inputs

-  name: Name of the custom time series.

-  data_type: The data type for the time series.

-  func: Routine to call when writing data for the custom output.

-  used_states_list: List of the all indices in the state vector used by *func*.

-  Return value

-  0 if the output is set, 1 if an error occurred.

This routine sets a custom output time series. A function is required that will be called every time output data is to be written for a link. The function set in this routine should have the specification:

::

C Interface:
data_type func(double t,VEC* y,
VEC* global_params,VEC* params,
IVEC* iparams,int state,void* user)

-  Inputs

-  double t: The current time.

-  VEC\* y: The current state vector at time *t* for the link.

-  VEC\* global_params: The vector of parameters uniform amongst all links.

-  VEC\* params: The vector of parameters for the link.

-  IVEC\* iparams: The vector of integer parameters for the link.

-  int state: The current state of the state vector.

-  void\* user: User defined data.

-  Return value

-  data_type: Returns the data to be written as output.

::

Python Interface:
func(t,y,global_params,params,iparams,state,user)

-  Inputs

-  t: The current time.

-  y: The current state vector at time *t* for the link.

-  global_params: The vector of parameters uniform between all links.

-  params: The vector of parameters for the link.

-  iparams: The vector of integer parameters for the link.

-  state: The current state of the state vector.

-  user: User defined data.

-  Return value

-  Returns the data to be written as output.

*Asynch_Set_Output* must also be given an array *used_states*. This array contains the indices of the states in the state vectors which are needed by the user specified routine *func*. All states listed in this array are guaranteed to be available for the routine *func*. The Python version of this function should be of type *ASYNCH_OUTPUT_<output type>_DATATYPE*, with <output type> either *INT* or *DOUBLE*.

Asynch_Set_Peakflow_Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

C Interface:
int Asynch_Set_Peakflow_Output(
asynchsolver* asynch,
char* name,
void (*func)(unsigned int,double,VEC*,VEC*,
VEC*,double,unsigned int,void*,char*))

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  char\* name: Name of the custom peakflow function.

-  void (\*func): Pointer to the routine to call when writing peakflow data for the custom output.

-  Return value

-  int: 1 if the output is set, 0 if an error occurred.

::

Python Interface:
Set_Peakflow_Output(name,func)

-  Inputs

-  name: Name of the custom peakflow function.

-  func: Routine to call when writing peakflow data for the custom output.

-  Return value

-  1 if the output is set, 0 if an error occurred.

This routine sets a custom peakflow function. A function is required that will be called every time peakflow data is to be written for a link. The function set in this routine should have the specification:

::

C Interface:
void func(unsigned int ID,
double peak_time,
VEC* peak_value,
VEC* params,
VEC* global_params,
double conversion,
unsigned int area_idx,
void* user,
char* buffer)

-  Inputs

-  unsigned int ID: Link ID of the link to write data.

-  double peak_time: The local time at which the peak flow was found.

-  VEC\* peak_value: The vector of states at the time when the peakflow was determined.

-  VEC\* params: The vector of parameters for the link.

-  VEC\* global_params: The vector of parameters uniform between all links.

-  double conversion: A unit conversion factor to apply to the upstream area, if desired.

-  unsigned int area_idx: The index in *params* containing the upstream area of the link.

-  void\* user: User defined data.

-  char\* buffer: Buffer to write the output data.

::

Python Interface:
func(ID,peak_time,peak_value,params,global_params,
conversion,area_idx,user,outputbuffer)

-  Inputs

-  ID: Link ID of the link to write data.

-  peak_time: The local time at which the peak flow was found.

-  peak_value: The vector of states at the time when the peak flow was determined.

-  params: The vector of parameters for the link.

-  global_params: The vector of parameters uniform between all links.

-  conversion: A unit conversion factor to apply to the upstream area, if desired.

-  area_idx: The index in *params* containing the upstream area of the link.

-  user: User defined data.

-  outputbuffer: String buffer to write the output data.

This user-defined function fills *outputbuffer* with the desired peakflow information. Whatever contents are placed into the output buffer will appear as a row or tuple in the final peakflow data. As an example, the following function is used when the peakflow function “Classic” is selected:

::

void OutputPeakflow_Classic_Format(unsigned int ID,
double peak_time,VEC* peak_value,VEC* params,
VEC* global_params,double conversion,
unsigned int area_idx,void* user,char* buffer)
{
sprintf(buffer,"%u %.4f %.8f %.8f\n",ID,
conversion*params->ve[area_idx],peak_time,
peak_value->ve[0]);
}

Each line of the output peakflow file will contain: link id, link upstream area (from the parameter list), the time to peak, and the value of the peakflow.

Custom Models
-------------

ASYNCH allows users to add their own models by modifying the routines described in Section [sec: model descriptions]. However, the ASYNCH user interface allows users to create their own routines to define a model, and pass them to the solvers. This allows users to create their own models without modifying the underlying solver code. This can be done with a call to the following routine:

::

C Interface:
int Asynch_Custom_Model(asynchsolver* asynch,
void (*SetParamSizes)(UnivVars*,void*),
void (*Convert)(VEC*,unsigned int,void*),
void (*Routines)(Link*,unsigned int,
unsigned int,unsigned short int,void*),
void (*Precalculations)(Link*,VEC*,VEC*,IVEC*,
unsigned int,unsigned int,unsigned short int,
unsigned int,void*),
unsigned int (*InitializeEqs)(VEC*,VEC*,IVEC*,
QVSData*,unsigned short int,VEC*,
unsigned int,void*))

-  Inputs

-  asynchsolver\* asynch: Pointer to the *asynchsolver* object to use.

-  void (\*SetParamSizes): Routine to set value needed to describe the model.

-  void (\*Convert): Routine to apply unit conversions to link parameters.

-  void (\*Routines): Routine to specify model routines.

-  void (\*Precalculations): Routine for model precalculations.

-  unsigned int (\*InitializeEqs): Routine to set initial conditions not determined through the global file.

-  Return Value

-  int: Returns 1 if an error occurred, 0 if the model was set successfully.

::

Python Interface:
Custom_Model(SetParamSizes,Convert,Routines,
Precalculations,InitializeEqs)

-  Inputs

-  SetParamSizes: Routine to set value needed to describe the model.

-  Convert: Routine to apply unit conversions to link parameters.

-  Routines: Routine to specify model routines.

-  Precalculations: Routine for model precalculations.

-  InitializeEqs: Routine to set initial conditions not determined through the global file.

-  Return Value

-  Returns 1 if an error occurred, 0 if the model was set successfully.

The five routines needed to create a model are described in Section [sec: model definition]. The routines for the differential and algebraic equations are set in the *Routines* routine, and do not need to be passed through the *Asynch_Custom_Model* routine directly.

Structures and Objects
----------------------

Several structures are defined for the ASYNCH solvers, and are available for use in user routines. The exact structure definitions can be found in the source files:

C Interface
: structs.h and mathmethods.h

Python Interface
: asynch_py/asynch_interface.py

Below we provide a description of structures and their members which are useful for the user to manipulate. Other structures are defined, but not likely to be utilized by a typical user. Descriptions of these structures can be found in the source code. Although the C versions of the structures are given here, the structures are similar for other languages.

Vector and Matrix Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vectors and matrices have the following members:

-  VEC

-  double\* ve

-  unsigned int dim

-  IVEC

-  int\* ve

-  unsigned int dim

-  MAT

-  double\* array

-  double\*\* me

-  unsigned int m,n

For the two vector structures *VEC* and *IVEC*, the array of *dim* vector elements is stored in the member *ve*. Similarly, the matrix structure *MAT* has its *m* by *n* elements stored in the member *me*. The elements of a *MAT* object can also be accessed linearly through the member *array*.

Typical users of ASYNCH will only need to modify elements of existing vectors (if that). Developers may need to use these structures more heavily. Routines for creating, destroying, and performing operations on vectors and matrices can be found in the source code *mathmethods.c*.

UnivVars Structure
~~~~~~~~~~~~~~~~~~

A single instance of a universal variables structure is maintained by the ASYNCH solver. This object holds much of the information that describes the general structure of the network. A user who creates a custom model will need to set several members of the universal variables object. Below are the members a typical user may encounter. Many, many more members exist, and a developer is strongly encouraged to view these members in *structs.h*.

unsigned short int type
: The index for the model used

int iter_limit
: The maximum number of time steps stored per link

int max_transfer_steps
: Maximum number of steps to communicate at once between processes

double maxtime
: The final local time for the simulation

double t_0
: Initial local time to start integration

unsigned int dim
: The dimension of the state vectors at each link

unsigned int discont_size
: Size of discontinuity buffers at each link

unsigned int diff_start
: Starting index of differential variables in solution vectors

unsigned int no_ini_start
: Starting index of differential variables not read from disk

unsigned short int uses_dam
: 1 if this type can use dams, 0 else

VEC\* global_params
: List of global parameters

unsigned int params_size
: Size of params at each link without a dam

unsigned int iparams_size
: Size of iparams at each link

unsigned int dam_params_size
: Size of params at each link with a dam

unsigned int disk_params
: Number of parameters to read from disk

unsigned int area_idx
: Index of upstream area (A_i) in params

unsigned int areah_idx
: Index of hillslope area (A_h) in params

unsigned int num_dense
: Number of states where dense output is calculated

unsigned int\* dense_indices
: List of indices in solution where dense output is needed

unsigned short int convertarea_flag
: 1 if hillslope and upstream areas are converted from :math:`km^2` to :math:`m^2`, 0 if not

unsigned short int template_flag
: 1 if XML templates are to be used, 0 if not

unsigned int num_forcings
: The number forcings the model requires

Link Structure
~~~~~~~~~~~~~~

Every link in the network is represented by an instance of a link structure. Users developing their own models need to set a few members of each link. Many more members exist, and developers are recommended to view *structs.h*.

unsigned int ID
: The unique ID of the link

void \*f
: The function for evaluating an ODE

void \*alg
: The function for evaluating algebraic equations

int \*state_check
: The function to check the discontinuity state of the state variables (for discontinuities)

void \*Jacobian
: The function for evaluating the Jacobian matrix of an ODE

int \*RKSolver
: The routine for the Runge Kutta solver

void \*CheckConsistency
: The function that alters a state vector to be within the constraints of the ODEs

The requirements for what tasks each of these functions perform is covered in Section [sec: model equations definition]. Each function should take the following arguments:

-
-  void alg(VEC\*,VEC\*,VEC\*,QVSData\*,int,VEC\*)

-  int state_check(VEC\*,VEC\*,VEC\*,QVSData\*,unsigned int)

-  void Jacobian(double,VEC\*,VEC\*\*,unsigned short int,VEC\*, double\*,VEC\*,MAT\*)

-  int RKSolver(Link\*,UnivVars\*,int\*,short int,FILE\*, ConnData\*,Forcing\*\*,TempStorage\*)

-  void CheckConsistency(VEC\* y,VEC\* params,VEC\* global_params)

If these functions are set with the Python interface, and are written in Python, then each must be decorated as an appropriate data type. See Section [sec: api for python].
