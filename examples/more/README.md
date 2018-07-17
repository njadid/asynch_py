# ASYNCH - Examples - More

These are additional examples of input files for models considered secondary or still under testing phase.

## Usage

Each folder contains the specific files for each model profile.

*Asynch* is expected to be called having the model folder as the working directory.

Output files will be written in the respective *results* folder and can be compared with the content in *results\_benchmark* for reference.

## For developers

It is considered a good practice that every new model created presents its own set of example input files. It is strongly recommended that realistic global parameters are given within the input files.

Also, every time an existing model is modified it should be tested using currently existing input files by comparing the results from a new run after compilation with respective *results\_benchmark* folder.