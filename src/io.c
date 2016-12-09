#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <string.h>

#include "io.h"

//Creates an io object
void OutputFunc_Init(unsigned short hydros_loc_flag, unsigned short peaks_loc_flag, unsigned short dump_loc_flag, OutputFunc* output_func)
{
    //Temporary Calculations
    output_func->PrepareTempOutput = &PrepareTempFiles;

    //Prepare Final Time Series Output
    if (hydros_loc_flag == 3)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->PrepareOutput = &PrepareDatabaseTable;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else
        output_func->PrepareOutput = NULL;

    //Prepare Peakflow Output
    if (peaks_loc_flag == 1)
        output_func->PreparePeakflowOutput = &PreparePeakFlowFiles;
    else
        output_func->PreparePeakflowOutput = NULL;

    //Create Final Time Series Output
    if (hydros_loc_flag == 1 || hydros_loc_flag == 2 || hydros_loc_flag == 4)
        output_func->CreateOutput = &Process_Data;
    else if (hydros_loc_flag == 3)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreateOutput = &DumpTimeSerieDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }       
    else
        output_func->CreateOutput = NULL;

    //Create Peakflow Output
    if (peaks_loc_flag == 1)
        output_func->CreatePeakflowOutput = &DumpPeakFlowText;
    else if (peaks_loc_flag == 2)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreatePeakflowOutput = &DumpPeakFlowDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }        
    else
        output_func->CreatePeakflowOutput = NULL;

    //Set data dump routines
    if (dump_loc_flag == 1)
        output_func->CreateSnapShot = &DumpStateText;
    else if (dump_loc_flag == 2)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreateSnapShot = &DumpStateDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else if (dump_loc_flag == 3)
        output_func->CreateSnapShot = &DumpStateH5;
    else
        output_func->CreateSnapShot = NULL;
}


//Reads a .dbc file and creates a corresponding database connection.
//This does NOT connect to the database.
//Returns NULL if there was an error.
void ReadDBC(char* filename, ConnData* const conninfo)
{
    bool res = true;
    unsigned int i = 0, j = 0;
    char connstring[ASYNCH_MAX_CONNSTRING_LENGTH];
    char c;

    //if(my_rank == 0)
    //{
    FILE* input = fopen(filename, "r");

    if (!input)
    {
        printf("Error opening .dbc file %s.\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (CheckWinFormat(input))
    {
        printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", filename);
        fclose(input);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //Read connection information
    //Currently, this expects 4 things like:
    fgets(connstring, ASYNCH_MAX_CONNSTRING_LENGTH, input);
    ConnData_Init(conninfo, connstring);
    if (!conninfo)
    {
        fclose(input);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //Get number of queries
    if (fscanf(input, "%u", &(conninfo->num_queries)) == EOF)
    {
        printf("[%i]: Error: failed to parse file.\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (i = 0; i < conninfo->num_queries; i++)
        conninfo->queries[i] = (char*)malloc(ASYNCH_MAX_QUERY_LENGTH * sizeof(char));

    //Get queries. They are delineated by a ;
    for (j = 0; j < conninfo->num_queries; j++)
    {
        //Get rid of whitespace
        c = fgetc(input);
        while (c != EOF && (c == ' ' || c == '\n' || c == '\t'))	c = fgetc(input);
        if (c == EOF)
        {
            printf("[%i]: Warning: did not see %u queries in %s.\n", my_rank, conninfo->num_queries, filename);
            break;
        }

        //Read in query
        for (i = 0; i < ASYNCH_MAX_QUERY_LENGTH - 2 && c != ';' && c != EOF; i++)
        {
            conninfo->queries[j][i] = c;
            c = fgetc(input);
        }

        //Check for problems and put stuff on the end
        if (i == ASYNCH_MAX_QUERY_LENGTH)
            printf("[%i]: Warning: query %u is too long in %s.\n", my_rank, j, filename);
        else if (c == EOF)
        {
            printf("[%i]: Warning: did not see %u queries in %s.\n", my_rank, conninfo->num_queries, filename);
            break;
        }
        else
        {
            conninfo->queries[j][i] = ';';
            conninfo->queries[j][i + 1] = '\0';
        }
    }

    fclose(input);

    ////Get string length for other procs
    //j = (unsigned int) strlen(connstring) + 1;
//}

////Check if an error occurred
//finish:
//MPI_Bcast(&has_error,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

////Transfer info from .dbc file
//if(!errorcode)
//{
//	MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//	MPI_Bcast(connstring,j,MPI_CHAR,0,MPI_COMM_WORLD);
//	if(my_rank != 0)
//		ConnData_Init(conninfo, connstring);
//	MPI_Bcast(&(conninfo->num_queries),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//	if(my_rank == 0)
//	{
//		for(i=0;i<conninfo->num_queries;i++)
//		{
//			j = (unsigned int) strlen(conninfo->queries[i]) + 1;
//			MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//			MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
//		}
//	}
//	else
//	{
//		conninfo->queries = (char**) malloc(conninfo->num_queries*sizeof(char*));
//		for(i=0;i<conninfo->num_queries;i++)
//		{
//			conninfo->queries[i] = (char*) malloc(ASYNCH_MAX_QUERY_LENGTH*sizeof(char));
//			MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//			MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
//		}
//	}
//}


}

//!!!! This really sucks. Is there any way to improve it? !!!!
//Writes a value to an ASCII file
void WriteValue(FILE* outputfile, char* specifier, char* data_storage, short int data_type, char* delim)
{
    switch (data_type)
    {
    case ASYNCH_DOUBLE:
        fprintf(outputfile, specifier, *(double*)data_storage);
        break;
    case ASYNCH_INT:
        fprintf(outputfile, specifier, *(int*)data_storage);
        break;
    case ASYNCH_FLOAT:
        fprintf(outputfile, specifier, *(float*)data_storage);
        break;
    case ASYNCH_SHORT:
        fprintf(outputfile, specifier, *(short int*)data_storage);
        break;
    case ASYNCH_CHAR:
        fprintf(outputfile, specifier, *(char*)data_storage);
        break;
    default:
        printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n", my_rank, data_type);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    fprintf(outputfile, delim);
}

unsigned int WriteStep(FILE* outputfile, unsigned int id, double t, VEC y, GlobalVars* GlobalVars, VEC params, unsigned int state, void* user, long int* pos_offset)
{
    unsigned int i;

    double output_d;
    //float output_f;
    //short int output_s;
    int output_i;
    //char output_c;
    long int total_written = 0;

    //Set file to current position
    //fsetpos(outputfile,pos);
    if (pos_offset)	fseek(outputfile, *pos_offset, SEEK_SET);

    //Write the step
    for (i = 0; i < GlobalVars->num_print; i++)
    {
        switch (GlobalVars->output_types[i])	//!!!! Get rid of this. Try char[] and output_sizes. !!!!
        {
        case ASYNCH_DOUBLE:
            output_d = (GlobalVars->outputs_d[i])(id, t, y, GlobalVars->global_params, params, state, user);
            fwrite(&output_d, sizeof(double), 1, outputfile);
            break;
        case ASYNCH_INT:
            output_i = (GlobalVars->outputs_i[i])(id, t, y, GlobalVars->global_params, params, state, user);
            fwrite(&output_i, sizeof(int), 1, outputfile);
            break;
        default:
            printf("[%i]: Error: Invalid output %s (%hi).\n", my_rank, GlobalVars->output_specifiers[i], GlobalVars->output_types[i]);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        total_written += GlobalVars->output_sizes[i];
        //if(pos_offset)	*pos_offset += GlobalVars->output_sizes[i];
    }

    if (pos_offset)	*pos_offset += total_written;
    return total_written;
}

unsigned int CatBinaryToString(char* submission, char* specifier, char* data_storage, short int data_type, char* delim)
{
    unsigned int written = 0;

    switch (data_type)
    {
    case ASYNCH_DOUBLE:
        written = sprintf(submission, specifier, *(double*)data_storage);
        break;
    case ASYNCH_INT:
        written = sprintf(submission, specifier, *(int*)data_storage);
        break;
    default:
        printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n", my_rank, data_type);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    sprintf(&(submission[written++]), delim);

    return written;
}


