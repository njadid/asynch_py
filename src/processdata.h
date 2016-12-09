#ifndef PROCESSDATA_H
#define PROCESSDATA_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include "structs.h"
#include "sort.h"
#include "comm.h"
#include "rkmethods.h"
#include "data_types.h"

#define DB_CONNS_AT_ONCE 10

extern int np;
extern int my_rank;

int Process_Data(Link* sys, GlobalVars* GlobalVars, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, unsigned int** id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile);
int UploadHydrosDB(Link* sys, GlobalVars* GlobalVars, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, unsigned int** id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile);
//void DataDump(Link* sys,unsigned int N,int* assignments,UnivVars* GlobalVars,unsigned int last_file);
int DataDump2(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);
int DataDumpH5(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);
int UploadDBDataDump(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);
void PrepareDatabaseTable(GlobalVars* GlobalVars, ConnData* conninfo);
int PreparePeakFlowFiles(GlobalVars* GlobalVars, unsigned int peaksave_size);
int DumpPeakFlowData(Link* sys, GlobalVars* GlobalVars, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, unsigned int** id_to_loc, ConnData* conninfo);
int UploadPeakFlowData(Link* sys, GlobalVars* GlobalVars, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, unsigned int** id_to_loc, ConnData* conninfo);
int RemoveTemporaryFiles(GlobalVars* GlobalVars, unsigned int my_save_size, char* additional_temp);
FILE* PrepareTempFiles(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, char* additional, unsigned int** id_to_loc);
int ResetTempFiles(double set_time, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* GlobalVars, unsigned int my_save_size, unsigned int** id_to_loc);
int SetTempFiles(double set_time, void* set_value, enum AsynchTypes data_type, unsigned int component_idx, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* GlobalVars, unsigned int my_save_size, unsigned int** id_to_loc);
void LoadRecoveryFile(char* filename, Link* sys, unsigned int N, unsigned int my_N, unsigned int* assignments, GlobalVars* GlobalVars);
int ConvertBinaryToString(double* data_storage, char* submission, unsigned int blocks, unsigned int dimp1, unsigned int id);
int overwrite_last_step(Link* link_i, GlobalVars *GlobalVars, FILE* outputfile);

#endif

