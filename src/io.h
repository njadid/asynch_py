#ifndef IO_H
#define IO_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdbool.h>

#include "structs.h"
#include "comm.h"
#include "modeloutputs.h"
#include "processdata.h"

extern int np;
extern int my_rank;

void OutputFunc_Init(unsigned short hydros_loc_flag, unsigned short peaks_loc_flag, unsigned short dump_loc_flag, OutputFunc* output_func);

void ReadDBC(char* filename, ConnData* const conninfo);
void WriteValue(FILE* outputfile,char* specifier,char* data_storage,short int data_type,char* delim);
unsigned int WriteStep(FILE* outputfile, unsigned int id, double t, VEC y, GlobalVars* GlobalVars, VEC params, unsigned int state, void* user, long int* pos_offset);
unsigned int CatBinaryToString(char* submission,char* specifier,char* data_storage,short int data_type,char* delim);

#endif
