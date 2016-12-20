#ifndef SOLVERS_H
#define SOLVERS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"
#include "rainfall.h"
#include "comm.h"
#include "rkmethods.h"
#include "system.h"
#include "forcings.h"

extern int np;
extern int my_rank;

void Advance(
    Link* sys,
    unsigned int N,
    unsigned int* my_sys,
    unsigned int my_N,
    GlobalVars* globals,
    int* assignments,
    short int* getting,
    unsigned int* res_list,
    unsigned int res_size,
    unsigned int** id_to_loc,
    TempStorage* workspace,
    Forcing* forcings,
    ConnData* db_connections,
    TransData* my_data,
    bool print_flag,
    FILE* outputfile);

#endif
