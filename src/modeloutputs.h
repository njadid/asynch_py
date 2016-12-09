#ifndef MODELOUTPUTS_H
#define MODELOUTPUTS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"
#include "data_types.h"

extern int np;
extern int my_rank;

void SetOutputFunctions(char* outputname, char* specifier, unsigned int* states_used, unsigned int* num_states_used, short int* output_size, enum AsynchTypes* output_type, OutputIntCallback **output_i, OutputDoubleCallback **output_d);
void SetPeakflowOutputFunctions(char* outputname,void (**peak_output)(unsigned int,double,VEC,VEC,VEC,double,unsigned int,void*,char*));
short int GetByteSize(short int type);
void GetSpecifier(char* specifier,short int type);
unsigned int CalcTotalOutputSize(GlobalVars* GlobalVars);
int OutputsSet(GlobalVars* GlobalVars);
int PeakflowOutputsSet(GlobalVars* GlobalVars);

//Output functions *****************************************************************************
int Output_LinkID(unsigned int id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
double Output_Time(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State0(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State1(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State2(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State3(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State4(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State5(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State6(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
double Output_State7(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
int Output_Time_Int(unsigned int id, double t,VEC y_i,VEC global_params,VEC params,int state,void* user);

//Peakflow output functions***********************************************************************************
void OutputPeakflow_Classic_Format(unsigned int ID,double peak_time,VEC peak_value,VEC params,VEC global_params,double conversion,unsigned int area_idx,void* user,char* buffer);
void OutputPeakflow_Forecast_Format(unsigned int ID,double peak_time,VEC peak_value,VEC params,VEC global_params,double conversion,unsigned int area_idx,void* user,char* buffer);

#endif

