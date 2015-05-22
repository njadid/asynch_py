#include "forcings.h"


Forcing* InitializeForcings()
{
	Forcing* forcing = (Forcing*) malloc(sizeof(Forcing));
	forcing->filename = NULL;
	forcing->GlobalForcing = NULL;
	return forcing;
}


void FreeForcing(Forcing** forcings)
{
	if(forcings && *forcings)
	{
		if((*forcings)->filename)	free((*forcings)->filename);
		Destroy_ForcingData(&((*forcings)->GlobalForcing));
		free(*forcings);
		*forcings = NULL;
	}
}

//GetPasses ************************************************************************************
//Forcings (0 = none, 1 = .str, 2 = binary, 3 = database, 4 = .ustr, 5 = forcasting, 6 = .gz binary, 7 = recurring)

//For flag = 0,1,4
//!!!! Added maxtime on 02/17/15. I think this is needed for when calls to Asynch_Set_Total_Simulation_Time occur. !!!!
unsigned int PassesOther(Forcing* forcing,double maxtime)
{
	if(forcing->maxtime < maxtime)	return 2;
	else				return 1;
	//return 1;
}

//For flag = 2,6
unsigned int PassesBinaryFiles(Forcing* forcing,double maxtime)
{
	unsigned int passes = (forcing->last_file - forcing->first_file + 1) / forcing->increment;
	if((forcing->last_file - forcing->first_file + 1)%forcing->increment != 0)	passes++;
	return passes;
}

//For flag = 3
unsigned int PassesDatabase(Forcing* forcing,double maxtime)
{
	return (unsigned int) ((forcing->last_file - forcing->first_file + 1)/(60.0*forcing->increment*forcing->file_time)) + 1;
}

//For flag = 7
unsigned int PassesRecurring(Forcing* forcing,double maxtime)
{
	return (forcing->last_file - forcing->first_file)/31536000 + 1;
}

//GetNextForcing ********************************************************************************
//Forcings (0 = none, 1 = .str, 2 = binary, 3 = database, 4 = .ustr, 5 = forcasting, 6 = .gz binary, 7 = recurring)

//For flag = 0,1,4
double NextForcingOther(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	return GlobalVars->maxtime;
}

//For flag = 2
double NextForcingBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration;
	double maxtime;
	if(iteration == passes-1)	maxtime = GlobalVars->maxtime;
	else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
	int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*forcing->increment,(double) (forcing->last_file + 1));

	Create_Rain_Data_Par(sys,N,my_N,GlobalVars,my_sys,assignments,forcing->filename,forcing->first_file+iteration*forcing->increment,maxfileindex,iteration*forcing->file_time*forcing->increment,forcing->file_time,forcing,id_to_loc,forcing->increment+1,forcing_idx);

	(forcing->iteration)++;
	return maxtime;
}

//For flag = 6
double NextForcingGZBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration;
	double maxtime;
	if(iteration == passes-1)	maxtime = GlobalVars->maxtime;
	else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
	int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*forcing->increment,(double) (forcing->last_file + 1));

	Create_Rain_Data_GZ(sys,N,my_N,GlobalVars,my_sys,assignments,forcing->filename,forcing->first_file+iteration*forcing->increment,maxfileindex,iteration*forcing->file_time*forcing->increment,forcing->file_time,forcing,id_to_loc,forcing->increment+1,forcing_idx);

	(forcing->iteration)++;
	return maxtime;
}

//For flag = 3
double NextForcingDatabase(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration,first_timestamp = 0;
	double maxtime;

//printf("!!!! In here: %i %i\n",(int)(sys[my_sys[0]]->last_t*60 + .001) + forcing->raindb_start_time,(int)(forcing->first_file+iteration*forcing->file_time*60.0*forcing->increment+0.01));
//printf("!!!! increment = %u iteration = %u, first_file = %u\n",forcing->increment,iteration,forcing->first_file);

	if( (int)(sys[my_sys[0]]->last_t*60 + .001) + forcing->raindb_start_time == (int)(forcing->first_file+iteration*forcing->file_time*60.0*forcing->increment+0.01) )
	{
		if(iteration == passes-1)	maxtime = GlobalVars->maxtime;	//!!!! Is this really needed? !!!!
		else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
		int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*60*forcing->file_time*forcing->increment,(double) forcing->last_file);

		first_timestamp = forcing->first_file + (unsigned int)(iteration*forcing->file_time*60.0*forcing->increment + 0.01);

		//If first_timestamp is off from the database increment, then read one previous timestamp
		if( (int)(first_timestamp - forcing->good_timestamp) % (int)(forcing->file_time*60+.01))
			first_timestamp -= (unsigned int)(forcing->file_time*60.0 + 0.01);

		Create_Rain_Database(sys,N,my_N,GlobalVars,my_sys,assignments,db_connections[ASYNCH_DB_LOC_FORCING_START + forcing_idx],first_timestamp,maxfileindex,forcing,id_to_loc,GlobalVars->maxtime,forcing_idx);
		(forcing->iteration)++;
	}
	else
		maxtime = min(GlobalVars->maxtime,(iteration)*forcing->file_time*forcing->increment);

	return maxtime;
}

//For flag = 7
double NextForcingRecurring(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	double maxtime;
	struct tm next_time,*first_time;
	time_t casted_first_file = (time_t) forcing->first_file;

	//Set current_epoch
	if(forcing->iteration)
	{
		//Get the year for the initial timestamp
		first_time = gmtime(&casted_first_file);
		next_time.tm_sec = 0;
		next_time.tm_min = 0;
		next_time.tm_hour = 0;
		next_time.tm_mday = 1;
		next_time.tm_mon = 0;
		next_time.tm_year = first_time->tm_year + forcing->iteration;
		next_time.tm_isdst = 0;
	}
	else
	{
		first_time = gmtime(&casted_first_file);
		copy_tm(first_time,&next_time);	//Yeah, this is lazy. But it's only done once...
	}

	maxtime = CreateForcing_Monthly(sys,my_N,my_sys,GlobalVars,forcing->GlobalForcing,forcing_idx,&next_time,forcing->first_file,forcing->last_file,sys[my_sys[0]]->last_t);
	(forcing->iteration)++;
	return maxtime;
}


