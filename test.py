from py.asynch_interface import asynchsolver
import sys

#Prepare system
asynch = asynchsolver()
asynch.Parse_GBL(sys.argv[1])
asynch.Load_System()

#Prepare outputs
asynch.Prepare_Temp_Files()
asynch.Write_Current_Step()
asynch.Prepare_Peakflow_Output()
asynch.Prepare_Output()

#Advance solver
asynch.Advance(1)

#Take a snapshot
asynch.Take_System_Snapshot(None)

#Create output files
asynch.Create_Output(None)
asynch.Create_Peakflows_Output()

#Cleanup
asynch.Delete_Temporary_Files()