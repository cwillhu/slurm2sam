#!/bin/bash
#
# For documentation on input parameters, see sacct man page
#
#Valid time formats:
#HH:MM[:SS] [AM|PM]
#MMDD[YY] or MM/DD[/YY] or MM.DD[.YY]
#MM/DD[/YY]-HH:MM[:SS]
#YYYY-MM-DD[THH:MM[:SS]]

sacct $@ \
      --format="Partition,JobID,State,NNodes,AllocCPUS,JobName,Submit,Start,End,ExitCode,MaxRSS,MaxVMSize,ReqMem,Suspended,Timelimit,Elapsed,TotalCPU,CPUTime,User" \
      --parsable2 



