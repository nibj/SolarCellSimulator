#!/bin/bash
echo Starting Computation
nohup matlab -nodesktop -nosplash < RunOptimization34reverse.m > logfile0.txt 2> err0.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile1.txt 2> err1.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile2.txt 2> err2.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile3.txt 2> err3.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile4.txt 2> err4.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile5.txt 2> err5.txt &
nohup matlab -nodesktop -nosplash < RunSlave.m > logfile6.txt 2> err6.txt &