#pragma once

#include <memory>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <iomanip>


#ifdef LIKWID_PERFMON
	#include <likwid.h>
#else
	#define LIKWID_MARKER_INIT
	#define LIKWID_MARKER_THREADINIT
	#define LIKWID_MARKER_SWITCH
	#define LIKWID_MARKER_REGISTER(regionTag)
	#define LIKWID_MARKER_START(regionTag)
	#define LIKWID_MARKER_STOP(regionTag)
	#define LIKWID_MARKER_CLOSE
	#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif


#include "mpi.h"

#define STD_TAG 0
#define CONVN 3

#define masterID 0
#define logWorkerID 0
#define LOG false

// MPI
int workerID, nWorkers;
MPI_Datatype row_t;
MPI_Datatype column_t;
MPI_Datatype block_t;
MPI_Datatype block2_t;
MPI_Datatype submatrix_t;
