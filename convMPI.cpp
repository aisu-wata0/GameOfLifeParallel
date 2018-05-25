#include <memory>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include "mpi.h"
#define STD_TAG 0
#define CONVN 3

#define masterID 0
#define logWorkerID 0
#define LOG 0

// MPI
int workerID, nWorkers;
MPI_Datatype row_t;
MPI_Datatype column_t;
MPI_Datatype block_t;
MPI_Datatype block2_t;
MPI_Datatype submatrix_t;

// Neighbours
class Neigh {
public:
	int north, south;
	int west, east;

	Neigh()
		: north(-1), south(-1), west(-1), east(-1)
	{
	}
};

template<typename number>
class Block : public Neigh {
public:
	int workerID;
// Memory arrays
	std::unique_ptr<number[]> src;
	std::unique_ptr<number[]> dst;
// Dimensions
	int rows, cols;
// Bounds of what should be processed, existance of neighbours change bounds
	int rowMin, rowMax;
	int colMin, colMax;
// Memory dimensions
	int rows_m, cols_m;
// Index of the first element of the block in the global matrix
	int rowStart, colStart;

	// MPI requests
	MPI_Request reqSendNorth;
	MPI_Request reqSendSouth;
	MPI_Request reqSendWest;
	MPI_Request reqSendEast;
	/* add Corners if necessary */
	MPI_Request reqRecvNorth;
	MPI_Request reqRecvSouth;
	MPI_Request reqRecvWest;
	MPI_Request reqRecvEast;

	Block (int height, int width, int rowD, int colD, int workerID)
		: Neigh()

		, rows(height / rowD), cols(width / colD)
	// Bounds of what should be processed, existance of neighbours change bounds
		, rowMin(2), rowMax(rows-1)
		, colMin(2), colMax(cols-1)

		, rows_m(rows+2), cols_m(cols+2)

		, workerID(workerID)

		, rowStart((workerID / colD) * rows), colStart((workerID % colD) * cols)

		, src(std::make_unique<number[]>(rows_m*cols_m))
		, dst(std::make_unique<number[]>(rows_m*cols_m))
	{
		if (src == NULL || dst == NULL) {
			std::cerr << "Not enough memory" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			return EXIT_FAILURE;
		}
	// Neighbours
		// not at the top
		if (rowStart != 0){
			north = workerID - colD;
			--rowMin;
		}
		// not at the bottom
		if (rowStart + rows != height){
			south = workerID + colD;
			++rowMax;
		}
		// not at the right
		if (colStart != 0){
			west = workerID - 1;
			++colMax;
		}
		// not at the left
		if (colStart + cols != width){
			east = workerID + 1;
			++colMin;
		}
	}

	number maxDelta(){
		number maxDelta
			= abs(src[(rowMin)*cols_m + colMin] - dst[(rowMin)*cols_m + colMin]);

		for(int i = rowMin ; i <= rowMax; ++i){
			for(int j = colMin ; j <= colMax ; ++j){
				number delta = fabs(src[(i)*cols_m + j] - dst[(i)*cols_m + j]);
				if(delta > maxDelta){
					maxDelta = delta;
				}
			}
		}
		return maxDelta;
	}

	void tradeBorders(){
		if (north != -1) { // send my north row
			if(LOG && workerID == logWorkerID) printf("North trade start\n");
			MPI_Isend(&(*this)[1][1], 1, row_t,
				north, STD_TAG, MPI_COMM_WORLD, &reqSendNorth);
			// request his border
			MPI_Irecv(&(*this)[0][1], 1, row_t,
				north, STD_TAG, MPI_COMM_WORLD, &reqRecvNorth);
			// Corners
			if (west != -1) { // send my corner pixel
				// receive his corner pixel
			}
			if (east != -1) { // send my corner pixel
				// receive his corner pixel
			}
		}
		if (south != -1) { // send my south row
			if(LOG && workerID == logWorkerID) printf("S trade start\n");
			MPI_Isend(&(*this)[block.rows][1], 1, row_t,
				south, STD_TAG, MPI_COMM_WORLD, &reqSendSouth);
			// request border
			MPI_Irecv(&(*this)[block.rows+1][1], 1, row_t,
				south, STD_TAG, MPI_COMM_WORLD, &reqRecvSouth);
			// Corners
			if (west != -1) { // send my corner pixel
				// receive his corner pixel
			}
			if (east != -1) { // send my corner pixel
				// receive his corner pixel
			}
		}
		if (west != -1) { // send my west col
			if(LOG && workerID == logWorkerID) printf("W trade start\n");
			MPI_Isend(&(*this)[1][1], 1, column_t,
				 west, STD_TAG, MPI_COMM_WORLD, &reqSendWest);
			// request his border
			MPI_Irecv(&(*this)[1][0], 1, column_t,
				 west, STD_TAG, MPI_COMM_WORLD, &reqRecvWest);
		}
		if (east != -1) { // send my east col
			if(LOG && workerID == logWorkerID) printf("E trade start\n");
			MPI_Isend(&(*this)[1][block.cols], 1, column_t,
				  east, STD_TAG, MPI_COMM_WORLD, &reqSendEast);
			// request his border
			MPI_Irecv(&(*this)[1][block.cols+1], 1, column_t,
				east, STD_TAG, MPI_COMM_WORLD, &reqRecvEast);
		}
	}
	// Applies filter to a single element on position x,y
	void convoluteElem(int x, int y, int cols_m, number filter[CONVN][CONVN]) {
		number value = 0;
		for (int i = x-1, k=0;   i <= x+1;   ++i, ++k){
			for (int j = y-1, l=0;   j <= y+1;   ++j, ++l){
				value += (*this)[i][j] * filter[k][l];
			}
		}
		mod(x, y) = value;
	}

	void convolute(int rowMin, int rowMax, int colMin, int colMax, number filter[CONVN][CONVN]) {
		for(int i = rowMin ; i <= rowMax; ++i)
			for(int j = colMin ; j <= colMax ; ++j)
				convoluteElem(i, j, filter);
	}

	void convoluteLocal(number filter[CONVN][CONVN]){
	//  convolute(*src,*dst, rowMin, rowMax, colMin, colMax, cols_m, **filter)
		convolute(2, rows-1,    2, cols-1,  filter);
	}

	void convoluteBorders(number filter[CONVN][CONVN]){
		MPI_Status status;
		if (north != -1) {
			if(LOG && workerID == logWorkerID) printf("N Wait recv\n");
			MPI_Wait(&reqRecvNorth, &status);
			convolute(1, 1, 2, cols-1,  filter);
		}
		if (south != -1) {
			if(LOG && workerID == logWorkerID) printf("S Wait recv\n");
			MPI_Wait(&reqRecvSouth, &status);
			convolute(rows, rows,   2, cols-1,  filter);
		}
		if (west != -1) {
			if(LOG && workerID == logWorkerID) printf("W Wait recv\n");
			MPI_Wait(&reqRecvWest, &status);
			convolute(2, rows-1,    1, 1,   filter);
		}
		if (east != -1) {
			if(LOG && workerID == logWorkerID) printf("E Wait recv\n");
			MPI_Wait(&reqRecvEast, &status);
			convolute(2, rows-1,    cols, cols, filter);
		}

		convoluteCorners(filter);
	}

	void convoluteCorners(number filter[CONVN][CONVN]){
		if (north != -1 && east != -1){
			convolute(1, 1, cols, cols, filter);
			if(LOG && workerID == logWorkerID) printf("ne\n");
		}
		if (south != -1 && east != -1){
			convolute(rows, rows, cols, cols, filter);
			if(LOG && workerID == logWorkerID) printf("se\n");
		}
		if (south != -1 && west != -1){
			convolute(rows, rows, 1, 1, filter);
			if(LOG && workerID == logWorkerID) printf("sw\n");
		}
		if (north != -1 && west != -1){
			convolute(1, 1, 1, 1, filter);
			if(LOG && workerID == logWorkerID) printf("nw\n");
		}
	}

	void waitSendBorders(){
		MPI_Status status;
		if (north != -1){
			if(LOG && workerID == logWorkerID) printf("N Wait send\n");
			MPI_Wait(&reqSendNorth, &status);
		}
		if (south != -1){
			if(LOG && workerID == logWorkerID) printf("S Wait send\\n");
			MPI_Wait(&reqSendSouth, &status);
		}
		if (west != -1){
			if(LOG && workerID == logWorkerID) printf("W Wait send\\n");
			MPI_Wait(&reqSendWest, &status);
		}
		if (east != -1){
			if(LOG && workerID == logWorkerID) printf("E Wait send\\n");
			MPI_Wait(&reqSendEast, &status);
		}
	}

	void print(){
		for(int i = 1; i <= block.rows; ++i){
			for(int j = 1; j <= block.cols; ++j){
				printf("%.0lf  ", block[i][j]);
			}
			printf("\n");
		}
	}

	void swapPointers(){
		src.swap(dst);
	}

	number & mod(int i, int j){
		return dst[i*cols_m + j];
	}
	number * operator[](int i){
		return &src[i*cols_m];
	}
};

// Copies submatrix src to matrix dst. src start at dst at index[rowStart][colStart]
// src will end at index[rowStart+rowsB-1][colStart+colsB-1]
// Assumes src has border of 1 col and row (on top/down/left/right)
void copyMatrix(double *dst, int dstCols_m, int rowStart, int colStart, double *src, int rowsB, int colsB){
	for(int i = 1 ; i <= rowsB; ++i){
		for(int j = 1 ; j <= colsB ; ++j){
			dst[(i+rowStart-1)*dstCols_m  +  j+colStart-1] = src[(i)*(colsB+2) + j];
		}
	}
}
void copyMatrix_Test(){
	double *dst;
	int dstCols_m = 6;
	dst = calloc(dstCols_m*dstCols_m, sizeof(double));

	double *src;
	int rowsB = 3;
	int colsB = 3;
	int colsB_m = colsB +2;
	src = calloc(rowsB*colsB, sizeof(double));
	for(int i = 1; i <= rowsB; ++i)
		for(int j = 1; j <= colsB; ++j)
			src[i*colsB_m +j] = i + j;

	copyMatrix(dst, dstCols_m, 2, 2, src, rowsB, colsB);
	for(int i = 0 ; i < dstCols_m; ++i){
		for(int j = 0 ; j < dstCols_m ; ++j){
			printf("%.0lf\t", dst[i*dstCols_m + j]);
		}
		printf("\n");
	}
}

void parseArgs(int argc, char **argv,/* char **filename, int *width, int *height,*/ int *loops) {
	int ai = 1, minArgs = 1;
	if (argc >= minArgs) {
/*      *filename = malloc((strlen(argv[ai])+1) * sizeof(char));
		strcpy(*filename, argv[ai]);
		++ai;

		*width = atoi(argv[ai]);
		++ai;

		*height = atoi(argv[ai]);
		++ai;
 */
		if (argc >= minArgs+1){
			*loops = atoi(argv[ai]);
			++ai;
		}
	} else {
		fprintf(stderr, "\nUsage: %s [FILE] [width] [height] [loops :optional]\n", argv[0]);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		exit(EXIT_FAILURE);
	}
}

/* Decompose matrix dimensions as to minimize perimeter
 * Returns the height of the decomposed grid (row decomposition)
 * workers/rowD == colD
 */
int divideBlocks(int rows, int cols, int nWorkers) {
	int perim, rowD, bestrowD = 0;
	int perimMin = rows + cols + 1;
	for (rowD = 1; rowD <= nWorkers; ++rowD) {
		// check decomposition remainders
		if ((nWorkers % rowD != 0) || (rows % rowD != 0))
			continue;
// to properly divide work, columns have to be divided in (colD) blocks
// colD is such that rowD*colD = nWorkers = number of blocks in total
		int colD = nWorkers / rowD;
		if (cols % colD != 0)
			continue;
		// perimiter, each block will have (rows/rowD) rows
		// perimiter to the upper, left, down, right blocks
		// doesn't count twice (up+down), as its used only to compare
		perim = rows/rowD + cols/colD;
		if (perim <= perimMin) { // <= gives preference to higher rowD
				perimMin = perim;
				bestrowD = rowD;
		}
	}

	return bestrowD;
}

int main(int argc, char** argv) {
	int width = 8, height = width, loops = INT_MAX;
	//char *filename;

	// MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nWorkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &workerID);

	MPI_Status status;

	int rowD, colD;
	if (workerID == masterID) {
		//parseArgs(argc, argv, &filename, &width, &height, &loops);
// Each worker has a block
// The matrix is decomposed into a matrix of rowD * colD blocks
// rowD * colD == # of workers == # of blocks
// each block has rowsB * colsB elements
// 0 1 2    // where the numbers are enumerated blocks (same as workerID)
// 3 4 5    // colD = 3, 3 columns of blocks
// 6 7 8    // if workerID = 4, the north worker is 4 - colD = 1
		rowD = divideBlocks(height, width, nWorkers);

		if(rowD == 0){
			fprintf(stderr, "%s: Could not divide data to workers\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			return EXIT_FAILURE;
		}

		colD = nWorkers / rowD;

		if(LOG)
		{
			printf("loops = %d; rowD = %d; colD = %d\n", loops, rowD, colD);
			printf("nWorkers = %d;\n", nWorkers);
		}
	}
/**
	if (workerID != masterID) {
		filename = malloc((strlen(argv[1])+1) * sizeof(char));
		strcpy(filename, argv[1]);
	}
/**/
// Broadcast parameters
	// MPI_Bcast(&filename, strlen(filename)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&loops, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rowD, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&colD, 1, MPI_INT, 0, MPI_COMM_WORLD);
/**/
	Block<double> block(height, width, rowD, colD, workerID);

// Memory size of block
// block's first & last lines are the borders handled by other workers
// same for columns
// the other worker's elements are in the block row/col [0] (and row/col [row/colB+2])
// My first row is [1][*]
// row data type
	MPI_Type_contiguous(block.cols, MPI_DOUBLE, &row_t);
	MPI_Type_commit(&row_t);
// column data type
	MPI_Type_vector(block.rows, 1, block.cols_m, MPI_DOUBLE, &column_t);
	MPI_Type_commit(&column_t);
// Block data type
	MPI_Type_vector(block.rows, block.cols, block.cols_m, MPI_DOUBLE, &block_t);
	MPI_Type_commit(&block_t);
// Submatrix data type
	// Dimensions
	int globalDim[2] = {height,width};
	int subDim[2] = {block.rows,block.cols};
	int startIn[2] = {0,0}; // Start Index
	MPI_Type_create_subarray(2, globalDim, subDim, startIn,
		MPI_ORDER_C, MPI_DOUBLE, &block2_t);
	MPI_Type_create_resized(block2_t, 0, block.cols*sizeof(double), &submatrix_t);
	MPI_Type_commit(&submatrix_t);
/**/
// Scatter/Gather Setup
	// How many blocks to each worker
	int counts[rowD*colD];
	// displacement of each block
	int displs[rowD*colD];

	if(workerID == masterID){
		if (workerID == 0) {
			for (int i=0; i<rowD; i++) {
				for (int j=0; j<colD; j++) {
					displs[i*colD+j] = i*block.rows*colD  +  j;
					counts[i*colD+j] = 1;
				}
			}
		}
	}

/**/
// Allocate global matrix
	double *globalMatrix;
	if(workerID == masterID){
		globalMatrix = (double*)malloc(height*width * sizeof(double));
	}
// Initialize Block
	for(int i = 1; i <= block.rows; ++i){
		if(block.colStart == 0){
			block[i][1] = 10;
		} else {
			block[i][1] = 5;
		}
		block.mod(i, 1) = block[i][1];
		for(int j = 2; j <= block.cols-1; ++j){
			block[i][j] = 5;
			block.mod(i, j) = block[i][j];
		}
		if(block.colStart + block.cols == width){
			block[i][block.cols] = 0;
		} else {
			block[i][block.cols] = 5;
		}
		block.mod(i, block.cols) = block[i][block.cols];
	}
/**
	MPI_Scatterv(matrix, counts, displs, submatrix_t,
		&block[1][1], 1, block_t,
		0, MPI_COMM_WORLD);
/**/
	if(LOG && workerID == logWorkerID){
		printf("max loops = %d;\n", loops);
		printf("block.rows = %d; block.cols = %d;\n", block.rows, block.cols);
		printf("rowStart = %d; colStart = %d;\n", block.rowStart, block.colStart);
		printf("north = %d; south = %d;\n", block.north, block.south);
		printf("west = %d; east = %d;\n", block.west, block.east);
	}

// Filter
	double filter[CONVN][CONVN] = {{0,1,0},{1,1,1},{0,1,0}};
	double totalWeight = 0;
	for(int i = 0; i < CONVN; ++i)
		for(int j = 0; j < CONVN; ++j)
			totalWeight += filter[i][j];

	for(int i = 0; i < CONVN; ++i)
		for(int j = 0; j < CONVN; ++j)
			filter[i][j] /= totalWeight;

	MPI_Barrier(MPI_COMM_WORLD);

	double timeStart = MPI_Wtime();

	int ended = 0; // if loop should end
// Convolute Loop
	int t;
	for (t = 0; t < loops && !ended; ++t) {
	// Border trading
		if(LOG) printf("worker %d: start loop: %d\n", workerID, t);
	/**/
		block.tradeBorders();
	/**/
	// Local data
		if(LOG && workerID == logWorkerID) printf("Convolute\n");

		block.convoluteLocal(filter);
	/**/
	// Wait remotes and compute borders
		block.convoluteBorders(filter);
	/**/
		double maxDelta = block.maxDelta();
		if(LOG && workerID == logWorkerID) printf("max delta = %lf\n", maxDelta);
	/**/
	// Sync maxDelta to see if continues to loop
		if(LOG && workerID == logWorkerID) printf("Reducing maxDelta\n");

		double maxDeltaGlobal;
		MPI_Allreduce(&maxDelta, &maxDeltaGlobal, 1,
			MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		if(maxDeltaGlobal < 0.001){
			ended = 1;
			if(LOG && workerID == logWorkerID) printf("+ All workers finished\n");
		}
	/**
	// Old way to finish sync state
		int finishedWorkers = 0;
		if(maxDelta < 0.001){
			finishedWorkers = 1;
		}
		MPI_Request requests[2];
	#define FINISH_TAG 4
		if (workerID != masterID) {
			MPI_Isend(&finishedWorkers, 1, MPI_INT,  0,
				FINISH_TAG, MPI_COMM_WORLD, &requests[1]);
		} else {
		// Master
			for (int i = 1; i < nWorkers; ++i) {
				int finished = 0;
				MPI_Recv(&finished, 1, MPI_INT, MPI_ANY_SOURCE,
					FINISH_TAG, MPI_COMM_WORLD, &status);

				finishedWorkers += finished;
				if(LOG) printf("Worker %d finished = %d\n", status.MPI_SOURCE, finishedWorkers);
			}

			if(finishedWorkers == nWorkers){
				ended = 1;
				if(LOG) printf("+ All workers finished\n");
			}
		}

		// Broadcast end state from master
		MPI_Bcast(&ended, 1, MPI_INT, 0, MPI_COMM_WORLD);
	/**
	// Wait borders sent
	// Maybe useless because of the sync above
	// test performance MPI_Barrier  vs MPI_Wait
		MPI_Barrier(MPI_COMM_WORLD);
	/**/
		block.waitSendBorders();
		// Wait corners sent
	/**/
	// Swap pointers
		block.swapPointers();

		if(LOG) printf("worker %d: end loop: %d;\n", workerID, t);
	}
// Loop End


/**/
	if(workerID == masterID) printf("%d iterations\n", t);
	// if(LOG && workerID == logWorkerID) {
	//  block.print();
	// }
/**
// Print local finished matrices
	if(LOG)
	for (int p=0; p<nWorkers; p++) {
		if (workerID == p) {
			printf("Local process on workerID %d is:\n", workerID);
			for (int i=1; i<=size/sizeD; i++) {
				for (int j=1; j<=size/sizeD; j++) {
					printf("%.0lf ", block[i][j]);
				}
				printf("\n");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
/**/
// Gather matrix
	MPI_Gatherv(&block[1][1], 1,  block_t,
				 globalMatrix, counts, displs, submatrix_t,
				 0, MPI_COMM_WORLD);
/**
// Old way of gathering
	if (workerID != masterID) {
		MPI_Send(&block[1][1], 1, block_t, 0, 2, MPI_COMM_WORLD);
	} else {
	// Master
		globalMatrix = malloc(width*height * sizeof(double));
		copyMatrix(globalMatrix, width, rowStart, colStart, src, block.rows, block.cols);
		for (int i = 1 ; i != nWorkers ; ++i) {
			MPI_Recv(&block[1][1], 1, block_t, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
			int workerID = status.MPI_SOURCE;
			int rowStart = (workerID / colD) * block.rows;
			int colStart = (workerID % colD) * block.cols;
if(LOG)
{
	printf("received from worker %d\n", workerID);
	printf("rowStart = %d; colStart = %d;\n", rowStart, colStart);
	block.print();
}
			copyMatrix(globalMatrix, width, rowStart, colStart, src, block.rows, block.cols);
		}
	}
/**/
// Print Matrix
//  if(LOG)
	if(workerID == masterID){
		for(int i = 0; i < height; ++i){
			for(int j = 0; j < width; ++j){
				printf("%.0lf ", globalMatrix[i*width + j]);
			}
			printf("\n");
		}

		//free(globalMatrix);
	}
/**/
// Poll max timeElapsed
	double timeElapsed = MPI_Wtime() - timeStart;

	if(LOG && workerID == logWorkerID) printf("Reducing timeElapsed\n");

	double timeElapsedGlobal;
	MPI_reduce(&timeElapsed, &timeElapsedGlobal, 1,
		MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (workerID != masterID) {
	//  if(LOG)
		printf("timeElapsed %f\n", timeElapsedGlobal);
	}
/**
// Old way
	double timeElapsed = MPI_Wtime() - timeStart;

	if(LOG && workerID == logWorkerID) printf("Reducing timeElapsed\n");

#define TIME_TAG 3
	if (workerID != masterID) {
		MPI_Send(&timeElapsed, 1, MPI_DOUBLE, 0,
			TIME_TAG, MPI_COMM_WORLD);
	} else {
		double timeRemote;
		for (int i = 1 ; i != nWorkers ; ++i) {
			MPI_Recv(&timeRemote, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
				TIME_TAG, MPI_COMM_WORLD, &status);
			if (timeRemote > timeElapsed)
				timeElapsed = timeRemote;
		}
//      if(LOG)
		printf("timeElapsed %f\n", timeElapsed);
	}
/**/
	// Free memory
	MPI_Type_free(&column_t);
	MPI_Type_free(&row_t);
	MPI_Type_free(&block_t);
	MPI_Type_free(&submatrix_t);
/**/
	MPI_Finalize();
	return EXIT_SUCCESS;
}
