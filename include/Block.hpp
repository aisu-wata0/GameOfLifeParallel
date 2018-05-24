#pragma once

#include "MPI.hpp"

// Neighbours
class Neigh {
public:
	int north, south;
	int west, east;

	int northwest, southwest;
	int northeast, southeast;

	Neigh()
		: north(-1), south(-1)
		, west(-1), east(-1)
		, northwest(-1), southwest(-1)
		, northeast(-1), southeast(-1)
	{
	}
};

template<typename number_t>
class Block : public Neigh {
public:
	int workerID;
// Memory arrays
	std::unique_ptr<number_t[]> src;
	std::unique_ptr<number_t[]> dst;
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
		// send
	MPI_Request reqSendNorth, reqSendSouth;
	MPI_Request reqSendWest, reqSendEast;
	// corners
	MPI_Request reqSendNorthwest, reqSendSouthwest;
	MPI_Request reqSendNortheast, reqSendSoutheast;
		// recv
	MPI_Request reqRecvNorth, reqRecvSouth;
	MPI_Request reqRecvWest, reqRecvEast;
	// corners
	MPI_Request reqRecvNorthwest, reqRecvSouthwest;
	MPI_Request reqRecvNortheast, reqRecvSoutheast;

	Block (int height, int width, int rowD, int colD, int workerID)
		: Neigh()

		, rows(height / rowD), cols(width / colD)
	// Bounds of what should be processed, existance of neighbours expand bounds
		, rowMin(2), rowMax(rows-1)
		, colMin(2), colMax(cols-1)

		, rows_m(rows+2), cols_m(cols+2)

		, workerID(workerID)
	// calculate position of this block in the global matrix
		, rowStart((workerID / colD) * rows), colStart((workerID % colD) * cols)
	// allocate local matrices
		, src(std::make_unique<number_t[]>(rows_m*cols_m))
		, dst(std::make_unique<number_t[]>(rows_m*cols_m))
	{
		if (src == NULL || dst == NULL) {
			std::cerr << "Not enough memory" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			return EXIT_FAILURE;
		}
	// Neighbours
		// not at the north border
		if (rowStart != 0){
			north = workerID - colD;
			--rowMin;
		}
		// not at the south border
		if (rowStart + rows != height){
			south = workerID + colD;
			++rowMax;
		}
		// not at the west border
		if (colStart != 0){
			west = workerID - 1;
			--colMin;
			// corners
			if(north != -1)
				northwest = west - colD;
			if(south != -1)
				southwest = west + colD;
		}
		// not at the east border
		if (colStart + cols != width){
			east = workerID + 1;
			++colMax;
			// corners
			if(north != -1)
				northeast = east - colD;
			if(south != -1)
				southeast = east + colD;
		}
	}

	inline void tradeBorders(){
		if (north != -1) { // send my north row
			if(LOG && workerID == logWorkerID) printf("North trade start\n");
			MPI_Isend(&(*this)[1][1], 1, row_t,
				north, STD_TAG, MPI_COMM_WORLD, &reqSendNorth);
			// request his border
			MPI_Irecv(&(*this)[0][1], 1, row_t,
				north, STD_TAG, MPI_COMM_WORLD, &reqRecvNorth);
		}
		if (south != -1) { // send my south row
			if(LOG && workerID == logWorkerID) printf("S trade start\n");
			MPI_Isend(&(*this)[block.rows][1], 1, row_t,
				south, STD_TAG, MPI_COMM_WORLD, &reqSendSouth);
			// request border
			MPI_Irecv(&(*this)[block.rows+1][1], 1, row_t,
				south, STD_TAG, MPI_COMM_WORLD, &reqRecvSouth);
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
		// TODO: Corners
		if (northwest != -1) { // send my corner pixel
			// receive his corner pixel
		}
		if (southwest != -1) { // send my corner pixel
			// receive his corner pixel
		}
		// Corners
		if (northeast != -1) { // send my corner pixel
			// receive his corner pixel
		}
		if (southeast != -1) { // send my corner pixel
			// receive his corner pixel
		}
	}
	
	// Applies filter to a single element on position x,y
	inline void convoluteElem(int x, int y, int cols_m, number_t filter[CONVN][CONVN]) {
		number_t value = 0;
		for (int i = x-1, k=0;   i <= x+1;   ++i, ++k){
			for (int j = y-1, l=0;   j <= y+1;   ++j, ++l){
				value += (*this)[i][j] * filter[k][l];
			}
		}
		
		number_t neib = value;
		
		if(cell(i,j) == LIVE)
		{
			if(neib < 2 || neib > 3)
			{
				cellTmp(i,j) = DEAD;
			} else {
				cellTmp(i,j) = LIVE;
			}
		}
		else // if cell(i,j) == DEAD
		{
			if(neib == 3) {
				cellTmp(i,j) = LIVE;
			} else {
				cellTmp(i,j) = DEAD;
			}
		}
	}

	inline void convolute(int rowMin, int rowMax, int colMin, int colMax, number_t filter[CONVN][CONVN]) {
		for(int i = rowMin ; i <= rowMax; ++i)
			for(int j = colMin ; j <= colMax ; ++j)
				convoluteElem(i, j, filter);
	}

	inline void convoluteLocal(number_t filter[CONVN][CONVN]){
	//  convolute(*src,*dst, rowMin, rowMax, colMin, colMax, cols_m, **filter)
		convolute(2, rows-1,    2, cols-1,  filter);
	}

	inline void convoluteBorders(number_t filter[CONVN][CONVN]){
		MPI_Status status;
		if (north != -1) {
			if(LOG && workerID == logWorkerID) printf("north Wait recv\n");
			MPI_Wait(&reqRecvNorth, &status);
			convolute(1, 1, 2, cols-1,  filter);
		}
		if (south != -1) {
			if(LOG && workerID == logWorkerID) printf("south Wait recv\n");
			MPI_Wait(&reqRecvSouth, &status);
			convolute(rows, rows,   2, cols-1,  filter);
		}
		if (west != -1) {
			if(LOG && workerID == logWorkerID) printf("west Wait recv\n");
			MPI_Wait(&reqRecvWest, &status);
			convolute(2, rows-1,    1, 1,   filter);
		}
		if (east != -1) {
			if(LOG && workerID == logWorkerID) printf("east Wait recv\n");
			MPI_Wait(&reqRecvEast, &status);
			convolute(2, rows-1,    cols, cols, filter);
		}

		convoluteCorners(filter);
	}

	inline void convoluteCorners(number_t filter[CONVN][CONVN]){
		if (northwest != -1){
			// if(LOG && workerID == logWorkerID) printf("northwest Wait recv\n");
			// MPI_Wait(&reqRecvNorthwest, &status);
			convolute(1, 1, 1, 1, filter);
			if(LOG && workerID == logWorkerID) printf("nw\n");
		}
		if (southwest != -1){
			// if(LOG && workerID == logWorkerID) printf("southwest Wait recv\n");
			// MPI_Wait(&reqRecvSouthwest, &status);
			convolute(rows, rows, 1, 1, filter);
			if(LOG && workerID == logWorkerID) printf("sw\n");
		}
		if (northeast != -1){
			// if(LOG && workerID == logWorkerID) printf("northeast Wait recv\n");
			// MPI_Wait(&reqRecvNortheast, &status);
			convolute(1, 1, cols, cols, filter);
			if(LOG && workerID == logWorkerID) printf("ne\n");
		}
		if (southeast != -1){
			// if(LOG && workerID == logWorkerID) printf("southeast Wait recv\n");
			// MPI_Wait(&reqRecvSoutheast, &status);
			convolute(rows, rows, cols, cols, filter);
			if(LOG && workerID == logWorkerID) printf("se\n");
		}
	}

	inline void waitSendBorders(){
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
		// TODO: Wait corners
	}

	void print(){
		for(int i = rowMin; i <= rowMax; ++i){
			for(int j = colMin; j <= colMax; ++j){
				std::cout << block[i][j]) << " ";
			}
			std::cout << std::endl;
		}
	}

	void swapPointers(){
		src.swap(dst);
	}
	
	inline char& cell(long i, long j){
		return (*this)[i][j];
	}
	
	inline char& cellTmp(long i, long j){
		return at_dst(i,j);
	}

	number_t & at_dst(int i, int j){
		return dst[i*cols_m + j];
	}
	number_t * operator[](int i){
		return &src[i*cols_m];
	}
};
