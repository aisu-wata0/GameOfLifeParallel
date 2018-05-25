
#include "main.hpp"

using namespace G;

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

int main(int argc, char **argv)
{
	int loops = 1;
	int size = atoi(argv[1]);
	int width = size, height = size;
	float chanceAlive = atof(argv[2]);

	LIKWID_MARKER_INIT;

	// MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nWorkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &workerID);

	srand(0);

	int rowD, colD;
	if (workerID == masterID) {
		if(LOG)
		{
		std::cout << "loops=" << loops << std::endl;
		std::cout << "nWorkers=" << nWorkers << std::endl;
		std::cout << "height=" << height << std::endl;
		std::cout << "width=" << width << std::endl;
		}
		//parseArgs(argc, argv, &filename, &width, &height, &loops);
// Each worker has a block
// The matrix is decomposed into a matrix of rowD * colD blocks
// rowD * colD == # of workers == # of blocks
// each block has rowsB * colsB elements
// 0 1 2	// where the numbers are enumerated blocks (same as workerID)
// 3 4 5	// colD = 3, 3 columns of blocks
// 6 7 8	// if workerID = 4, the north worker is 4 - colD = 1
		rowD = divideBlocks(height, width, nWorkers);

		if(rowD == 0){
			fprintf(stderr, "%s: Could not divide data to workers\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			return EXIT_FAILURE;
		}

		colD = nWorkers / rowD;

		if(LOG)
		{
			std::cout << "rowD=" << rowD << std::endl;
			std::cout << "colD=" << colD << std::endl;
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
	Block<char> block(height, width, rowD, colD, workerID);

// Memory size of block
// block's first & last lines are the borders handled by other workers
// same for columns
// the other worker's elements are in the block row/col [0] (and row/col [row/colB+2])
// My first row is [1][*]
// row data type
	MPI_Type_contiguous(block.cols, MPI_CHAR, &row_t);
	MPI_Type_commit(&row_t);
// column data type
	MPI_Type_vector(block.rows, 1, block.cols_m, MPI_CHAR, &column_t);
	MPI_Type_commit(&column_t);
// Block data type
	MPI_Type_vector(block.rows, block.cols, block.cols_m, MPI_CHAR, &block_t);
	MPI_Type_commit(&block_t);
// Submatrix data type
	// Dimensions
	int globalDim[2] = {size,size};
	int subDim[2] = {block.rows,block.cols};
	int startIn[2] = {0,0}; // Start Index
	MPI_Type_create_subarray(2, globalDim, subDim, startIn,
		MPI_ORDER_C, MPI_CHAR, &block2_t);
	MPI_Type_create_resized(block2_t, 0, block.cols*sizeof(char), &submatrix_t);
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
					counts[i*colD+j] = 1;
					displs[i*colD+j] = i*block.rows*colD  +  j;
				}
			}
		}
	}

/**
	Game *g;

	if (workerID == masterID) {
		g = new Game(size);
		g->init(chanceAlive);
		g->print();
	} else {
		g = new Game(1);
	}
	MPI_Scatterv(&g->cell(0,0), counts, displs, submatrix_t,
		&block[1][1], 1, block_t,
		masterID, MPI_COMM_WORLD);
/**/
// Allocate global matrix
	char *globalMatrix;
	if(workerID == masterID){
		globalMatrix = (char*)malloc(height*width * sizeof(char));
	}
// Initialize it
	//#pragma omp for schedule(static)
	for(long i = 0; i < height; ++i)
	{
		for(long j = 0; j < width; ++j)
		{
			// globalMatrix[i*width + j] = (char)((i*width +j) % 99);
			double val = (double)rand() / RAND_MAX;
			if (val < chanceAlive)
				globalMatrix[i*width + j] = LIVE;
			else
				globalMatrix[i*width + j] = DEAD;
		}
	}

	for(long i = 0; i < height; ++i)
	{
		std::cout << "|";
		for(long j = 0; j < width; ++j)
		{
			// std::cout << std::setw(2) << (int)(*this)[i][j] << "  ";
			if(globalMatrix[i*width + j] == LIVE)
				std::cout << "O";
			else
				std::cout << ".";
		}
		std::cout << "|\n";
	}
	std::cout << std::flush;
/**/
	MPI_Scatterv(&globalMatrix[0], counts, displs, submatrix_t,
		&block[1][1], 1, block_t,
		masterID, MPI_COMM_WORLD);
/**/


	if(LOG && workerID == logWorkerID){
		printf("max loops = %d;\n", loops);
		printf("block.rows = %d; block.cols = %d;\n", block.rows, block.cols);
		printf("rowStart = %d; colStart = %d;\n", block.rowStart, block.colStart);
		printf("north = %d; south = %d;\n", block.north, block.south);
		printf("west = %d; east = %d;\n", block.west, block.east);
	}

	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (int p=0; p<nWorkers; p++) {
			if (workerID == p) {
				printf("Local process on workerID %d is:\n", workerID);
				block.print();
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
// Filter
	char filter[CONVN][CONVN] = {{1,1,1},{1,0,1},{1,1,1}};
//	double totalWeight = 0;
//	for(int i = 0; i < CONVN; ++i)
//		for(int j = 0; j < CONVN; ++j)
//			totalWeight += filter[i][j];
//
//	for(int i = 0; i < CONVN; ++i)
//		for(int j = 0; j < CONVN; ++j)
//			filter[i][j] /= totalWeight;

	MPI_Barrier(MPI_COMM_WORLD);

	double timeStart = MPI_Wtime();

	int ended = 0; // if loop should end
// Convolute Loop
	int t;
	for (t = 0; t < loops && !ended; ++t)
	{
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
	// Wait borders
		block.waitSendBorders();
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
/**/
// Print local finished matrices in worker order
//	if(LOG)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (int p=0; p<nWorkers; p++) {
			if (workerID == p) {
				printf("Local process on workerID %d is:\n", workerID);
				block.print();
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
/**/
// Gather matrix
	MPI_Gatherv(&block[1][1], 1,  block_t,
			&globalMatrix[0], counts, displs, submatrix_t,
			masterID, MPI_COMM_WORLD);
/**/
// Poll max timeElapsed
	double timeElapsed = MPI_Wtime() - timeStart;

	if(LOG && workerID == logWorkerID) printf("Reducing timeElapsed\n");

	double timeElapsedGlobal;
	MPI_Reduce(&timeElapsed, &timeElapsedGlobal, 1,
		MPI_DOUBLE, MPI_MAX, masterID, MPI_COMM_WORLD);

	if (workerID == masterID) {
	//  if(LOG)
		printf("timeElapsed %f\n", timeElapsedGlobal);
	}

/**/
// Print Matrix
//  if(LOG)
	if(workerID == masterID){
		for(int i = 0; i < height; ++i){
			std::cout << "|";
			for(int j = 0; j < width; ++j){
				std::cout << std::setw(2) << (int)globalMatrix[i*width + j] << "  ";
				// if((*this)[i][j] == 1)
				// 	std::cout << "O";
				// else
				// 	std::cout << ".";
			}
			std::cout << "|\n";
		}
		std::cout << std::flush;
	}
/**/
// Free memory
	MPI_Type_free(&column_t);
	MPI_Type_free(&row_t);
	MPI_Type_free(&block_t);
	MPI_Type_free(&submatrix_t);
/**/
	MPI_Finalize();

LIKWID_MARKER_CLOSE;

	return(EXIT_SUCCESS);
}
