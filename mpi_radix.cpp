
#include "mpi.h"
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <random>
#include <string>
#include <iomanip>

#define MASTER 0

enum SortLevel {
	SORTED,
	PERTURBED,
	RANDOM,
	REVERSED 
};

void initializeIntArray(int *toSort, int numElements, SortLevel level);
bool isSortedInt(int *toSort, int numElements);
void printArray(int *arrayToPrint, int numElements, int rank);

int main(int argc, char *argv[]) {

	if (argc != 3) {
		printf("Usage: ./mpi_radix <num_elements> <sort_level>");
		return 1;
	}

	int totalNumElements = std::stoi(argv[1]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[2]));
	int* toSort;
	int* sortedArray;

	// communication variables
	int numtasks,
		taskid,
		localNumElements;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (totalNumElements%numtasks != 0) {
		printf("Error: num_procs indivisible by num_elements");
		int rc;
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}

	localNumElements = totalNumElements/numtasks;
	int *recvElem = new int[localNumElements];

	if (taskid == MASTER) {
		// initialization code
		toSort = new int[totalNumElements];
		initializeIntArray(toSort, totalNumElements, level);
		// printArray(toSort, totalNumElements, -1);
	} else{
		// do nothing for now I guess
	}

	// Implement the algorithm here
	for (int place = 0; place < 4; place++) {
		sortedArray = new int[totalNumElements];
		MPI_Scatter(toSort, localNumElements, MPI_INT, 
					recvElem, localNumElements, MPI_INT, 
					MASTER, MPI_COMM_WORLD);
		
		// construct the local histogram
		int histogram[256];
		memset(histogram, 0, sizeof(int)*256);
		for (int i = 0; i < localNumElements; i++) {
			// peel the radix
			unsigned char* bytePtr = (unsigned char *)(&recvElem[i]);
			unsigned char radix = *(bytePtr + place);
			histogram[radix]++;
		}
		
		// construct the global histogram
		if (taskid == MASTER) {
			MPI_Reduce(MPI_IN_PLACE, &histogram, 256, MPI_INT, MPI_SUM , MASTER, MPI_COMM_WORLD);
		} else {
			MPI_Reduce(&histogram, nullptr, 256, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
		}

		// build global offset table, each index corresponds to 
		// radix value
		if (taskid == MASTER) {

			int offsetIdx[256];
			offsetIdx[0] = 0;
			
			for (int i = 1; i < 256; i++) {
				offsetIdx[i] = offsetIdx[i-1] + histogram[i-1];
			}

			// put the items in the right place
			for (int i = 0; i < totalNumElements; i++) {
				unsigned char* bytePtr = (unsigned char *)(&toSort[i]);
				unsigned char radix = *(bytePtr + place);
				int correctIdx = offsetIdx[radix]++;
				sortedArray[correctIdx] = toSort[i];
			}
			delete[] toSort;
			toSort = sortedArray;
		}

		// pseudocode to calculate sending elements correctly
		// Build local histogram
		// MPI_Bcast histogram to everyone
		// Calculate how many of each radix each process receives
		// Calculate how many of each radix each process sends to p1 ... pn
		// Calculate how many of each radix each process receives from p1 ... pn
		// Recv from everyone smaller than me, send to everyone bigger than me
		// Recv from everyone bigger than me, send to everyone smaller than me
	}
	if (taskid == MASTER) {
		std::cout << "Total num elements: " << totalNumElements << std::endl;
		std::cout << "Array is sorted: " << isSortedInt(sortedArray, totalNumElements) << std::endl;
	}
	
	if (recvElem != nullptr) {
		delete[] recvElem;
		recvElem = nullptr;
	}
	if (taskid == MASTER && toSort!=nullptr) {
		delete[] toSort;
	}

	MPI_Finalize();
	
}

void initializeIntArray(int *toSort, int numElements, SortLevel level) {
	int biggest_number = 2e9;
	int max_increment = floor((double)(biggest_number)/(double)(numElements));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> random_numbers(0, biggest_number);
	std::uniform_int_distribution<> increment(0, max_increment);
	std::uniform_int_distribution<> indexes(0, numElements-1);

	switch (level)
	{
	case SORTED:
		// initialize a sorted array
		for (int i = 0; i < numElements; i++) {
			if (i != 0) {
				toSort[i] = toSort[i-1]+increment(gen);
			} else {
				toSort[i] = 0;
			}
		}
		break;
	
	case PERTURBED: 
		{
			// initialize a sorted array
			for (int i = 0; i < numElements; i++) {
				if (i != 0) {
					toSort[i] = toSort[i-1]+increment(gen);
				} else {
					toSort[i] = 0;
				}
			}
			// randomly swap 1% of all indices
			double swap_percentage = 0.01;
			int swapCount = ceil((double)(numElements) * swap_percentage);

			for (int i = 0; i < swapCount; i++) {
				int indx1 = indexes(gen);
				int indx2 = indexes(gen);

				while (indx1 == indx2) {
					indx2 = indexes(gen);
				}

				std::swap(toSort[indx1], toSort[indx2]);
			}
		}
		break;
	case RANDOM:
		// initialize a randomly sorted array
		for (int i = 0; i < numElements; i++) {
			toSort[i] = random_numbers(gen);
		}
		break;
	case REVERSED:
		for (int i = numElements-1; i >= 0; i--) {
			if (i!=numElements-1) {
				toSort[i] = toSort[i+1]+increment(gen);
			} else {
				toSort[i] = 0;
			}
		}
		break;
	default:
		break;
	}
}


bool isSortedInt(int *toSort, int numElements) {

	bool sorted = true;
	
	for (int i = 1; i < numElements-1; i++) {
		if (toSort[i] > toSort[i+1]) {
		   sorted = false;
		}
	}
	 
	return sorted;
 }

void printArray(int *arrayToPrint, int numElements, int rank) {
	std::cout << "Rank " << rank << " Array: ";
	for (int i = 0; i < numElements; i++) {
		std::cout << arrayToPrint[i] << " ";
	}
	std::cout << std::endl;
}

