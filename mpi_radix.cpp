
#include "mpi.h"
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <random>
#include <string>
#include <iomanip>
#include <cassert>

#define GlOBAL_MAX_NUMBER 400e6

enum SortLevel {
	SORTED,
	PERTURBED,
	RANDOM,
	REVERSED 
};

void initializeIntArray(int *arrayToInit, int numElements, SortLevel level, int taskid, int numtasks);
bool isLocallySorted(int *arrayToCheck, int numElements);
void printArray(int *arrayToPrint, int numElements, int rank);

int main(int argc, char *argv[]) {

	if (argc != 3) {
		printf("Usage: ./mpi_radix <num_elements> <sort_level>");
		return 1;
	}

	int totalNumElements = std::stoi(argv[1]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[2]));

	// communication variables
	int numtasks,
		taskid,
		localNumElements;

	int *localArray;
	int **globalHistogram;
	int **radixDist;
	int ***sendRecvDist; // (radix, send, recv)
	int *localOffsetIdx;
	int *localSortedArray;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (totalNumElements%numtasks != 0) {
		printf("Error: num_procs indivisible by num_elements");
		int rc;
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
	// init the local array
	localNumElements = totalNumElements/numtasks;
	localArray = new int[localNumElements];
	initializeIntArray(localArray, localNumElements, level, taskid, numtasks);

	for (int place = 0; place < 4; place++) {
		localSortedArray = new int[localNumElements];
		
		// init global histogram
		globalHistogram = new int*[numtasks];
		radixDist = new int*[numtasks];
		for (int i = 0; i < numtasks; i++) {
			globalHistogram[i] = new int[256];
			radixDist[i] = new int[256];
			memset(globalHistogram[i], 0, sizeof(int)*256);
			memset(radixDist[i], 0, sizeof(int)*256);
		}

		// init send/recv histogram
		sendRecvDist = new int**[256];
		for (int i = 0; i < 256; i++) {
			sendRecvDist[i] = new int*[numtasks];
			for (int j = 0; j < numtasks; j++) {
				sendRecvDist[i][j] = new int[numtasks];
				memset(sendRecvDist[i][j], 0, sizeof(int)*numtasks);
			}
		}

		// Build local histogram
		int* localHist = globalHistogram[taskid];
		for (int i = 0; i < localNumElements; i++) {
			// peel the radix
			unsigned char* bytePtr = (unsigned char *)(&localArray[i]);
			unsigned char radix = *(bytePtr + place);
			localHist[radix]++;
		}

		// Create local offset table
		localOffsetIdx = new int[256];
		localOffsetIdx[0] = 0;
		for (int i = 1; i < 256; i++) {
			localOffsetIdx[i] = localOffsetIdx[i-1] + localHist[i-1];
		}

		// put local elements in right place
		for (int i = 0; i < localNumElements; i++) {
			unsigned char* bytePtr = (unsigned char *)(&localArray[i]);
			unsigned char radix = *(bytePtr + place);
			int correctIdx = localOffsetIdx[radix]++;
			localSortedArray[correctIdx] = localArray[i];
		}

		// broadcast histogram to everyone
		for (int i = 0; i < numtasks; i++) {
			MPI_Bcast(globalHistogram[i], 256, MPI_INT, i, MPI_COMM_WORLD);
		}

		// reduce the histogram
		int* histogramSum = new int[256];
		for (int i = 0; i < 256; i++) {
			int count = 0;
			for (int proc = 0; proc < numtasks; proc++) {
				count += globalHistogram[proc][i];
			}
			histogramSum[i] = count;
		}

		// Calculate how many of each radix each process receives
		int currCount = 0;
		int currProc = 0;
		for (int i = 0; i < 256; i++) {
			while (histogramSum[i] !=0) {
				int capacity = localNumElements - currCount;
				int load = std::min(capacity, histogramSum[i]);
				radixDist[currProc][i]+=load;
				currCount+=load;
				histogramSum[i] -= load;
				if (currCount == localNumElements) {
					currProc++;
					currCount = 0;
				}
			}
		}

		int recvProc = 0;
		// Calculate how many of each radix each process sends to p1 ... pn
		for (int radix = 0; radix < 256; radix++) {
			recvProc = 0;
			for (int sendProc = 0; sendProc < numtasks; sendProc++) {
				while (globalHistogram[sendProc][radix] != 0) {
					int numToSend = std::min(radixDist[recvProc][radix], globalHistogram[sendProc][radix]);
					radixDist[recvProc][radix] -= numToSend;
					globalHistogram[sendProc][radix] -= numToSend;
					sendRecvDist[radix][sendProc][recvProc] = numToSend;
					if (radixDist[recvProc][radix] == 0) {
						recvProc++;
					}
				}
			}
		}

		MPI_Status status;
		int *recvTracker = localArray;
		int *sendTracker = localSortedArray;

		// Recv from everyone smaller than me, then send to everyone smaller than me
		// Send to everyone bigger than me, then recv from everyone bigger than me
		for (int radix = 0; radix < 256; radix++) {
			for (int pivotProc = 0; pivotProc < numtasks; pivotProc++) {
				int numToSend = sendRecvDist[radix][taskid][pivotProc];
				int numToRecv = sendRecvDist[radix][pivotProc][taskid];
				if (pivotProc > taskid) {
					// send first, then recv
					if (numToSend > 0) {
						MPI_Send(sendTracker, numToSend, MPI_INT, pivotProc, radix, MPI_COMM_WORLD);
						sendTracker += numToSend;
					}
					if (numToRecv > 0) {
						MPI_Recv(recvTracker, numToRecv, MPI_INT, pivotProc, radix, MPI_COMM_WORLD, &status);
						recvTracker += numToRecv;
					}
				} else if (pivotProc < taskid) {
					// recv first, then send
					if (numToRecv > 0) {
						MPI_Recv(recvTracker, numToRecv, MPI_INT, pivotProc, radix, MPI_COMM_WORLD, &status);
						recvTracker += numToRecv;
					}
					if (numToSend > 0) {
						MPI_Send(sendTracker , numToSend, MPI_INT, pivotProc, radix, MPI_COMM_WORLD);
						sendTracker += numToSend;
					}
				} else {
					if (numToSend > 0) {
						std::memcpy(recvTracker, sendTracker, numToSend*sizeof(int));
						recvTracker += numToRecv;
						sendTracker += numToSend;
					} 
				}
			}
		}

		if (localSortedArray != nullptr) {
			delete[] localSortedArray;
			localSortedArray = nullptr;
		}

		if (globalHistogram != nullptr) {
			for (int i = 0; i < numtasks; i++) {
				delete[] globalHistogram[i];
				globalHistogram[i] = nullptr;
			}
			delete[] globalHistogram;
			globalHistogram = nullptr;
		}

		if (radixDist != nullptr) {
			for (int i = 0; i < numtasks; i++) {
				delete[] radixDist[i];
				radixDist[i] = nullptr;
			}
			delete[] radixDist;
			radixDist = nullptr;
		}


		if (sendRecvDist != nullptr) {
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < numtasks; j++) {
					delete[] sendRecvDist[i][j];
					sendRecvDist[i][j] = nullptr;
				}
				delete[] sendRecvDist[i];
				sendRecvDist[i] = nullptr;
			}
			delete[] sendRecvDist;
			sendRecvDist = nullptr;
		}

		if (localOffsetIdx != nullptr) {
			delete[] localOffsetIdx;
			localOffsetIdx = nullptr;

		}

		if (histogramSum != nullptr) {
			delete[] histogramSum;
		}

	}

	// check that the array is globally sorted
	int *firstElements = new int[numtasks];
	firstElements[taskid] = localArray[0];

	for (int proc = 0; proc < numtasks; proc++) {
		MPI_Bcast(&firstElements[proc], 1, MPI_INT, proc, MPI_COMM_WORLD);
	}
	
	assert(isLocallySorted(localArray, localNumElements));
	if (taskid != numtasks-1) {
		bool correctOrder = localArray[localNumElements-1] <= firstElements[taskid+1];
		assert(correctOrder);
	}

	if (localArray != nullptr) {
		delete[] localArray;
		localArray = nullptr;
	}

	if (firstElements != nullptr) {
		delete[] firstElements;
		firstElements = nullptr;
	}
	
	MPI_Finalize();
}

void initializeIntArray(int *arrayToInit, int numElements, SortLevel level, int taskid, int numtasks) {
	
	int localMaxNumber = ceil(GlOBAL_MAX_NUMBER * ((double)(taskid+1))/((double)numtasks));
	int localMinNumber;

	if (taskid > 0) {
		int prevMaxNumber = ceil(GlOBAL_MAX_NUMBER * ((double)(taskid))/((double)numtasks));
		localMinNumber = prevMaxNumber;
	} else {
		localMinNumber = 0;
	}

	if (level == REVERSED) {
		localMaxNumber = ceil(GlOBAL_MAX_NUMBER * ((double)(numtasks-taskid))/((double)numtasks));
		localMinNumber = ceil(GlOBAL_MAX_NUMBER * ((double)(numtasks-taskid-1))/((double)numtasks));
	}

	int maxIncrement = floor((double)(localMaxNumber-localMinNumber)/(double)(numElements));

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> random_numbers(0, GlOBAL_MAX_NUMBER);
	std::uniform_int_distribution<> increment(0, maxIncrement);
	std::uniform_int_distribution<> indexes(0, numElements-1);

	switch (level)
	{
	case SORTED:
		// initialize a sorted array
		for (int i = 0; i < numElements; i++) {
			if (i != 0) {
				arrayToInit[i] = arrayToInit[i-1]+increment(gen);
			} else {
				arrayToInit[i] = localMinNumber;
			}
		}
		break;
	
	case PERTURBED: 
		{
			// initialize a sorted array
			for (int i = 0; i < numElements; i++) {
				if (i != 0) {
					arrayToInit[i] = arrayToInit[i-1]+increment(gen);
				} else {
					arrayToInit[i] = localMinNumber;
				}
			}
			// randomly swap 1% of all indices
			double swap_percentage = 0.01;
			int swapCount = ceil((double)(numElements) * swap_percentage);

			for (int i = 0; i < swapCount; i++) {
				int indx1 = indexes(gen);
				int indx2 = indexes(gen);

				while (indx1 == indx2 && numElements > 1) {
					indx2 = indexes(gen);
				}

				std::swap(arrayToInit[indx1], arrayToInit[indx2]);
			}
		}
		break;
	case RANDOM:
		// initialize a randomly sorted array
		for (int i = 0; i < numElements; i++) {
			arrayToInit[i] = random_numbers(gen);
		}
		break;
	case REVERSED:
	{
		for (int i = numElements-1; i >= 0; i--) {
			if (i!=numElements-1) {
				arrayToInit[i] = arrayToInit[i+1]+increment(gen);
			} else {
				arrayToInit[i] = localMinNumber;
			}
		}
	}

		break;
	default:
		break;
	}
}


bool isLocallySorted(int *arrayToCheck, int numElements) {

	bool sorted = true;
	
	for (int i = 1; i < numElements-1; i++) {
		if (arrayToCheck[i] > arrayToCheck[i+1]) {
		   sorted = false;
		}
	}
	 
	return sorted;
 }

void printArray(int *arrayToPrint, int numElements, int rank) {
	std::cout << "Rank " << rank << " Array: " << std::endl;
	for (int i = 0; i < numElements; i++) {
		std::cout << arrayToPrint[i] << " " << std::flush;
	}
	std::cout << std::endl;
}

