
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
	initializeIntArray(localArray, localNumElements, level);
	printArray(localArray, localNumElements, taskid);

	for (int place = 0; place < 1; place++) {
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
		// int *recvTracker = localArray;
		// int *sendTracker = localSortedArray;
		int *recvOffset = new int[numtasks];
		recvOffset[numtasks-1] = 0;
		int *sendOffset = new int[numtasks];
		sendOffset[numtasks-1] = 0;
		// Recv from everyone bigger than me, send to everyone smaller than me
		for (int radix = 0; radix < 256; radix++) {
			if (radix > 0) {
				recvOffset[0] = recvOffset[numtasks-1] + sendRecvDist[radix-1][numtasks-1][taskid];
				sendOffset[0] = sendOffset[numtasks-1] + sendRecvDist[radix-1][taskid][numtasks-1];
			} else {
				recvOffset[0] = recvOffset[numtasks-1];
				sendOffset[0] = sendOffset[numtasks-1];
			}

			for (int proc = 1; proc < numtasks; proc++) {
				recvOffset[proc] = recvOffset[proc-1] + sendRecvDist[radix][proc-1][taskid];
			}
			
			for (int proc = 1; proc < numtasks; proc++) {
				sendOffset[proc] = sendOffset[proc-1] + sendRecvDist[radix][taskid][proc-1];
			}

			for (int proc = 0; proc < numtasks; proc++) {
				if (proc > taskid) {
					int numToReceive = sendRecvDist[radix][proc][taskid];
					if (numToReceive == 0) continue;
					MPI_Recv(localArray+recvOffset[proc], numToReceive, MPI_INT, proc, radix, MPI_COMM_WORLD, &status);
				} else if (proc < taskid) {
					int numToSend = sendRecvDist[radix][taskid][proc];
					if (numToSend == 0) continue;
					MPI_Send(localSortedArray+sendOffset[proc], numToSend, MPI_INT, proc, radix, MPI_COMM_WORLD);
				} else {
					int numToPlace = sendRecvDist[radix][taskid][taskid];
					if (numToPlace == 0) continue;
					for (int i = 0; i < numToPlace; i++) {
						int sendIdx = sendOffset[taskid] + i;
						int recvIdx = sendOffset[taskid] + i;
						localArray[recvIdx] = localSortedArray[sendIdx];
					}
				}
			}
		}
		
		recvOffset[numtasks-1] = 0;
		sendOffset[numtasks-1] = 0;
		// Recv from everyone smaller than me, send to everyone bigger than me
		for (int radix = 0; radix < 256; radix++) {
			if (radix > 0) {
				recvOffset[0] = recvOffset[numtasks-1] + sendRecvDist[radix-1][numtasks-1][taskid];
				sendOffset[0] = sendOffset[numtasks-1] + sendRecvDist[radix-1][taskid][numtasks-1];
			} else {
				recvOffset[0] = recvOffset[numtasks-1];
				sendOffset[0] = sendOffset[numtasks-1];
			}

			for (int proc = 1; proc < numtasks; proc++) {
				recvOffset[proc] = recvOffset[proc-1] + sendRecvDist[radix][proc-1][taskid];
			}
			
			for (int proc = 1; proc < numtasks; proc++) {
				sendOffset[proc] = sendOffset[proc-1] + sendRecvDist[radix][taskid][proc-1];
			}
			for (int proc = 0; proc < numtasks; proc++) {
				if (proc > taskid) {
					int numToSend = sendRecvDist[radix][taskid][proc];
					if (numToSend == 0) continue;
					MPI_Send(localSortedArray+sendOffset[proc], numToSend, MPI_INT, proc, radix, MPI_COMM_WORLD);
				} else if (proc < taskid) {
					int numToReceive = sendRecvDist[radix][proc][taskid];
					if (numToReceive == 0) continue;
					MPI_Recv(localArray+recvOffset[proc], numToReceive, MPI_INT, proc, radix, MPI_COMM_WORLD, &status);
				}
			}
		}
		delete[] recvOffset;
		recvOffset = nullptr;
		delete[] sendOffset;
		sendOffset = nullptr;

		// for (int proc = 0; proc < numtasks; proc++) {
		// 	if (proc > taskid) {
		// 		for (int radix = 0; radix < 256; radix++) {
		// 			int numToReceive = sendRecvDist[radix][proc][taskid];
		// 			if (numToReceive == 0) continue;
		// 			MPI_Recv(recvTracker, numToReceive, MPI_INT, proc, radix, MPI_COMM_WORLD, &status);
		// 			recvTracker += numToReceive;
		// 		}
		// 	} else if (proc < taskid) {
		// 		for (int radix = 0; radix < 256; radix++) {
		// 			int numToSend = sendRecvDist[radix][taskid][proc];
		// 			if (numToSend == 0) continue;
		// 			MPI_Send(sendTracker, numToSend, MPI_INT, proc, radix, MPI_COMM_WORLD);
		// 			sendTracker += numToSend;
		// 		}
		// 	} else {
		// 		for (int radix = 0; radix < 256; radix++) {
		// 			int numToPlace = sendRecvDist[radix][taskid][taskid];
		// 			if (numToPlace == 0) continue;
		// 			for (int i = 0; i < numToPlace; i++) {
		// 				*recvTracker = *sendTracker;
		// 				sendTracker++;
		// 				recvTracker++;
		// 			}
		// 		}
		// 	}
		// }

		// Recv from everyone smaller than me, send to everyone bigger than me
		// for (int proc = 0; proc < numtasks; proc++) {
		// 	if (proc < taskid) {
		// 		for (int radix = 0; radix < 256; radix++) {
		// 			int numToReceive = sendRecvDist[radix][proc][taskid];
		// 			if (numToReceive == 0) continue;
		// 			MPI_Recv(recvTracker, numToReceive, MPI_INT, proc, radix, MPI_COMM_WORLD, &status);
		// 			recvTracker += numToReceive;
		// 		}
		// 	} else if (proc > taskid) {
		// 		for (int radix = 0; radix < 256; radix++) {
		// 			int numToSend = sendRecvDist[radix][taskid][proc];
		// 			if (numToSend == 0) continue;
		// 			MPI_Send(sendTracker, numToSend, MPI_INT, proc, radix, MPI_COMM_WORLD);
		// 			sendTracker += numToSend;
		// 		}
		// 	}
		// }

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

	printArray(localArray, localNumElements, taskid);
	// std::cout << "Is Sorted: " << isSortedInt(localArray, localNumElements) << std::endl;

	if (localArray != nullptr) {
		delete[] localArray;
		localArray = nullptr;
	}
	
	MPI_Finalize();
	
}

void initializeIntArray(int *toSort, int numElements, SortLevel level) {
	int biggest_number = 100;
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
	std::cout << "Rank " << rank << " Array: " << std::endl;
	for (int i = 0; i < numElements; i++) {
		std::cout << arrayToPrint[i] << " " << std::flush;
	}
	std::cout << std::endl;
}

