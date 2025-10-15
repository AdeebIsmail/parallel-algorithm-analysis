
#include "mpi.h"
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <random>
#include <string>
#include <iomanip>

enum SortLevel {
	SORTED,
	PERTURBED,
	RANDOM,
	REVERSED 
};

void initializeIntArray(int *toSort, int numElements, SortLevel level);
bool isSortedInt(int *toSort, int numElements);
void sequentialRadixSortInt(int *&toSort, int numElements);


int main(int argc, char *argv[]) {
	if (argc != 3) {
		printf("Usage: ./mpi_radix <num_elements> <sort_level>");
		return 1;
	}
	MPI_Init(&argc, &argv);
	
	int numElements = std::stoi(argv[1]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[2]));

	int* toSort = new int[numElements];
	initializeIntArray(toSort, numElements, level);
	std::cout << "Before: ";
	for (int i = 0; i < numElements; i++) {
		std::cout << toSort[i] << " ";
	}
	std::cout << std::endl;
	sequentialRadixSortInt(toSort, numElements);
	std::cout << " issorted: " << isSortedInt(toSort, numElements) << std::endl;
	delete[] toSort;
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



void sequentialRadixSortInt(int *&toSort, int numElements) {
	// grace is little endian
	// std::cout << "Size of int: " << sizeof(int) << std::endl;
	// int testInt = 65536;
	// unsigned char *bytes = (unsigned char*)(&testInt);
	// printf("First byte: %i", *bytes);
	// printf("Second byte: %i", *(bytes+1));
	// printf("Third byte: %i", *(bytes+2));
	// printf("Fourth byte: %i", *(bytes+3));
	
	
	// Sort byte by byte
	for (int place = 0; place < 4; place++) {
		int *sortedArray = new int[numElements];
		// build the histogram
		int histogram[256];
		memset(histogram, 0, sizeof(int)*256);

		for (int i = 0; i < numElements; i++) {
			// peel the radix
			unsigned char* bytePtr = (unsigned char *)(&toSort[i]);
			unsigned char radix = *(bytePtr + place);
			histogram[radix]++;
		}

		// build offset table, each index corresponds to 
		// radix value
		int offsetIdx[256];
		offsetIdx[0] = 0;
		for (int i = 1; i < 256; i++) {
			offsetIdx[i] = offsetIdx[i-1] + histogram[i-1];
		}
		
		// swap each element into correct place
		for (int i = 0; i < numElements; i++) {
			unsigned char* bytePtr = (unsigned char *)(&toSort[i]);
			unsigned char radix = *(bytePtr + place);
			int correctIdx = offsetIdx[radix]++;
			sortedArray[correctIdx] = toSort[i];
		}
		// figure out mem leak
		delete[] toSort;
		toSort = sortedArray;
	}

	std::cout << "After: ";
	for (int i = 0; i < numElements; i++) {
		std::cout << toSort[i] << " ";
	}
	std::cout << std::endl;
}



