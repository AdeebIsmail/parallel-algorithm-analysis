
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

enum DataType {
	INT,
	DOUBLE
};

void initializeDoubleArray(double *toSort, int numElements, SortLevel level);
void initializeIntArray(int *toSort, int numElements, SortLevel level);
bool isSortedDouble(double *toSort, int numElements);
bool isSortedInt(int *toSort, int numElements);
void sequentialRadixSortDouble(double *toSort, int numElements);
void sequentialRadixSortInt(int *toSort, int numElements);


int main(int argc, char *argv[]) {
	if (argc != 5) {
		printf("Usage: ./mpi_radix <data_type> <num_elements> <num_procs> <sort_level>");
	}
	
	DataType type = static_cast<DataType>(std::stoi(argv[1]));
	int numElements = std::stoi(argv[2]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[4]));

	switch (type)
	{
	case INT:
		{
			int* toSort = new int[numElements];
			initializeIntArray(toSort, numElements, level);
			// for (int i = 0; i < numElements; i++) {
			// 	std::cout << toSort[i] << " ";
			// }
			std::cout << std::endl << " issorted: " << isSortedInt(toSort, numElements) << std::endl;
			sequentialRadixSortInt(toSort, numElements);
			delete[] toSort;
			
		}
		break;
	case DOUBLE:
		{
			double* toSort = new double[numElements];
			initializeDoubleArray(toSort, numElements, level);
			delete[] toSort;
		}
		break;
	default:
		break;
	}

}


void initializeDoubleArray(double *toSort, int numElements, SortLevel level) {

	double biggest_number = 2e6;
	double max_increment = biggest_number/(double)(numElements);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> small_numbers(0, 1);
	std::uniform_real_distribution<> increment(0, max_increment);
	std::uniform_int_distribution<> indexes(0, numElements-1);

	switch (level)
	{
	case SORTED:
		// initialize a sorted array
		for (int i = 0; i < numElements; i++) {
			if (i != 0) {
				toSort[i] = toSort[i-1]+increment(gen);
			} else {
				toSort[i] = small_numbers(gen);
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
					toSort[i] = small_numbers(gen);
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
			toSort[i] = biggest_number*small_numbers(gen);
		}
		break;
	case REVERSED:
		for (int i = numElements-1; i >= 0; i--) {
			if (i!=numElements-1) {
				toSort[i] = toSort[i+1]+increment(gen);
			} else {
				toSort[i] = small_numbers(gen);
			}
		}
		break;
	default:
		break;
	}

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


bool isSortedDouble(double *toSort, int numElements) {

   bool sorted = true;
   
   for (int i = 1; i < numElements-1; i++) {
       if (toSort[i] > toSort[i+1]) {
          sorted = false;
       }
   }
    
   return sorted;
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

void sequentialRadixSortDouble(double *toSort, int numElements) {
	
}

void sequentialRadixSortInt(int *toSort, int numElements) {
	std::cout << "Size of int: " << sizeof(int) << std::endl;
	
}



