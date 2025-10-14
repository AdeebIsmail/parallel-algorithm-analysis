
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

void initializeArray(double *toSort, int numElements, char *dataType, SortLevel level);
bool isSorted(double *toSort, int numElements);


int main(int argc, char *argv[]) {
	if (argc != 5) {
		printf("Usage: ./mpi_radix <data_type> <num_elements> <num_procs> <sort_level>");
	}
	
	int numElements = std::stoi(argv[2]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[4]));
	double* toSort = new double[numElements];
	initializeArray(toSort, numElements, argv[1], level);
        std::cout << "Is sorted: " << isSorted(toSort, numElements) << std::endl;
        for (int i = 0; i < numElements; i++) {
	    std::cout << toSort[i] << " ";
	
	}	
	delete[] toSort;

}


void initializeArray(double *toSort, int numElements, char *dataType, SortLevel level) {

	double biggest_number = 10;
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
bool isSorted(double *toSort, int numElements) {

   bool sorted = true;
   
   for (int i = 1; i < numElements-1; i++) {
       if (toSort[i] > toSort[i+1]) {
          sorted = false;
       }
   }
    
   return sorted;
}    		



