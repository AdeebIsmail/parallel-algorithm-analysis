
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

void initialize_array(double *toSort, int numElements, char *dataType, SortLevel level);

int main(int argc, char *argv[]) {
	if (argc != 5) {
		printf("Usage: ./mpi_radix <data_type> <num_elements> <num_procs> <sort_level>");
	}
	
	int numElements = std::stoi(argv[2]);
	SortLevel level = static_cast<SortLevel>(std::stoi(argv[4]));
	double* toSort = new double[numElements];
	initialize_array(toSort, numElements, argv[1], level);
	
	// for (int i = 0; i < numElements; i++) {
	// 	std::cout << toSort[i] << std::endl;
	// }
	// std::cout << toSort[numElements-1] << " " << toSort[0] << std::endl;
	std::cout << sizeof(double) << " " << sizeof(int) << std::endl;	
	delete[] toSort;

}

void initialize_array(double *toSort, int numElements, char *dataType, SortLevel level) {

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


