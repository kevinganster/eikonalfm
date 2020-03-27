// memory leak detection
//#include "debugtests.h"

#include "factoredmarcher.hpp"
#include <iostream>
#include <iomanip>
//#include <errcode.h>
#include "heap.cpp"

using namespace std;


int main()
{
	int ndim = 2;
	size_t size = 11 * 6;
	size_t* shape = new size_t[ndim] { 6, 11 };
	double* dx = new double[ndim]{ 0.8, 0.8 };
	int order = 2;

	/*double c_[] = {
		0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
		1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3,
		2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1,
		2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9,
		3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7,
		4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5
	};*/
	double c_[] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    };
	double* c = (double*)&c_;
	size_t x0 = 5;

	char str[21];
	sprintf(str, "a %f", 1.0f);
	cout << str << endl;

	double* tau = new double[size];
	Marcher* m = new Marcher(c, ndim, shape, dx, order);
	m->solve(x0, tau);

	int s = 0;
	cout << "marcher solution:" << endl;
	for (int i = 0; i < size; i++)
	{
		printf("%2.2f\t", tau[i]);
		if (++s >= shape[1])
		{
			s = 0;
			cout << endl;
		}
	}

	delete m;
	delete[] tau;


	tau = new double[size];
	m = new FactoredMarcher(c, ndim, shape, dx, order);
	m->solve(x0, tau);

	s = 0;
	cout << "factored marcher solution:" << endl;
	for (int i = 0; i < size; i++)
	{
		printf("%2.2f\t", tau[i]);
		if (++s >= shape[1])
		{
			s = 0;
			cout << endl;
		}
	}

	delete[] shape;
	delete[] dx;
	delete m;
	delete[] tau;

	//_CrtDumpMemoryLeaks();
	return 0;
}