// memory leak detection
#include "debugtests.h"

#include "factoredmarcher.h"
#include <iostream>
#include <iomanip>
#include <errcode.h>
#include "heap.cpp"

using namespace std;


int main()
{
	int ndim = 2;
	size_t size = 11 * 11;
	long* shape = new long[ndim] { 11, 11 };
	double* dx = new double[ndim]{ 0.4, 0.7 };
	int order = 2;

	double c_[] = {
		0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
	};
	double* c = (double*)&c_;
	size_t x0 = 0;

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