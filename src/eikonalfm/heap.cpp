// memory leak detection
//#ifdef _DEBUG
//	#include "debugtests.h"
//#endif

#include <vector>
#include <algorithm>

//#include <iostream>
//using namespace std;


template<typename _Comp>
class Heap
{
protected:
	std::vector<unsigned long> c;	// underlying container
	_Comp comp;				        // comparator functor
    unsigned long *m;				// mapping from position to index in the container

public:
	explicit Heap(const _Comp& _Comparator, const unsigned long max_size) :
		c(), comp(_Comparator), m(new unsigned long[max_size])
	{	// construct with specified comparator
	}

	~Heap()
	{
		delete[] m;
	}
	
	bool empty() const
	{	// test if the heap is empty
		return c.empty();
	}

    unsigned long size() const
	{	// return length of queue
		return (unsigned long)c.size();
	}

    unsigned long top() const
	{
		return c.front();
	}

	void push(const unsigned long item)
	{
		// append to the end of the container
		c.push_back(item);
		// save mapping
		m[item] = size() - 1;
		// rebuild heap property from bottom to top
		sift_up(size() - 1);
	}

	void pop()
	{
		m[c.back()] = 0;
		// swap root with the last element and remove it from the container
		std::swap(c.front(), c.back());
		c.pop_back();
		// rebuild the heap property
		sift_down(0);
	}

	void update(const unsigned long item)
	{
        unsigned long start = m[item];
		sift_up(start);
		sift_down(start);
	}

	//void print()
	//{
	//	cout << "container" << endl << "[";
	//	for (unsigned long i : c)
	//		cout << i << ", ";
	//	cout << "]" << endl;

	//	//std::cout << "mappings" << std::endl;
	//	//for (std::map<_ItemType, unsigned long>::iterator i=m.begin(); i != m.end(); i++)
	//	//	std::cout << i->first << ":" << i->second << ", ";
	//	//std::cout << std::endl;
	//}

private:
	inline unsigned long parent(const unsigned long i) const
	{
		return (i - 1) >> 1;
	}

	inline unsigned long left(const unsigned long i) const
	{
		return (i << 1) + 1;
	}

	inline unsigned long right(const unsigned long i) const
	{
		return (i << 1) + 2;
	}

	void sift_down(const unsigned long start)
	{
        unsigned long s = size();
        unsigned long i = start;

		while (true)
		{
            unsigned long min = i;
			// left child's value is smallest -> needs swap
			if (left(i) < s && comp(c[left(i)], c[min]))
				min = left(i);
			// right child's value is smallest -> needs swap
			if (right(i) < s && comp(c[right(i)], c[min]))
				min = right(i);
			// heap condition restored -> stop sifting
			if (min == i)
				break;

			std::swap(m[c[i]], m[c[min]]);
			std::swap(c[i], c[min]);
			i = min;
		}
	}

	void sift_up(const unsigned long start)
	{
        unsigned long i = start;
		while (true)
		{
			if (i == 0)
				break;

            unsigned long p = parent(i);
			// parent value is 'larger' than i's -> needs swap
			if (comp(c[i], c[p]))
			{
				std::swap(m[c[i]], m[c[p]]);
				std::swap(c[i], c[p]);
				i = p;
			}
			// heap condition restored -> stop sifting
			else
				break;
		}
	}
};


//int main()
//{
//	double* tau = new double[6]{ 1, 2, 3, 4, 5, 6 };
//	auto heap_comp = [&tau](const long e1, const long e2) { return tau[e1] < tau[e2]; };
//	Heap<long, decltype(heap_comp)> h(heap_comp);
//
//	h.push(0);
//	h.push(1);
//	h.push(2);
//
//	std::cout << "brefore update" << std::endl;
//	h.print();
//
//	tau[1] = 0;
//	h.update(1);
//	
//	std::cout << std::endl;
//	std::cout << "after update" << std::endl;
//	h.print();
//
//
//	std::cout << "size: " << h.size() << std::endl;
//	while (!h.empty())
//	{
//		std::cout << h.top() << ":" << tau[h.top()] << ", ";
//		h.pop();
//	}
//	std::cout << std::endl;
//
//	return 0;
//}
