#include <iostream>
#include <climits>
#include <cmath>
#include <time.h>
#include <fstream>
#include <random>
using namespace std;
//left = 2i + 1
//right = 2i + 2
//parent = (i-1)/2
struct heap 
{
	long long* elements; 
	int size; 
};
void max_heapify(long long* elements, int size, int n)
{
	int largest = n; 
	if(2*n + 1 < size && elements[2*n + 1] > elements[largest])
	{
		largest = 2*n + 1; 
	}
	if(2*n + 2 < size && elements[2*n + 2] > elements[largest])
	{
		largest = 2*n + 2; 
	}
	if(largest != n)
	{
		long long temp = elements[n]; 
		elements[n] = elements[largest]; 
		elements[largest] = temp; 
		max_heapify(elements, size, largest); 
	}
}
void build_heap(long long* a, int length)
{
	for(int i = length/2 - 1; i >= 0; i--)
	{
		max_heapify(a, length, i); 
	}
}
long long extract_max(heap* h)
{
	long long max = h->elements[0]; 
	h->elements[0] = h->elements[h->size - 1]; 
	h->size--; 
	max_heapify(h->elements, h->size, 0); 
	return max; 
}
void insert(heap* h, long long v)
{
	h->size++; 
	h->elements[h->size-1] = v;  
	int n = h->size-1; 
	while(n != 0 && h->elements[(n-1)/2] < h->elements[n])
	{
		long long temp = h->elements[(n-1)/2]; 
		h->elements[(n-1)/2] = h->elements[n]; 
		h->elements[n] = temp; 
		n = (n-1)/2; 
	}
}
long long karmakar_karp(long long* a, int length)
{
	heap h; 
	h.elements = a;
	h.size = length; 
	long long x; 
	long long y; 
	build_heap(h.elements, h.size); 
	while((h.size > 2 && (h.elements[1] != 0 || h.elements[2] != 0)) || (h.size == 2 && h.elements[1] != 0))
	{
		x = extract_max(&h);
		y = extract_max(&h);
		insert(&h, abs(x-y)); 
	}
	return h.elements[0];
}
void karmakar_karp_copy(long long* a, int length, long long* result)
{
	heap h; 
	h.elements = result;
	for (int i = 0; i < length; i++)
	{
		h.elements[i] = a[i];
	} 
	h.size = length; 
	long long x; 
	long long y; 
	build_heap(h.elements, h.size); 
	while((h.size > 2 && (h.elements[1] != 0 || h.elements[2] != 0)) || (h.size == 2 && h.elements[1] != 0))
	{
		x = extract_max(&h);
		y = extract_max(&h);
		insert(&h, abs(x-y)); 
	}
}
long long repeated_random(long long* a)
{
	long long best_residue = 0; 
	long long current_residue = 0; 
	//int* s = (int*)malloc(100*sizeof(int));
	for(int k = 0; k < 25000; k++)
	{
		current_residue = 0; 
		for(int i = 0; i < 100; i++)
		{
			if(((float)(rand()))/((float)(RAND_MAX)) < 0.5)
			{
				current_residue += a[i]; 
			}
			else
			{
				current_residue -= a[i];
			}
		}
		current_residue = abs(current_residue); 
		if(k == 1)
		{
			best_residue = current_residue; 
		}
		else if(current_residue < best_residue)
		{
			best_residue = current_residue;
		}
	}
	return best_residue;
}
long long calculate_residue(long long* a, int* s)
{
	long long residue = 0; 
	for(int i = 0; i < 100; i++)
	{
		residue += (a[i] * s[i]); 
	}
	return abs(residue); 
}
long long hill_climbing(long long* a)
{ 
	long long current_residue = 0; 
	long long new_residue = 0; 
	int* s = (int*)malloc(100*sizeof(int));
	int i; 
	int j; 
	bool swap_j;
	for(int k = 0; k < 25000; k++)
	{
		//current_residue = abs(current_residue); 
		if(k == 1)
		{
			for(int i = 0; i < 100; i++)
			{
				if(((float)(rand()))/((float)(RAND_MAX)) < 0.5)
				{
					current_residue += a[i]; 
					s[i] = 1; 
				}
				else
				{
					current_residue -= a[i];
					s[i] = -1;
				}
			}
		}
		else
		{
			i = (int)(rand() % 100); 
			j = (int)(rand() % 100);
			while(j == i)
			{
				j = rand() % 100;
			}
			if(((float)(rand()))/((float)(RAND_MAX)) < 0.5)
			{
				swap_j = true; 
			}
			else
			{
				swap_j = false; 
			}
			if(s[i] == 1)
			{
				if(swap_j && s[j] == 1)
				{
					if(abs(current_residue - 2*a[i] - 2*a[j]) < abs(current_residue))
					{
						current_residue = current_residue - 2*a[i] - 2*a[j];
						s[i] *= -1; 
						s[j] *= -1; 
					}
				}
				else if(swap_j && s[j] == -1)
				{
					if(abs(current_residue - 2*a[i] + 2*a[j]) < abs(current_residue))
					{
						current_residue = current_residue - 2*a[i] + 2*a[j];
						s[i] *= -1; 
						s[j] *= -1;
					}

				}
				else
				{
					if(abs(current_residue - 2*a[i]) < abs(current_residue))
					{
						current_residue = current_residue - 2*a[i];
						s[i] *= -1; 
						//s[j] *= -1; 
					}
				}
			}
			else
			{
				if(swap_j && s[j] == 1)
				{
					if(abs(current_residue + 2*a[i] - 2*a[j]) < abs(current_residue))
					{
						current_residue = current_residue + 2*a[i] - 2*a[j];
						s[i] *= -1; 
						s[j] *= -1;
					}
				}
				else if(swap_j && s[j] == -1)
				{
					if(abs(current_residue + 2*a[i] + 2*a[j]) < abs(current_residue))
					{
						current_residue = current_residue + 2*a[i] + 2*a[j];
						s[i] *= -1; 
						s[j] *= -1;
					}
				}
				else
				{
					if(abs(current_residue + 2*a[i]) < abs(current_residue))
					{
						current_residue = current_residue + 2*a[i];
						s[i] *= -1; 
					}
				}
			}
		}
	}
	return abs(current_residue);
}
double temp(long long iter)
{
	long long div = iter/300; 
	return (pow(10,10)*pow(0.8,div));
}
long long simulated_annealing(long long* a)
{
	long long best_residue = 0; 
	long long current_residue = 0; 
	long long new_residue = 0; 
	double power; 
	int* s = (int*)malloc(100*sizeof(int));
	int i; 
	int j; 
	bool swap_j;
	for(int k = 0; k < 25000; k++)
	{
		//current_residue = abs(current_residue); 
		if(k == 1)
		{
			for(int i = 0; i < 100; i++)
			{
				if(((float)(rand()))/((float)(RAND_MAX)) < 0.5)
				{
					current_residue += a[i]; 
					s[i] = 1; 
				}
				else
				{
					current_residue -= a[i];
					s[i] = -1;
				}
			}
			best_residue = abs(current_residue); 
		}
		else
		{
			i = rand() % 100; 
			j = rand() % 100;
			while(j == i)
			{
				j = rand() % 100;
			}
			if(((float)(rand()))/((float)(RAND_MAX)) < 0.5)
			{
				swap_j = true; 
			}
			else
			{
				swap_j = false; 
			}
			if(s[i] == 1)
			{
				if(swap_j && s[j] == 1)
				{
					power = -1 * (double)(abs(current_residue - 2*a[i] - 2*a[j])-abs(current_residue))/temp(k);
					if(abs(current_residue - 2*a[i] - 2*a[j]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue - 2*a[i] - 2*a[j];
						s[i] = -1; 
						s[j] = -1; 
					}
				}
				else if(swap_j && s[j] == -1)
				{
					power = -1 * (double)(abs(current_residue - 2*a[i] + 2*a[j])-abs(current_residue))/temp(k);
					if(abs(current_residue - 2*a[i] + 2*a[j]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue - 2*a[i] + 2*a[j];
						s[i] = -1; 
						s[j] = 1; 
					}

				}
				else
				{
					power = -1 * (double)(abs(current_residue - 2*a[i])-abs(current_residue))/temp(k);
					if(abs(current_residue - 2*a[i]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue - 2*a[i];
						s[i] = -1; 
					}
				}
			}
			else
			{
				if(swap_j && s[j] == 1)
				{
					power = -1 * (double)(abs(current_residue + 2*a[i] - 2*a[j])-abs(current_residue))/temp(k);
					if(abs(current_residue + 2*a[i] - 2*a[j]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue + 2*a[i] - 2*a[j];
						s[i] = 1; 
						s[j] = -1; 
					}
				}
				else if(swap_j && s[j] == -1)
				{
					power = -1 * (double)(abs(current_residue + 2*a[i] + 2*a[j])-abs(current_residue))/temp(k);
					if(abs(current_residue + 2*a[i] + 2*a[j]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue + 2*a[i] + 2*a[j];
						s[i] = 1; 
						s[j] = 1; 
					}
				}
				else
				{
					power = -1 * (double)(abs(current_residue + 2*a[i])-abs(current_residue))/temp(k);
					if(abs(current_residue + 2*a[i]) < abs(current_residue) || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
					{
						current_residue = current_residue + 2*a[i];
						s[i] = 1; 
					}
				}
			}
		}
		if(abs(current_residue) < best_residue)
		{
			best_residue = abs(current_residue);
		}
}
	return best_residue;
}
void partition_to_arr(long long* a_prime, long long* a, int* p)
{
	for(int i = 0; i < 100; i++)
	{
		a_prime[i] = 0; 
	}
	for(int i = 0; i < 100; i++)
	{
		a_prime[p[i]] += a[i]; 
	}
}
long long repeated_random_pre(long long* a)
{
	int* partition = (int*)malloc(sizeof(int)*100);
	long long* a_prime = (long long*)malloc(sizeof(long long)*100);
	long long best_residue;
	long long current_residue;
	for(int k = 0; k < 25000; k++)
	{
		for(int i = 0; i < 100; i++)
		{
			partition[i] = (int)(rand() % 100); 
		}
		partition_to_arr(a_prime, a, partition); 
		current_residue = karmakar_karp(a_prime,100); 
		if(k == 1 || current_residue < best_residue)
		{
			best_residue = current_residue; 
		}
	}
	return best_residue; 
}
long long simulated_annealing_pre(long long* a)
{
	int* p = (int*)malloc(sizeof(int)*100);
	int old_val; 
	int i; 
	int j; 
	long long* a_prime = (long long*)malloc(sizeof(long long)*100);
	long long* res = (long long*)malloc(sizeof(long long)*100);
	long long best_residue;
	long long new_residue;
	long long current_residue;
	double power; 
	for(int k = 0; k < 25000; k++)
	{
		if(k == 1)
		{
			for(int i = 0; i < 100; i++)
			{
				p[i] = rand() % 100; 
			}
			partition_to_arr(a_prime, a, p); 
			karmakar_karp_copy(a_prime,100,res); 
			best_residue = current_residue = res[0];
		}
		else
		{
			i = rand() % 100; 
			j = rand() % 100;
			while(j == p[i])
			{
				j = rand() % 100;
			}
			a_prime[p[i]] -= a[i]; 
			a_prime[j] += a[i]; 
			karmakar_karp_copy(a_prime,100,res); 
			new_residue = res[0];
			power = -1 * (double)(abs(new_residue - current_residue))/temp(k);
			if(new_residue < current_residue || ((float)(rand()))/((float)(RAND_MAX)) < exp(power))
			{
				p[i] = j; 
				current_residue = new_residue; 
			}
			else
			{
				a_prime[p[i]] += a[i]; 
				a_prime[j] -= a[i]; 
			}
			if(current_residue < best_residue)
			{
				best_residue = current_residue;
			}
		}
	}
	return best_residue; 
}
long long hill_climbing_pre(long long* a)
{
	int* p = (int*)malloc(sizeof(int)*100);
	int old_val; 
	int i; 
	int j; 
	long long* a_prime = (long long*)malloc(sizeof(long long)*100);
	long long* res = (long long*)malloc(sizeof(long long)*100);
	long long new_residue;
	long long best_residue;
	long long current_residue;
	for(int k = 0; k < 25000; k++)
	{
		if(k == 1)
		{
			for(int i = 0; i < 100; i++)
			{
				p[i] = rand() % 100; 
			}
			partition_to_arr(a_prime, a, p); 
			karmakar_karp_copy(a_prime,100,res); 
			current_residue = res[0]; 
		}
		else
		{
			i = rand() % 100; 
			j = rand() % 100;
			while(j == p[i])
			{
				j = rand() % 100;
			}
			a_prime[p[i]] -= a[i]; 
			a_prime[j] += a[i]; 
			karmakar_karp_copy(a_prime,100,res); 
			new_residue = res[0]; 
			if(new_residue < current_residue)
			{
				p[i] = j; 
				current_residue = new_residue; 
			}
			else
			{
				a_prime[p[i]] += a[i]; 
				a_prime[j] -= a[i]; 
			}
		}
	}
	return current_residue; 
}
int main(int argc, char **argv) 
{

   if(argc > 4)
   {
   	cout << "Too many arguments" << endl; 
   }
   else
   {
     srand(time(NULL));
     int version = atoi(argv[1]); 
     int alg = atoi(argv[2]); 
     char* filename;
     ifstream test_input; 
     clock_t begin_time;
     float time = 0; 
     if(argc == 4)
     {
     	filename = argv[3];
     	test_input.open(filename); 

     }
     long long* elements = (long long*)malloc(100*sizeof(long long)); 
     if(version != 0)
     {
     	random_device rd;
	 	mt19937 generator(rd());
     	uniform_int_distribution<long long> distribution(0,1000000000000);
     	long long total = 0; 
     	for(int k = 0; k < 100; k++)
     	{
     	 	for(int i = 0; i < 100; i++)
     	 	{
     			elements[i] = distribution(generator);   
     	 	} 
     	 	begin_time = clock();
     	 	karmakar_karp(elements,100);
     	 	time += float(clock () - begin_time)/CLOCKS_PER_SEC;
 	 	}
 	 	cout << "Total Time: " << (double)time/100.0 << endl;
 	}
 	else
 	{
 		for(int i = 0; i < 100; i++)
 		{
 			test_input >> elements[i];
 		}
 		if(alg == 0)
 		{
 			cout << karmakar_karp(elements,100); 
 		}
 		else if(alg == 1)
 		{
 			cout << repeated_random(elements); 
 		}
 		else if(alg == 2)
 		{
 			cout << hill_climbing(elements); 
 		}
 		else if(alg == 3)
 		{
 			cout << simulated_annealing(elements); 
 		}
 		else if(alg == 11)
 		{
 			cout << repeated_random_pre(elements);
 		}
 		else if(alg == 12)
 		{
 			cout << hill_climbing_pre(elements);
 		}
 		else if(alg == 13)
 		{
 			cout << simulated_annealing_pre(elements);
 		}
 	}
     
     return 0; 
	}
}