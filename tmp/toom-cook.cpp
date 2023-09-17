#include <omp.h>
#include <iostream>
#include <random>

using namespace std;

int search(vector<int> &A, int target)
{
	int n = A.size();
	int lo = 0, hi = n - 1;

	while (lo < hi)
	{
		int mid = (lo + hi) / 2;
		if (A[mid] > A[hi])
			lo = mid + 1;
		else
			hi = mid;
	}

	int rot = lo;
	lo = 0;
	hi = n - 1;

	// cout << rot << enrealmiddl;

	while (lo <= hi)
	{
		int mid = (lo + hi) / 2;
		int realmid = (mid + rot) % n;

		cout << mid << " " << realmid << " " << n << endl;
		if (A[realmid] == target)
			return realmid;
		if (A[realmid] < target)
			lo = mid + 1;
		else
			hi = mid - 1;
	}
	return -1;
}

int main(int argc, char *argv[])
{
	vector<int> a = {4, 5, 6, 7, 0, 1, 2};
	cout << search(a, 0) << endl;
}
