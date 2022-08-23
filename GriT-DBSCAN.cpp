#include "ClusterAnalysis.h"

int main()
{
	ClusterAnalysis myClusterAnalysis;       //Clustering algorithm object declaration.
	myClusterAnalysis.Init("data\\farm.txt", 2000, 100);      //Algorithm initialization.
	printf("clusting the data...\n");
	double start = getCurrentTime();
	myClusterAnalysis.DoDBSCANRecursive();                    //Perform GriT-DBSCAN.
	start = getCurrentTime() - start;
	myClusterAnalysis.printMessage();
	printf("the running time is %.4f\n", start);
    //myClusterAnalysis.WriteToFile("data\\result.txt");      //Save the result.
	system("pause");

	return 0;
}
