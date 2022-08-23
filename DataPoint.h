#pragma once
#include <cmath>
using namespace std;

const int DIME_NUM = 5+1;        //set the DIME_NUME as d+1

class DataPoint
{
private:
	float dimension[DIME_NUM];  // the coordinate of this point
	int oldIndex; // the index of this point in the original file
	int clusterId;                    //the cluster ID
	int iscore; // indicate whether this point a core point
public:

	DataPoint();
	DataPoint(float* dimension, int old);

	float* GetDimension();
	void SetDimension(float* dimension);
	int GetOldIndex();
	int GetClusterId();
	void SetClusterId(int classId);
	int isCore();
	void SetCore();
};

