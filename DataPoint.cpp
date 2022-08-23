#include "DataPoint.h"

//Default constructor
DataPoint::DataPoint()
{
}

//constructor
DataPoint::DataPoint(float* dimension, int old) :
	oldIndex(old), clusterId(-1), iscore(0)
{
	for (int i = 0; i<DIME_NUM; i++)
	{
		this->dimension[i] = dimension[i];
	}
}

//set the data dim
void DataPoint::SetDimension(float* dimension)
{
	for (int i = 0; i<DIME_NUM; i++)
	{
		this->dimension[i] = dimension[i];
	}
}

int DataPoint::GetOldIndex()
{
	return oldIndex;
}

//get the data dim
float* DataPoint::GetDimension()
{
	return this->dimension;
}

int DataPoint::GetClusterId()
{
	return this->clusterId;
}

void DataPoint::SetClusterId(int clusterId)
{
	this->clusterId = clusterId;
}

int DataPoint::isCore()
{
	return this->iscore;
}

void DataPoint::SetCore()
{
	this->iscore = 1;
}
