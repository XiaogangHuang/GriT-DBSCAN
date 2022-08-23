#pragma once
#include <vector>
#include "Toolfunction.h"

class GridCell {
protected:
	int gridCoords[DIME_NUM - 1]; //The coordinates of this grid cell
	std::vector<int> points; //The points in this grid cell
	int size;
	int clust;
	int coreNum;
	int diffPos;
public:
	GridCell() {}
	GridCell(int* coord, vector<int>& a, int s, int d) :
		points(a), size(s), clust(-1), coreNum(0), diffPos(d)
	{
		for (int i = 0; i < DIME_NUM-1; i++)
		{
			this->gridCoords[i] = coord[i];
		}
	}

	std::vector<int>& getPoints(); // return the points in this grid cell
	int* getCoords(); // return the coordinates of this grid cell
	int getSize(); // get the size of this grid cell
	void setClust(int c);	// set the cluster id for this grid cell
	int getClust();	// return the cluster id for this grid cell
	void setCoreNum(int c);	// set the number of core points in this grid cell
	void incCoreNum();	// increment the number of core points in this grid cell
	int getCoreNum();	// return the number of core points in this grid cell
	int getDiffPos();	// return the position of different coordinate 
};