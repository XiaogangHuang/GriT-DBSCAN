#include "GridCell.h"

std::vector<int>& GridCell::getPoints()
{
	return points;
}

int* GridCell::getCoords()
{
	return gridCoords;
}

int GridCell::getSize()
{
	return this->size;
}

void GridCell::setClust(int c)
{
	this->clust = c;
}

int GridCell::getClust()
{
	return this->clust;
}


void GridCell::setCoreNum(int c)
{
	this->coreNum = c;
}

void GridCell::incCoreNum()
{
	coreNum++;
}

int GridCell::getCoreNum()
{
	return this->coreNum;
}

int GridCell::getDiffPos()
{
	return this->diffPos;
}