#pragma once
#include <cstdio>
#include <vector>
#include "DataPoint.h"

class TreeNodeVec {
protected:
	std::vector<int> prefix;
	int offset;
	int next;
	int header;
	int size;
	int key;
	int gridID;
public:
	TreeNodeVec() {}
	TreeNodeVec(int h, int key, int g) :
		offset(0), next(-1), header(h), size(1), key(key), gridID(g) { }

	void print();
	std::vector<int>& getPrefix();
	int getOffset();
	void setOffset(int o);
	int getNext();
	void setNext(int n);
	int getHeader();
	int getSize();
	void incSize();
	int getKey();
	int getGridID();
};