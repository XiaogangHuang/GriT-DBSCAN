#pragma once
#include <cstdio>
#include <vector>
#include "DataPoint.h"

class TreeNode {
protected:
	int offset;
	int next;
	int tail;
	int key;
	int gridID;
public:
	TreeNode() {}
	TreeNode(int t, int key, int g) :
		offset(0), next(-1), tail(t), key(key), gridID(g) {}
	
	void print();
	int getOffset();
	void setOffset(int o);
	int getNext();
	void setNext(int n);
	int getTail();
	void setTail(int t);
	int getKey();
	int getGridID();
};