#include "TreeNode.h"

void TreeNode::print()
{
	printf(" next: %d", next); 
	printf(" tail: %d", tail);
	printf(" key: %d", key);
	printf(" gridID: %d\n", gridID);
}

int TreeNode::getOffset()
{
	return this->offset;
}

void TreeNode::setOffset(int o)
{
	this->offset = o;
}

int TreeNode::getNext()
{
	return this->next;
}

void TreeNode::setNext(int n)
{
	this->next = n;
}

int TreeNode::getTail()
{
	return this->tail;
}

void TreeNode::setTail(int t)
{
	this->tail = t;
}

int TreeNode::getKey()
{
	return this->key;
}

int TreeNode::getGridID()
{
	return this->gridID;
}
