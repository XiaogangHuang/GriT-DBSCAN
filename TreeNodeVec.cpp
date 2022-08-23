#include "TreeNodeVec.h"

/** Print the tree node information **/
void TreeNodeVec::print()
{
	printf("next: %d header: %d size: %d key: %d gridID: %d\n", next,
		header, size, key, gridID);
	if (size > 10)
	{
		for (int i = 0; i < prefix.size(); i++)
		{
			printf("%d ", prefix[i]);
		}
		printf("\n");
	}
}

/** return the vector of prefix **/
std::vector<int>& TreeNodeVec::getPrefix()
{
	return this->prefix;
}


/** return current offset **/
int TreeNodeVec::getOffset()
{
	return this->offset;
}


/** set current offset **/
void TreeNodeVec::setOffset(int o)
{
	this->offset = o;
}


/** return next node **/
int TreeNodeVec::getNext()
{
	return this->next;
}


/** set next node **/
void TreeNodeVec::setNext(int n)
{
	this->next = n;
}


/** return header of next level **/
int TreeNodeVec::getHeader()
{
	return this->header;
}


/** return the size **/
int TreeNodeVec::getSize()
{
	return this->size;
}

/** increment the size **/
void TreeNodeVec::incSize()
{
	this->size++;
}

/** return the key **/
int TreeNodeVec::getKey()
{
	return this->key;
}

/** return the corresponding grid id **/
int TreeNodeVec::getGridID()
{
	return this->gridID;
}
