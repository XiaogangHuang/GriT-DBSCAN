#pragma once
#include <cstdio>

class Matrix {
public:
	// constructor with number of rows and columns
	Matrix(): rows(0), columns(0),p(NULL) {

	}
	void alloc(int rows, int columns){
		this->rows = rows;
		this->columns = columns;
		// allocate memory
		p = (int*)malloc(rows * columns * sizeof(int));
		if (p == 0) {
			printf("No enough memory in alloc matrix!\n");
		}
		else {
			// initialize all elements to zero
			memset(p, 0, rows * columns * sizeof(int));
		}
	}
	void realloc(int rows, int columns) {
		this->rows = rows;
		this->columns = columns;
		memset(p, 0, rows * columns * sizeof(int));
	}
	// get number of rows
	int nrows() const {
		return rows;
	}
	// get number of columns
	int ncols() const {
		return columns;
	}
	// destructor
	~Matrix() {
		if (p)
		{
			free(p);
		}
	}
	// Operator [] gives a pointer to row number r.
	// A second [] with the column number gives a single element
	int* operator[] (int r) {
		// no bounds check
		return p + r * columns;
	}
	void my_swap(Matrix& a)
	{
		int* t = p;
		rows = a.rows;
		columns = a.columns;
		p = a.p;
		a.p = t;
	}

protected:
	int rows; // number of rows
	int columns; // number of columns
	int* p; // pointer to memory
};