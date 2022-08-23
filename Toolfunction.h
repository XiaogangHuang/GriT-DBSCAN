#pragma once
#include <chrono>
#include <vector>
#include "TreeNode.h"

double getCurrentTime();  // return the current time

float distance_1(float* p1, float* p2); // calculate Euclidean distance

float sqrDistance_1(float* p1, float* p2); // calculate the square of the Euclidean distance

int myAverage(vector< vector<int> >& vec);	// calculate the mean of a vector

int myMin(int a, int b);	// return the smaller of two numbers
