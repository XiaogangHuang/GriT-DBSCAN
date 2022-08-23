#include "Toolfunction.h"

double getCurrentTime() {
    long long time = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now().time_since_epoch()).count();
    return time / 1000000.0;
}

float distance_1(float* p1, float* p2)
{
    //Euclidean distance
    float sum = 0.0;
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        float d1 = p1[i] - p2[i];
        d1 *= d1;
        sum += d1;
    }
    return sqrt(sum);

    //Chebyshev distance
    //    float sum =0;
    //    for(int i = 0; i < DIME_NUM-1;i++)
    //    {
    //        float d2 =fabsf(p2[i] -p1[i]);
    //        if(d2 > sum)
    //            sum = d2;
    //    }
    //    return sum;
}

float sqrDistance_1(float* p1, float* p2)
{
    float sum = 0.0;
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        float d1 = p1[i] - p2[i];
        d1 *= d1;
        sum += d1;
    }
    return sum;
}

int myAverage(vector< vector<int> >& vec)
{
    long long sum = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        sum += vec[i].size();
    }
    return sum/vec.size();
}

int myMin(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}
