#include <algorithm>
#include <windows.h>
#include <cstdlib>
#include <fstream>
#include <iosfwd>
#include "TreeNodeVec.h"
#include "GridCell.h"
#include "Matrix.h"
using namespace std;

#define FLT_MX 3.402823466e+38F

//Cluster analysis type.
class ClusterAnalysis
{
private:
    vector< DataPoint > dadaSets;
    vector< GridCell > gridCells;
    vector< vector<TreeNodeVec> > gridTreeVec;
    vector< vector<int> > neiGrid;
    float* raw_data;
    float radius;
    int coordMax[DIME_NUM - 1];
    int dataNum;
    int minPTs;
    int clusterId;
public:
    ClusterAnalysis() {}                //Default constructor.
    void Init(const char* fileName, float radius, int minPTs);    //Initialization operation

    void DoDBSCANRecursive();
    void constructGridCell();
    void getDataMin(float* result);
    void RaidxSort(Matrix& gridCoords);
    void ConstructGridTree();
    void FindNeighborGrid();
    void FindNeighborGridHalf();
    void RaidxSort_Nei(vector<int>& vec);
    void RaidxSort_Nei_Half(vector<int>& vec);
    void identifyCorePts();
    void constructGridGraph();
    void RaidxSort_Grid(vector<int>& vec);
    bool Merge_fast(int a, int b);
    bool Merge(int a, int b);
    bool Filter(int a, int b, int* visited_a, int* visited_b, int& query);
    void assignNonCorePts();
    void printMessage();

    int get_dim(char* s, const char* delims);
    float* get_data(char* s, int dim, const char* delims);
    void read_data_dim_size(const char* filename, int* data_dim, int* data_size, const char* delims);
    float* read_data(const char* filename, const char* delims);
    float* read_data(const char* filename, const char* delims, int* dim, int* data_size);
    bool WriteToFile(const char* fileName);    //save results
};
