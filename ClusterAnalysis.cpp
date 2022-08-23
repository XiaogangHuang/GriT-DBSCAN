#pragma warning(disable:4996)
#include "ClusterAnalysis.h"

// read the data set from fileName
void ClusterAnalysis::Init(const char* fileName, float radius, int minPTs) {
    this->radius = radius;        //set the radius
    this->minPTs = minPTs;        //set the minPTs
    this->clusterId = 0;
    int dim;
    int data_size;
    printf("loading the file...\n");
    printf("reading data ...\n");
    raw_data = read_data(fileName, " ", &dim, &data_size);
    dataNum = data_size;
    printf("dimension is: %d\n", dim - 1);
    printf("dataNum is: %d\n", dataNum);
}

// run GriT-DBSCAN
void ClusterAnalysis::DoDBSCANRecursive() {
    // construct grid structure
    double start = getCurrentTime();
    constructGridCell();
    printf("the constructGridCell time is %f\n", getCurrentTime() - start);
    printf("total number of cells: %zd\n", gridCells.size());

    // construct grid tree
    start = getCurrentTime();
    ConstructGridTree();
    printf("the ConstructGridTree time is %f\n", getCurrentTime() - start);
    // find the non-empty neighboring grids for each grid
    start = getCurrentTime();
    FindNeighborGridHalf();
    /*if (gridCells.size() > 9e5)
        FindNeighborGridHalf();
    else
        FindNeighborGrid();*/
    printf("the FindNeighborGrid time is %f\n", getCurrentTime() - start);
    printf("average neighbor number: %d\n", myAverage(neiGrid));
    //step 1: identifyCorePts
    start = getCurrentTime();
    identifyCorePts();
    printf("the identifyCorePts time is %f\n", getCurrentTime() - start);
    //step 2: Creating a grid graph
    start = getCurrentTime();
    constructGridGraph();
    printf("the constructGridGraph time is %f\n", getCurrentTime() - start);
    //step 3: assignNonCorePts
    start = getCurrentTime();
    assignNonCorePts();
    printf("the assignNonCorePts time is %f\n", getCurrentTime() - start);
}

// find the minimum value in each dimension
void ClusterAnalysis::getDataMin(float* result)
{
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        result[i] = raw_data[i];
    }
    for (int i = 1; i < dataNum; i++) {
        float* temp = raw_data + i * DIME_NUM;
        for (int j = 0; j < DIME_NUM - 1; j++)
        {
            if (temp[j] < result[j])
                result[j] = temp[j];
        }
    }
}

// Sorting the grid based on the grid's coordinates
void ClusterAnalysis::RaidxSort(Matrix& gridCoords)
{
    Matrix gridTemp;
    gridTemp.alloc(dataNum, DIME_NUM);
    int mx = *max_element(gridCoords[0], gridCoords[dataNum - 1] + DIME_NUM - 1) + 1;
    for (int i = 0; i < dataNum; i++)
    {
        gridCoords[i][DIME_NUM - 1] = i;
    }
    int* count = (int*)malloc(mx * sizeof(int));
    for (int j = DIME_NUM - 2; j >= 0; j--)
    {
        memset(count, 0, mx * sizeof(int));
        for (int i = 0; i < dataNum; i++)
        {
            count[gridCoords[i][j]] += 1;
        }
        for (int i = 1; i < mx; i++)
            count[i] += count[i - 1];
        for (int i = dataNum - 1; i >= 0; i--)
        {
            int pos = --count[gridCoords[i][j]];
            for (int d = 0; d < DIME_NUM; d++)
            {
                gridTemp[pos][d] = gridCoords[i][d];
            }
        }
        gridCoords.my_swap(gridTemp);
    }
    free(count);
}

// function for constructing grid structure
void ClusterAnalysis::constructGridCell()
{
    // calculate the grid coordinates
    Matrix gridCoords;
    gridCoords.alloc(dataNum, DIME_NUM);
    float dataMin[DIME_NUM - 1];
    getDataMin(dataMin);
    for (int i = 0; i < dataNum; i++)
    {
        float* temp = raw_data + i * DIME_NUM;
        for (int j = 0; j < DIME_NUM - 1; j++)
        {
            gridCoords[i][j] = (int)((temp[j] - dataMin[j]) * sqrt(DIME_NUM - 1) / radius);
        }
    }
    for (int j = 0; j < DIME_NUM - 1; j++)
        coordMax[j] = gridCoords[0][j] + 1;
    for (int i = 1; i < dataNum; i++)
    {
        for (int j = 0; j < DIME_NUM - 1; j++)
        {
            if (gridCoords[i][j] >= coordMax[j])
            {
                coordMax[j] = gridCoords[i][j] + 1;
            }
        }
    }
    RaidxSort(gridCoords);
    for (int i = 0; i < dataNum; i++) {
        DataPoint tempDP(raw_data + (gridCoords[i][DIME_NUM - 1]) * (DIME_NUM), 
            gridCoords[i][DIME_NUM - 1]);
        dadaSets.push_back(tempDP);
    }
    free(raw_data);

    // construct the grid structure
    vector<int> gridcell;
    gridcell.push_back(0);
    int TNCounts[DIME_NUM - 1];
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        TNCounts[i] = 1;
    }
    for (int i = 1; i < dataNum; i++)
    {
        int cond = 0;
        for (int j = 0; j < DIME_NUM - 1; j++)
        {
            if (cond == 0)
            {
                if (gridCoords[i][j] != gridCoords[i - 1][j])
                {
                    GridCell temp(gridCoords[i - 1], gridcell, gridcell.size(), j);
                    gridCells.push_back(temp);
                    gridcell.clear();
                    gridcell.push_back(i);
                    TNCounts[j]++;
                    cond = 1;
                }
            }
            else
            {
                TNCounts[j]++;
            }
        }
        if (gridcell[gridcell.size() - 1] != i)
        {
            gridcell.push_back(i);
        }
    }
    GridCell temp(gridCoords[dataNum - 1], gridcell, gridcell.size(), DIME_NUM);
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        vector<TreeNodeVec> TreeVec;
        TreeVec.reserve(TNCounts[i]);
        gridTreeVec.push_back(TreeVec);
    }
    gridCells.push_back(temp);
}

// function for constructing grid tree
void ClusterAnalysis::ConstructGridTree()
{
    int* coord = gridCells[0].getCoords();
    int pre[DIME_NUM];
    memset(pre, 0, DIME_NUM * sizeof(int));
    for (int i = 0; i < DIME_NUM - 1; i++)
    {
        TreeNodeVec node(0, coord[i], 0);
        gridTreeVec[i].push_back(node);
    }
    for (int i = 1; i < gridCells.size(); i++)
    {
        coord = gridCells[i].getCoords();
        int j = gridCells[i - 1].getDiffPos();
        TreeNodeVec node(pre[j + 1] + 1, coord[j], i);
        gridTreeVec[j].push_back(node);
        gridTreeVec[j][pre[j]].setNext(pre[j] + 1);
        if (j > 0)
        {
            gridTreeVec[j - 1][pre[j - 1]].incSize();
        }
        pre[j]++;
        for (j++; j < DIME_NUM - 1; j++)
        {
            TreeNodeVec node(pre[j + 1] + 1, coord[j], i);
            gridTreeVec[j].push_back(node);
            pre[j]++;
        }
    }
    for (int i = 0; i < DIME_NUM - 2; i++)
    {
        for (int j = 0; j < gridTreeVec[i].size(); j++)
        {
            if (gridTreeVec[i][j].getSize() > minPTs)
            {
                vector<int>& prefix = gridTreeVec[i][j].getPrefix();
                vector<int> pre(coordMax[i + 1], -1);
                int cur = gridTreeVec[i][j].getHeader();
                while (cur != -1)
                {
                    int temp = gridTreeVec[i + 1][cur].getKey();
                    pre[temp] = cur;
                    cur = gridTreeVec[i + 1][cur].getNext();
                }
                int ps = 0;
                int maxDiff = ceil(sqrt(DIME_NUM));
                for (int t = 0; t < coordMax[i + 1]; t++)
                {
                    if (ps < t - maxDiff)
                        ps = t - maxDiff;
                    while (ps < coordMax[i + 1] - 1 && pre[ps] == -1 && ps <= t + maxDiff)
                        ps++;
                    prefix.push_back(pre[ps]);
                }
            }
        }
    }
}

// function for finding non-empty neighboring grids
void ClusterAnalysis::FindNeighborGrid()
{
    Matrix cacheQueue;
    cacheQueue.alloc(DIME_NUM - 1, gridCells.size() + 1);
    neiGrid.reserve(gridCells.size());
    int pos = 0;
    int maxDiff = ceil(sqrt(DIME_NUM));
    for (int i = 0; i < gridCells.size(); i++)
    {
        int* coord = gridCells[i].getCoords();
        if (pos == 0)
        {
            int cur = cacheQueue[0][1];
            while (coord[0] - gridTreeVec[0][cur].getKey() > maxDiff && cur != -1)
            {
                cur = gridTreeVec[0][cur].getNext();
            }
            int j = 1;
            while (cur != -1 && gridTreeVec[0][cur].getKey() - coord[0] <= maxDiff)
            {
                cacheQueue[0][j++] = cur;
                int diff = gridTreeVec[0][cur].getKey() - coord[0];
                if (diff < -1)
                    gridTreeVec[0][cur].setOffset((diff + 1) * (diff + 1));
                else if (diff > 1)
                    gridTreeVec[0][cur].setOffset((diff - 1) * (diff - 1));
                else
                    gridTreeVec[0][cur].setOffset(0);
                cur = gridTreeVec[0][cur].getNext();
            }
            cacheQueue[0][0] = j - 1;
            pos = 1;
        }
        for (int j = pos; j < DIME_NUM - 1; j++)
        {
            int key = coord[j];
            int sz = cacheQueue[j - 1][0];
            int szj = 1;
            for (int t = 1; t <= sz; t++)
            {
                int que = -1;
                TreeNodeVec& treeNodeTemp =
                    gridTreeVec[j - 1][cacheQueue[j - 1][t]];
                int offset = treeNodeTemp.getOffset();
                if (treeNodeTemp.getSize() > minPTs)
                    que = treeNodeTemp.getPrefix()[key];
                else
                {
                    que = treeNodeTemp.getHeader();
                    if (key < gridTreeVec[j][que].getKey() - maxDiff ||
                        gridTreeVec[j][que + treeNodeTemp.getSize() - 1].getKey() < 
                        key - maxDiff)
                        continue;
                    while (que != -1 && gridTreeVec[j][que].getKey() < key - maxDiff)
                    {
                        que = gridTreeVec[j][que].getNext();
                    }
                }
                int diff = gridTreeVec[j][que].getKey() - key;
                while (que != -1 && diff <= maxDiff)
                {
                    int offsetTemp = offset;
                    if (diff < -1)
                        offsetTemp += (diff + 1) * (diff + 1);
                    else if (diff > 1)
                        offsetTemp += (diff - 1) * (diff - 1);
                    if (offsetTemp < DIME_NUM - 1)
                    {
                        gridTreeVec[j][que].setOffset(offsetTemp);
                        cacheQueue[j][szj++] = que;
                    }
                    que = gridTreeVec[j][que].getNext();
                    diff = gridTreeVec[j][que].getKey() - key;
                }
            }
            cacheQueue[j][0] = szj - 1;
        }
        int* queue = cacheQueue[DIME_NUM - 2];
        vector<int> nei_sort;
        nei_sort.reserve(queue[0]);
        for (int j = 1; j <= queue[0]; j++)
        {
            nei_sort.push_back(move(queue[j]));
        }
        RaidxSort_Nei(nei_sort);
        for (int j = 0; j < nei_sort.size(); j++)
        {
            nei_sort[j] = gridTreeVec[DIME_NUM - 2][nei_sort[j]].getGridID();
        }
        neiGrid.push_back(move(nei_sort));
        nei_sort.clear();
        pos = gridCells[i].getDiffPos();
    }
}

// function for finding non-empty neighboring grids by search half of the space
void ClusterAnalysis::FindNeighborGridHalf()
{
    Matrix cacheQueue;
    cacheQueue.alloc(DIME_NUM - 1, gridCells.size() + 1);
    neiGrid.reserve(gridCells.size());
    int pos = 0;
    int maxDiff = ceil(sqrt(DIME_NUM));
    for (int i = 0; i < gridCells.size(); i++)
    {
        int* coord = gridCells[i].getCoords();
        if (pos == 0)
        {
            int cur = cacheQueue[0][1];
            while (coord[0] - gridTreeVec[0][cur].getKey() > maxDiff && cur != -1)
            {
                cur = gridTreeVec[0][cur].getNext();
            }
            int j = 1;
            while (cur != -1 && gridTreeVec[0][cur].getKey() <= coord[0])
            {
                cacheQueue[0][j++] = cur;
                int diff = gridTreeVec[0][cur].getKey() - coord[0];
                if (diff < -1)
                    gridTreeVec[0][cur].setOffset((diff + 1) * (diff + 1));
                else if (diff > 1)
                    gridTreeVec[0][cur].setOffset((diff - 1) * (diff - 1));
                else
                    gridTreeVec[0][cur].setOffset(0);
                cur = gridTreeVec[0][cur].getNext();
            }
            cacheQueue[0][0] = j - 1;
            pos = 1;
        }
        for (int j = pos; j < DIME_NUM - 1; j++)
        {
            int key = coord[j];
            int sz = cacheQueue[j - 1][0];
            int szj = 1;
            for (int t = 1; t <= sz; t++)
            {
                int que = -1;
                TreeNodeVec& treeNodeTemp =
                    gridTreeVec[j - 1][cacheQueue[j - 1][t]];
                int offset = treeNodeTemp.getOffset();
                if (treeNodeTemp.getSize() > minPTs)
                    que = treeNodeTemp.getPrefix()[key];
                else
                {
                    que = treeNodeTemp.getHeader();
                    if (key < gridTreeVec[j][que].getKey() - maxDiff ||
                        gridTreeVec[j][que + treeNodeTemp.getSize() - 1].getKey() <
                        key - maxDiff)
                        continue;
                    while (que != -1 && gridTreeVec[j][que].getKey() < key - maxDiff)
                    {
                        que = gridTreeVec[j][que].getNext();
                    }
                }
                int diff = gridTreeVec[j][que].getKey() - key;
                while (que != -1 && diff <= maxDiff && gridTreeVec[j][que].getGridID() <= i)
                {
                    int offsetTemp = offset;
                    if (diff < -1)
                        offsetTemp += (diff + 1) * (diff + 1);
                    else if (diff > 1)
                        offsetTemp += (diff - 1) * (diff - 1);
                    if (offsetTemp < DIME_NUM - 1)
                    {
                        gridTreeVec[j][que].setOffset(offsetTemp);
                        cacheQueue[j][szj++] = que;
                    }
                    que = gridTreeVec[j][que].getNext();
                    diff = gridTreeVec[j][que].getKey() - key;
                }
            }
            cacheQueue[j][0] = szj - 1;
        }
        int* queue = cacheQueue[DIME_NUM - 2];
        vector<int> nei_sort;
        nei_sort.reserve(2 * queue[0]);
        nei_sort.push_back(i);
        for (int j = 1; j <= queue[0] - 1; j++)
        {
            int gridId = gridTreeVec[DIME_NUM - 2][queue[j]].getGridID();
            nei_sort.push_back(move(queue[j]));
            neiGrid[gridId].push_back(i);
        }
        RaidxSort_Nei_Half(nei_sort);
        neiGrid.push_back(move(nei_sort));
        nei_sort.clear();
        pos = gridCells[i].getDiffPos();
    }
}

// sort the non-empty neighboring grids
void ClusterAnalysis::RaidxSort_Nei(vector<int>& vec)
{
    int* count = (int*)calloc(DIME_NUM, sizeof(int));
    for (int i = 0; i < vec.size(); i++)
    {
        count[gridTreeVec[DIME_NUM - 2][vec[i]].getOffset()] += 1;
    }
    for (int i = 1; i < DIME_NUM; i++)
        count[i] += count[i - 1];
    vector<int> vec_new(vec);
    for (int i = vec_new.size() - 1; i >= 0; i--)
    {
        int pos = --count[gridTreeVec[DIME_NUM - 2][vec_new[i]].getOffset()];
        vec[pos] = vec_new[i];
    }
    free(count);
}

void ClusterAnalysis::RaidxSort_Nei_Half(vector<int>& vec)
{
    int* count = (int*)calloc(DIME_NUM, sizeof(int));
    for (int i = 0; i < vec.size(); i++)
    {
        count[gridTreeVec[DIME_NUM - 2][vec[i]].getOffset()] += 1;
    }
    for (int i = 1; i < DIME_NUM; i++)
        count[i] += count[i - 1];
    vector<int> vec_new(vec);
    for (int i = vec_new.size() - 1; i >= 0; i--)
    {
        int pos = --count[gridTreeVec[DIME_NUM - 2][vec_new[i]].getOffset()];
        vec[pos] = gridTreeVec[DIME_NUM - 2][vec_new[i]].getGridID();
    }
    free(count);
}

// identify core points in the data set
void ClusterAnalysis::identifyCorePts()
{
    float sqr_r = radius * radius;
    for (int i = 0; i < gridCells.size(); i++)
    {
        int sz = gridCells[i].getSize();
        vector<int>& pts = gridCells[i].getPoints();
        if (sz >= minPTs)
        {
            gridCells[i].setCoreNum(sz);
            for (int query = 0; query < sz; query++)
                dadaSets[pts[query]].SetCore();
        }
        else
        {
            for (int query = 0; query < sz; query++)
            {
                int cnt = 0;
                float* point = dadaSets[pts[query]].GetDimension();
                for (int neiCell = 0; neiCell < neiGrid[i].size(); neiCell++)
                {
                    vector<int>& neiPts = gridCells[neiGrid[i][neiCell]].getPoints();
                    for (int j = 0; j < neiPts.size(); j++)
                    {
                        float* neipoint = dadaSets[neiPts[j]].GetDimension();
                        if (sqrDistance_1(point, neipoint) <= sqr_r) {
                            cnt++;
                            if (cnt >= minPTs)
                                break;
                        }
                    }
                    if (cnt >= minPTs)
                    {
                        dadaSets[pts[query]].SetCore();
                        gridCells[i].incCoreNum();
                        break;
                    }
                }
            }
        }
    }
    int count = 0;
    for (int i = 0; i < gridCells.size(); i++)
    {
        count += gridCells[i].getCoreNum();
    }
    printf("total number of core points : %d\n", count);
}

// construct the graph
void ClusterAnalysis::constructGridGraph()
{
    vector<int> clu(gridCells.size(), -1);
    for (int i = 0; i < gridCells.size(); i++)
    {
        if (gridCells[i].getCoreNum() == 0)
            clu[i] = -2;
    }
    for (int i = 0; i < gridCells.size(); i++)
    {
        if (clu[i] == -1)
        {
            clu[i] = clusterId;
            vector<int> temp;
            temp.push_back(i);
            int pos = 0;
            while (pos < temp.size())
            {
                int current = temp[pos++];
                for (int j = 0; j < neiGrid[current].size(); j++)
                {
                    int neigrid = neiGrid[current][j];
                    if (clu[neigrid] == -1 && Merge_fast(current, neigrid) == true)
                    {
                        clu[neigrid] = clusterId;
                        temp.push_back(neigrid);
                    }
                }
            }
            clusterId += 1;
        }
    }
    for (int jt = 0; jt < gridCells.size(); jt++)
    {
        int clust = clu[jt];
        if (clust >= 0)
        {
            vector<int>& pts = gridCells[jt].getPoints();
            for (int i = 0; i < pts.size(); i++)
            {
                dadaSets[pts[i]].SetClusterId(clust);
            }
        }
    }
}

// sort the grids based on the number of core points
void ClusterAnalysis::RaidxSort_Grid(vector<int>& vec)
{
    int mx = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        if (gridCells[vec[i]].getCoreNum() >= mx)
            mx = gridCells[vec[i]].getCoreNum() + 1;
    }
    int* count = (int*)calloc(mx, sizeof(int));
    for (int i = 0; i < vec.size(); i++)
    {
        count[gridCells[vec[i]].getCoreNum()] += 1;
    }
    for (int i = 1; i < mx; i++)
        count[i] += count[i - 1];
    vector<int> vec_new(vec);
    for (int i = vec_new.size() - 1; i >= 0; i--)
    {
        int pos = --count[gridCells[vec_new[i]].getCoreNum()];
        vec[pos] = vec_new[i];
    }
    free(count);
}

// fast merging algorithm
bool ClusterAnalysis::Merge_fast(int a, int b) {
    int n1 = gridCells[a].getSize();
    int n2 = gridCells[b].getSize();
    if (gridCells[a].getCoreNum() < 20 || gridCells[a].getCoreNum() < 20)
    {
        return Merge(a, b);
    }
    int* visited_a = (int*)calloc(n1 + 1, sizeof(int));
    visited_a[n1] = gridCells[a].getCoreNum();
    vector<int>& pts_a = gridCells[a].getPoints();
    int query;
    for (int i = 0; i < n1; i++)
    {
        if (dadaSets[pts_a[i]].isCore() == 0)
        {
            visited_a[i] = 1;
        }
        else
        {
            query = i;
        }
    }
    int* visited_b = (int*)calloc(n2 + 1, sizeof(int));
    visited_b[n2] = gridCells[b].getCoreNum();
    vector<int>& pts_b = gridCells[b].getPoints();
    for (int i = 0; i < n2; i++)
    {
        if (dadaSets[pts_b[i]].isCore() == 0)
        {
            visited_b[i] = 1;
        }
    }
    while (visited_a[n1] * visited_b[n2] != 0)
    {
        if (Filter(a, b, visited_a, visited_b, query) ||
            Filter(b, a, visited_b, visited_a, query))
        {
            free(visited_a);
            free(visited_b);
            return true;
        }
    }
    free(visited_a);
    free(visited_b);
    return false;
}

// merging algorithm using brute force
bool ClusterAnalysis::Merge(int a, int b) {
    int n1 = gridCells[a].getSize();
    int n2 = gridCells[b].getSize();
    float sqr_r = radius * radius;
    vector<int>& pts = gridCells[a].getPoints();
    for (int i = 0; i < n1; i++)
    {
        if (dadaSets[pts[i]].isCore())
        {
            float* pt = dadaSets[pts[i]].GetDimension();
            vector<int>& pts_b = gridCells[b].getPoints();
            for (int j = 0; j < n2; j++)
            {
                if (dadaSets[pts_b[j]].isCore() &&
                    sqrDistance_1(pt, dadaSets[pts_b[j]].GetDimension()) <= sqr_r)
                {
                    return true;
                }
            }
        }
    }
    return false;
}

// prune trivial points in grid a
bool ClusterAnalysis::Filter(int a, int b, int* visited_a, int* visited_b, int& query)
{
    int q = 0;
    int n1 = gridCells[a].getSize();
    int n2 = gridCells[b].getSize();
    vector<int>& coreGrid_a = gridCells[a].getPoints();
    vector<int>& coreGrid_b = gridCells[b].getPoints();
    float* distTemp = (float*)calloc(n2, sizeof(float));
    float minidist = FLT_MX;
    visited_a[query] = 1;
    visited_a[n1] -= 1;
    float* queryPoint = dadaSets[coreGrid_a[query]].GetDimension();
    for (int i = 0; i < n2; i++) {
        if (visited_b[i] == 0)
        {
            float* point = dadaSets[coreGrid_b[i]].GetDimension();
            distTemp[i] = distance_1(queryPoint, point);
            if (distTemp[i] <= radius)
                return true;
            else if (distTemp[i] < minidist)
            {
                minidist = distTemp[i];
                q = i;
            }
        }
    }
    float pt[DIME_NUM - 1];
    float temp = 0;
    for (int i = 0; i < DIME_NUM - 1; i++) {
        float a = queryPoint[i];
        pt[i] = a - dadaSets[coreGrid_b[q]].GetDimension()[i];
        temp += a * pt[i];
    }
    float theta = 1;
    for (int i = 0; i < n2; i++) {
        if (visited_b[i] == 0)
        {
            float alpha = temp;
            for (int j = 0; j < DIME_NUM - 1; j++) {
                alpha -= dadaSets[coreGrid_b[i]].GetDimension()[j] * pt[j];
            }
            alpha = alpha / (distTemp[i] * minidist);
            float beta = radius / distTemp[i];
            float arc = alpha * pow(1 - beta * beta, 0.5F);
            if (1 - alpha * alpha > 0)
                arc -= beta * pow(1 - alpha * alpha, 0.5F);
            if (arc < theta)
            {
                theta = arc;
            }
        }
    }
    for (int i = 0; i < n1; i++) {
        if (visited_a[i] == 0)
        {
            float* point = dadaSets[coreGrid_a[i]].GetDimension();
            float dist = distance_1(queryPoint, point);
            float alpha = temp;
            for (int j = 0; j < DIME_NUM - 1; j++) {
                alpha -= point[j] * pt[j];
            }
            if (alpha < theta * dist * minidist || dist < minidist - radius)
            {
                visited_a[i] = 1;
                --visited_a[n1];
            }
        }
    }
    query = q;
    free(distTemp);
    return false;
}

// assign each non core point to corresponding cluster or noise
void ClusterAnalysis::assignNonCorePts()
{
    float sqr_r = radius * radius;
    for (int i = 0; i < gridCells.size(); i++)
    {
        int coreNum = gridCells[i].getCoreNum();
        int sz = gridCells[i].getSize();
        if (coreNum == 0)
        {
            vector<int>& pts = gridCells[i].getPoints();
            for (int query = 0; query < pts.size(); query++)
            {
                int cond = 0;
                float* point = dadaSets[pts[query]].GetDimension();
                for (int neiCell = 0; neiCell < neiGrid[i].size(); neiCell++)
                {
                    if (gridCells[neiGrid[i][neiCell]].getCoreNum() > 0)
                    {
                        vector<int>& neiPts = gridCells[neiGrid[i][neiCell]].getPoints();
                        for (int j = 0; j < neiPts.size(); j++)
                        {
                            if (dadaSets[neiPts[j]].isCore())
                            {
                                float* neipoint = dadaSets[neiPts[j]].GetDimension();
                                if (sqrDistance_1(point, neipoint) <= sqr_r) {
                                    dadaSets[pts[query]].SetClusterId(dadaSets[neiPts[j]].GetClusterId());
                                    break;
                                }
                            }
                        }
                        if (dadaSets[pts[query]].GetClusterId() != -1)
                            break;
                    }
                }
            }
        }
    }
}

// print the clustering results
void ClusterAnalysis::printMessage()
{
    printf("There are %d clusters\n", clusterId);
    int* temp = (int*)calloc(clusterId, sizeof(int));
    int noisePts = 0;
    for (int i = 0; i < dataNum; i++)
    {
        if (dadaSets[i].GetClusterId() != -1)
        {
            temp[dadaSets[i].GetClusterId()] += 1;
        }
        else
            noisePts++;
    }
    for (int i = 0; i < clusterId; i++)
    {
        printf("Cluster %d: %d pts\n", i + 1, temp[i]);
    }
    printf("There are %d noise points\n", noisePts);
    free(temp);
}


/********************************/
/*    functions for file I/O    */
/********************************/

int ClusterAnalysis::get_dim(char* s, const char* delims) {
    char* val_str = NULL;
    val_str = strtok(s, delims);
    int dim = 0;
    while (val_str != NULL) {
        dim++;
        val_str = strtok(NULL, delims);
    }
    return dim;
}

float* ClusterAnalysis::get_data(char* s, int dim, const char* delims) {
    float* temp = (float*)malloc(dim * sizeof(float));
    char* val_str = NULL;
    val_str = strtok(s, delims);
    int counter = 0;
    while (val_str != NULL) {
        temp[counter] = atof(val_str);
        counter++;
        val_str = strtok(NULL, delims);
    }
    return temp;
}

void ClusterAnalysis::read_data_dim_size(const char* filename, int* data_dim, int* data_size, const char* delims) {
    int n_size = 0;
    int dim = 0;
    char s[10000];
    freopen(filename, "r", stdin);
    while (gets_s(s))
    {
        if (dim == 0)
            dim = get_dim(s, delims);
        n_size++;
    }
    *data_dim = dim;
    *data_size = n_size;
    fclose(stdin);
}

float* ClusterAnalysis::read_data(const char* filename, const char* delims) {
    int dim, n_size;
    read_data_dim_size(filename, &dim, &n_size, delims);

    float* data = (float*)malloc(n_size * dim * sizeof(float));
    freopen(filename, "r", stdin);
    int counter = 0;
    char s[10000];
    while (gets_s(s))
    {
        float* tmp_data = get_data(s, dim, delims);
        memcpy(data + counter * dim, tmp_data, dim * sizeof(float));
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

float* ClusterAnalysis::read_data(const char* filename, const char* delims, int* dim, int* data_size) {
    read_data_dim_size(filename, dim, data_size, delims);

    float* data = (float*)malloc((*data_size) * (*dim) * sizeof(float));
    freopen(filename, "r", stdin);
    int counter = 0;
    char s[10000];
    while (gets_s(s))
    {
        float* tmp_data = get_data(s, *dim, delims);
        memcpy(data + counter * (*dim), tmp_data, (*dim) * sizeof(float));
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

// write the results to fileName
bool ClusterAnalysis::WriteToFile(const char* fileName) {
    ofstream of1(fileName);
    for (int i = 0; i < dataNum; i++)
    {
        of1 << dadaSets[i].GetOldIndex() << " " <<
            dadaSets[i].GetClusterId() << " " << dadaSets[i].isCore() << endl;
    }
    of1.close();
    return true;
}
