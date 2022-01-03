#include "../include/algorithms.hpp"

using namespace std;


// Convert ClassCurve to Fred's curve to calculate continuous frechet distance between curves.
Curve cont_convert_curve(ClassCurve curve)
{
    Points FredPoints(curve.cpoints[0].vpoint.size());
    for (int i = 0; i < curve.cpoints.size(); i++)
    {
        Point FredPoint(1);
        FredPoint.assign(1,curve.cpoints[i].vpoint[0]);
        //std::cout << FredPoint << endl;
        FredPoints.push_back(FredPoint);
    }
    Curve FredCurve(FredPoints,"0");
    return FredCurve;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TRUE //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pair<ClassPoint,double> true_NN(ClassPoint q, Vector_of_points inputData, double &time)
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    ClassPoint b;
    pair<ClassPoint,double> best;
    best.first = b; //best point-candidate
    best.second = DBL_MAX; //best distance of best candidate
    
    for (int i = 0; i < inputData.points.size(); i++)
    {
        double dist = distance(q.vpoint,inputData.points[i].vpoint, 2);
        if (dist < best.second)
        {
            best.second = dist;
            best.first = inputData.points[i];
        }
        
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();
    return best;

}


set<pair<ClassPoint,double>, CompDist> true_nNN(ClassPoint q, int N, Vector_of_points inputData, double &time)
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    // Initialise a set two hold pairs of true best point/best distance.
    set<pair<ClassPoint,double>, CompDist> bestPointsDists;
    ClassPoint a;
    bestPointsDists.insert(make_pair(a,DBL_MAX));
    
    for (int i = 0; i < inputData.points.size(); i++)
    {
        double dist = distance(q.vpoint,inputData.points[i].vpoint, 2);
        if (bestPointsDists.size()==N) //if set is full
        {
            //if the biggest distance in set is equal/greater than current distance
            if (prev(bestPointsDists.end())->second >= dist)
            {
                bestPointsDists.erase(prev(bestPointsDists.end())); //pop biggest distance pair
                bestPointsDists.insert(make_pair(inputData.points[i],dist)); //insert new point/distance
            }
        }
        else if (bestPointsDists.size()<N) //if there is space in set insert pair
        {
            bestPointsDists.insert(make_pair(inputData.points[i],dist));
        }
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();
    return bestPointsDists;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LSH //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


pair<ClassPoint,double> lsh_approximate_NN(ClassPoint q, vector<HashTable> hashTables, LSH_hash_info *hInfo, double &time)
{
    ClassPoint b;
    pair<ClassPoint,double> best;
    best.first = b; //best point-candidate
    best.second = DBL_MAX; //best distance of best candidate

    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();

    int L = hInfo->get_L();
    for (int i = 0; i < L; i++) {
        // Update hinfo with the right vectors for every hash table, to compute query's g-value
        hInfo->update_v(hashTables[i].v);
        hInfo->update_t(hashTables[i].t);
        hInfo->update_r(hashTables[i].r);
        // Find g value for query point.
        vector<int> hValues;
        int k = hInfo->get_k();
        vector<double> vp = q.vpoint;
        for (int j = 0; j < k; j++)
        {
            hValues.push_back(compute_hValue(j, vp, hInfo));
            
        }
        long int ID = compute_IDvalue(hValues, hInfo);
        int g = compute_gValue(ID, hashTables[i].get_bucketsNumber());
        list<HashNode> listToSearch = hashTables[i].get_bucketList(g);
        typename list<HashNode>::iterator current;
        for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
            if (ID != current->ID)
            {
                continue;
            }

            double dist = distance(q.vpoint,current->point->vpoint, 2);
            if (dist < best.second)
            {
                best.second = dist;
                best.first = *(current->point);
            }
        }
    }

    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();

    return best;
    
}


set<pair<ClassPoint,double>, CompDist> lsh_approximate_nNN(ClassPoint q, int N, vector<HashTable> hashTables, LSH_hash_info *hInfo, double &time)
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    // Initialise a set two hold pairs of best point/best distance.
    set<pair<ClassPoint,double>, CompDist> bestPointsDists;
    ClassPoint a;
    bestPointsDists.insert(make_pair(a,DBL_MAX));
    int L = hInfo->get_L();
    for (int i = 0; i < L; i++) {
        // Update hinfo with the right vectors for every hash table, to compute query's g-value
        hInfo->update_v(hashTables[i].v);
        hInfo->update_t(hashTables[i].t);
        hInfo->update_r(hashTables[i].r);
        // Find g value for query point.
        vector<int> hValues;
        int k = hInfo->get_k();
        vector<double> vp = q.vpoint;
        for (int j = 0; j < k; j++)
        {
            hValues.push_back(compute_hValue(j, vp, hInfo));
            
        }
        long int ID = compute_IDvalue(hValues, hInfo);
        int g = compute_gValue(ID, hashTables[i].get_bucketsNumber());
        list<HashNode> listToSearch = hashTables[i].get_bucketList(g);
        typename list<HashNode>::iterator current;
        for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
            if (ID != current->ID)
            {
                continue;
            }
            
            double dist = distance(q.vpoint,current->point->vpoint, 2);
            if (bestPointsDists.size()==N) //if set is full
            {
                //if the biggest distance in set is equal/greater than current distance
                if (prev(bestPointsDists.end())->second >= dist)
                {
                    bestPointsDists.erase(prev(bestPointsDists.end())); //pop biggest distance pair
                    bestPointsDists.insert(make_pair(*(current->point),dist)); //insert new point/distance
                }
            }
            else if (bestPointsDists.size()<N) //if there is space in set insert pair
            {
                bestPointsDists.insert(make_pair(*(current->point),dist));
            }
        }
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();
    return bestPointsDists;
}


unordered_map<string,double> lsh_approximate_range_search(ClassPoint q, double R, vector<HashTable> hashTables, LSH_hash_info *hInfo)
{
    // Initialise an unordered map two hold points-distances inside radius r.
    unordered_map<string,double> rPoints;
    int L = hInfo->get_L();
    for (int i = 0; i < L; i++) {
        // Update hinfo with the right vectors for every hash table, to compute query's g-value
        hInfo->update_v(hashTables[i].v);
        hInfo->update_t(hashTables[i].t);
        hInfo->update_r(hashTables[i].r);
        // Find g value for query point.
        vector<int> hValues;
        int k = hInfo->get_k();
        vector<double> vp = q.vpoint;
        for (int j = 0; j < k; j++)
        {
            hValues.push_back(compute_hValue(j, vp, hInfo));
            
        }
        long int ID = compute_IDvalue(hValues, hInfo);
        int g = compute_gValue(ID, hashTables[i].get_bucketsNumber());
        list<HashNode> listToSearch = hashTables[i].get_bucketList(g);
        typename list<HashNode>::iterator current;
        for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
            if (ID != current->ID)
            {  
                continue;
            }

            double dist = distance(q.vpoint,current->point->vpoint, 2);
            if (dist < R && dist != 0)
            {
                rPoints.insert(make_pair(current->point->itemID,dist));
            }
        }
    }
    return rPoints;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CUBE //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int hammingDistance(int n1, int n2)
{
    int x = n1 ^ n2;
    int setBits = 0;
 
    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}


pair<ClassPoint,double> cube_approximate_NN(ClassPoint q, CubeTable cubeTable, CUBE_hash_info *hInfo)
{
    ClassPoint b;
    pair<ClassPoint,double> best;
    best.first = b; //best point-candidate
    best.second = DBL_MAX; //best distance of best candidate
    
    // Update hinfo with the right vectors for hash table, to compute query's g-value
    hInfo->update_v(cubeTable.v);
    hInfo->update_t(cubeTable.t);
    // Find g value for query point.
    vector<int> hValues;
    int k = hInfo->get_k();
    vector<double> vp = q.vpoint;
    for (int i = 0; i < k; i++)
    {
        hValues.push_back(compute_hValue(i, vp, hInfo));
    }
    vector<int> fValues;
    for (int i = 0; i < k; i++)
    {
        fValues.push_back(compute_fValue(i, hValues[i], hInfo));
    }
    int g = compute_gValue(fValues, hInfo);

    int M = hInfo->get_M(); //maximum number of points to search
    int currentM = 0; //keep track how many points we have searched till now

    // We finished checking vertice g but we have more points to check (because currentM < M).
    // We can't check more than 'probes minus the g' vertices and we will search until hamming distance = maxHD.
    int maxProbes = hInfo->get_probes();
    int maxHD = hInfo->maxHD;
    int numVertices = cubeTable.get_bucketsNumber();
    int currentProbes = maxProbes; //current probes will be reduced every time a new vertex is checked
    
    // Loop hamming distances till maximum hd defined in cube program.
    for (int hd = 0; hd <= maxHD; hd++) 
    {
        vector<int> HDhasProbes; //will consist of all vertices-indexes of current hamming distance
        
        // Loop all vertices and store only the ones with hamming distance = hd from g-index.
        for (int i = 0; i < numVertices; i++)
        {
            if(hammingDistance(g,i) == hd)
            {
                HDhasProbes.push_back(i);
            }
        }
        random_shuffle(HDhasProbes.begin(), HDhasProbes.end()); //randomize vertices
        
        // Loop from 0 to maxProbes and access current hd's probes.
        for (int i = 0; i < maxProbes; i++)
        {
            if (i==HDhasProbes.size()) //current hd has no more probes - break and increment hd
            {
                break;
            }
                        
            // Finally search vertice for best distance.
            list<Vertice> listToSearch = cubeTable.get_bucketList(HDhasProbes[i]);
            typename list<Vertice>::iterator current;
            for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
                double dist = distance(q.vpoint,current->point->vpoint, 2);
                currentM++;
                if (dist < best.second)
                {
                    best.second = dist;
                    best.first = *(current->point);
                }
                if (currentM == M) //maximum points to check reached
                {
                    return best;
                }   
            }

            currentProbes--;
            // If no more probes to check return.
            if (currentProbes == 0)
            {
                return best;
            }
        }
        HDhasProbes.clear(); //clear table to get next hd's probes

    }
    // Max HD reached.
    return best;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


set<pair<ClassPoint,double>, CompDist> cube_approximate_nNN(ClassPoint q, int N, CubeTable cubeTable, CUBE_hash_info *hInfo, double &time)
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    // Initialise a set two hold pairs of best point/best distance.
    set<pair<ClassPoint,double>, CompDist> bestPointsDists;
    ClassPoint a;
    bestPointsDists.insert(make_pair(a,DBL_MAX));
    
    // Update hinfo with the right vectors for hash table, to compute query's g-value
    hInfo->update_v(cubeTable.v);
    hInfo->update_t(cubeTable.t);
    // Find g value for query point.
    vector<int> hValues;
    int k = hInfo->get_k();
    vector<double> vp = q.vpoint;
    for (int i = 0; i < k; i++)
    {
        hValues.push_back(compute_hValue(i, vp, hInfo));
    }
    vector<int> fValues;
    for (int i = 0; i < k; i++)
    {
        fValues.push_back(compute_fValue(i, hValues[i], hInfo));
    }
    int g = compute_gValue(fValues, hInfo);

    int M = hInfo->get_M(); //maximum number of points to search
    int currentM = 0; //keep track how many points we have searched till now

    // We finished checking vertice g but we have more points to check (because currentM < M).
    // We can't check more than 'probes minus the g' vertices and we will search until hamming distance = maxHD.
    int maxProbes = hInfo->get_probes();
    int maxHD = hInfo->maxHD;
    int numVertices = cubeTable.get_bucketsNumber();
    int currentProbes = maxProbes; //current probes will be reduced every time a new vertex is checked
    
    // Loop hamming distances till maximum hd defined in cube program.
    for (int hd = 0; hd <= maxHD; hd++) 
    {
        vector<int> HDhasProbes; //will consist of all vertices-indexes of current hamming distance
        
        // Loop all vertices and store only the ones with hamming distance = hd from g-index.
        for (int i = 0; i < numVertices; i++)
        {
            if(hammingDistance(g,i) == hd)
            {
                HDhasProbes.push_back(i);
            }
        }
        random_shuffle(HDhasProbes.begin(), HDhasProbes.end()); //randomize vertices
        
        // Loop from 0 to maxProbes and access current hd's probes.
        for (int i = 0; i < maxProbes; i++)
        {
            if (i==HDhasProbes.size()) //current hd has no more probes - break and increment hd
            {
                break;
            }

            // Finally search vertice for best distance.
            list<Vertice> listToSearch = cubeTable.get_bucketList(HDhasProbes[i]);
            typename list<Vertice>::iterator current;
            for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
                double dist = distance(q.vpoint,current->point->vpoint, 2);
                currentM++;
                if (bestPointsDists.size()==N) //if set is full
                {
                    //if the biggest distance in set is equal/greater than current distance
                    if (prev(bestPointsDists.end())->second >= dist)
                    {
                        bestPointsDists.erase(prev(bestPointsDists.end())); //pop biggest distance pair
                        bestPointsDists.insert(make_pair(*(current->point),dist)); //insert new point/distance
                    }
                }
                else if (bestPointsDists.size()<N) //if there is space in set insert pair
                {
                    bestPointsDists.insert(make_pair(*(current->point),dist));
                }
                if (currentM == M) //maximum points to check reached
                {
                    high_resolution_clock::time_point end = high_resolution_clock::now();
                    duration<double> time_span = duration_cast<duration<double>>(end - start);
                    time = time_span.count();
                    return bestPointsDists;
                }   
            }

            currentProbes--;
            // If no more probes to check return.
            if (currentProbes == 0)
            {
                high_resolution_clock::time_point end = high_resolution_clock::now();
                duration<double> time_span = duration_cast<duration<double>>(end - start);
                time = time_span.count();
                return bestPointsDists;
            }
        }
        HDhasProbes.clear(); //clear table to get next hd's probes

    }
    // Max HD reached.
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();
    return bestPointsDists;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unordered_map<string,double> cube_approximate_range_search(ClassPoint q, double R, CubeTable cubeTable, CUBE_hash_info *hInfo)
{
    // Initialise an unordered map two hold points-distances inside radius r.
    unordered_map<string,double> rPoints;

    // Update hinfo with the right vectors for hash table, to compute query's g-value
    hInfo->update_v(cubeTable.v);
    hInfo->update_t(cubeTable.t);
    // Find g value for query point.
    vector<int> hValues;
    int k = hInfo->get_k();
    vector<double> vp = q.vpoint;
    for (int i = 0; i < k; i++)
    {
        hValues.push_back(compute_hValue(i, vp, hInfo));
    }
    vector<int> fValues;
    for (int i = 0; i < k; i++)
    {
        fValues.push_back(compute_fValue(i, hValues[i], hInfo));
    }
    int g = compute_gValue(fValues, hInfo);

    int M = hInfo->get_M(); //maximum number of points to search
    int currentM = 0; //keep track how many points we have searched till now

    // We finished checking vertice g but we have more points to check (because currentM < M).
    // We can't check more than 'probes minus the g' vertices and we will search until hamming distance = maxHD.
    int maxProbes = hInfo->get_probes();
    int maxHD = hInfo->maxHD;
    int numVertices = cubeTable.get_bucketsNumber();
    int currentProbes = maxProbes; //current probes will be reduced every time a new vertex is checked
    
    // Loop hamming distances till maximum hd defined in cube program.
    for (int hd = 0; hd <= maxHD; hd++) 
    {
        vector<int> HDhasProbes; //will consist of all vertices-indexes of current hamming distance
        
        // Loop all vertices and store only the ones with hamming distance = hd from g-index.
        for (int i = 0; i < numVertices; i++)
        {
            if(hammingDistance(g,i) == hd)
            {
                HDhasProbes.push_back(i);
            }
        }
        random_shuffle(HDhasProbes.begin(), HDhasProbes.end()); //randomize vertices
        
        // Loop from 0 to maxProbes and access current hd's probes.
        for (int i = 0; i < maxProbes; i++)
        {
            if (i==HDhasProbes.size()) //current hd has no more probes - break and increment hd
            {
                break;
            }
                        
            // Finally search vertice for best distance.
            list<Vertice> listToSearch = cubeTable.get_bucketList(HDhasProbes[i]);
            typename list<Vertice>::iterator current;
            for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
                double dist = distance(q.vpoint,current->point->vpoint, 2);
                currentM++;
                for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
                    double dist = distance(q.vpoint,current->point->vpoint, 2);
                    if (dist < R && dist !=0)
                    {
                        rPoints.insert(make_pair(current->point->itemID,dist));
                    }
                }
                if (currentM == M) //maximum points to check reached
                {
                    return rPoints;
                }   
            }

            currentProbes--;
            // If no more probes to check return.
            if (currentProbes == 0)
            {
                return rPoints;
            }
        }
        HDhasProbes.clear(); //clear table to get next hd's probes

    }
    // Max HD reached.
    return rPoints;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DISCRETE - CONTINUOUS FRECHET //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


pair<ClassCurve,double> lsh_approximate_NN(ClassCurve q, vector<GridTable> gridTables, LSH_hash_info *hInfo, double &time)
{
    ClassCurve c;
    pair<ClassCurve,double> best;
    best.first = c; //best point-candidate
    best.second = DBL_MAX; //best distance of best candidate

    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();

    int L = hInfo->get_L();

    if (gridTables[0].pointDim == 1) // CONTINUOUS FRECHET
    {
        for (int i = 0; i < L; i++) 
        {
            filtering(&q, gridTables[i].epsilon);
            ClassCurve gridCurve = snapTo1dGrid(q,gridTables[i].tShiftGrid,gridTables[i].delta);
            minima_maxima(&gridCurve);
            padding(&gridCurve, gridTables[i].curveDim);
            vector<double> LSHvector = keyLSHvector1D(gridCurve);
            // Update hinfo with the right vectors for every hash table, to compute query's g-value
            hInfo->update_v(gridTables[i].v);
            hInfo->update_t(gridTables[i].t);
            hInfo->update_r(gridTables[i].r);
            // Find g value for query point.
            vector<int> hValues;
            int k = hInfo->get_k();
            for (int j = 0; j < k; j++)
            {
                hValues.push_back(compute_hValue(j, LSHvector, hInfo));
                
            }
            long int ID = compute_IDvalue(hValues, hInfo);
            int g = compute_gValue(ID, gridTables[i].get_bucketsNumber());
            list<GridNode> listToSearch = gridTables[i].get_bucketList(g);
            typename list<GridNode>::iterator current;

            // Convert query curve q to be conpatible with Frechet::Continuous::distance
            Curve QueryFredCurve = cont_convert_curve(q);

            for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) 
            {
                if (ID != current->ID)
                {
                    continue;
                }
                
                // Convert current bucket curve to be conpatible with Frechet::Continuous::distance
                Curve CurrentFredCurve = cont_convert_curve(*current->curve);
                
                Frechet::Continuous::Distance d = Frechet::Continuous::distance(QueryFredCurve,CurrentFredCurve);
                double dist = d.value;
                if (dist < best.second)
                {
                    best.second = dist;
                    best.first = *(current->curve);
                }
            }
            if (best.second == DBL_MAX)
            {
                for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) 
                {
                    // Convert current bucket curve to be conpatible with Frechet::Continuous::distance
                    Curve CurrentFredCurve = cont_convert_curve(*current->curve);
                    
                    Frechet::Continuous::Distance d = Frechet::Continuous::distance(QueryFredCurve,CurrentFredCurve);
                    double dist = d.value;
                    if (dist < best.second)
                    {
                        best.second = dist;
                        best.first = *(current->curve);
                    }
                }
            }
        }
    }
    else if (gridTables[0].pointDim == 2) // DISCRETE FRECHET
    {
        for (int i = 0; i < L; i++) 
        {
            ClassCurve grid_curve = snapTo2dGrid(q,gridTables[i].tShiftGrid,gridTables[i].delta);
            padding(&grid_curve, gridTables[i].curveDim);
            vector<double> LSHvector = keyLSHvector2D(grid_curve);
            // Update hinfo with the right vectors for every hash table, to compute query's g-value
            hInfo->update_v(gridTables[i].v);
            hInfo->update_t(gridTables[i].t);
            hInfo->update_r(gridTables[i].r);
            // Find g value for query point.
            vector<int> hValues;
            int k = hInfo->get_k();
            for (int j = 0; j < k; j++)
            {
                hValues.push_back(compute_hValue(j, LSHvector, hInfo));
                
            }
            long int ID = compute_IDvalue(hValues, hInfo);
            int g = compute_gValue(ID, gridTables[i].get_bucketsNumber());
            list<GridNode> listToSearch = gridTables[i].get_bucketList(g);
            typename list<GridNode>::iterator current;
            for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current )
            {
                if (ID != current->ID)
                {
                    continue;
                }

                double dist = discrete_frechet_distance(q,*(current->curve));
                // cout << dist << endl;
                if (dist < best.second)
                {
                    best.second = dist;
                    best.first = *(current->curve);
                }
            }
            if (best.second == DBL_MAX)
            {
                for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current )
                {
                    double dist = discrete_frechet_distance(q,*(current->curve));
                    // cout << dist << endl;
                    if (dist < best.second)
                    {
                        best.second = dist;
                        best.first = *(current->curve);
                    }
                }
            }
            
        }
    }
    else
    {
        cout << "Unexpected error occured: ClassPoint dimension must be 1 for continuous / 2 for discrete frechet distance." << endl;
        exit (EXIT_FAILURE);
    }

    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();
    
    return best;
}


pair<ClassCurve,double> true_NN(ClassCurve q, Vector_of_curves inputData, double &time)
{
    ClassCurve b;
    pair<ClassCurve,double> best;
    best.first = b; //best Curve-candidate
    best.second = DBL_MAX; //best distance of best candidate

    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();

    if (q.cpoints[0].vpoint.size() == 1)
    {
        // Convert query curve q to be conpatible with Frechet::Continuous::distance
        Curve QueryFredCurve = cont_convert_curve(q);

        for (int i = 0; i < inputData.curves.size(); i++)
        {
            // Convert data curve to be conpatible with Frechet::Continuous::distance
            Curve DataFredCurve = cont_convert_curve(inputData.curves[i]);

            Frechet::Continuous::Distance d = Frechet::Continuous::distance(QueryFredCurve,DataFredCurve);
            double dist = d.value;
            if (dist < best.second)
            {
                best.second = dist;
                best.first = inputData.curves[i];
            }
            
        }
    }
    else if (q.cpoints[0].vpoint.size() == 2)
    {
        for (int i = 0; i < inputData.curves.size(); i++)
        {
            double dist = discrete_frechet_distance(q,inputData.curves[i]);
            if (dist < best.second)
            {
                best.second = dist;
                best.first = inputData.curves[i];
            }
            
        }
    }

    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    time = time_span.count();

    return best;

}


unordered_map<string,double> lsh_approximate_range_search(ClassCurve q, double R, vector<GridTable> gridTables, LSH_hash_info *hInfo)
{
    // Initialise an unordered map two hold points-distances inside radius r.
    unordered_map<string,double> rCurves;
    
    int L = hInfo->get_L();
    for (int i = 0; i < L; i++) {

        ClassCurve grid_curve = snapTo2dGrid(q,gridTables[i].tShiftGrid,gridTables[i].delta);
        padding(&grid_curve, gridTables[i].curveDim);
        vector<double> LSHvector = keyLSHvector2D(grid_curve);

        // Update hinfo with the right vectors for every hash table, to compute query's g-value
        hInfo->update_v(gridTables[i].v);
        hInfo->update_t(gridTables[i].t);
        hInfo->update_r(gridTables[i].r);
        // Find g value for query point.
        vector<int> hValues;
        int k = hInfo->get_k();
        for (int j = 0; j < k; j++)
        {
            hValues.push_back(compute_hValue(j, LSHvector, hInfo));
            
        }
        long int ID = compute_IDvalue(hValues, hInfo);
        int g = compute_gValue(ID, gridTables[i].get_bucketsNumber());
        list<GridNode> listToSearch = gridTables[i].get_bucketList(g);
        typename list<GridNode>::iterator current;
        for (current = listToSearch.begin() ; current != listToSearch.end() ; ++current ) {
            if (ID != current->ID)
            {  
                continue;
            }

            double dist = discrete_frechet_distance(q,*(current->curve));
            if (dist < R && dist != 0)
            {
                rCurves.insert(make_pair(current->curve->curveID,dist));
            }
        }
    }
    return rCurves;
}
