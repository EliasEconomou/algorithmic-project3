#include "../include/search_methods.hpp"
#define MAX_HAMMING_DIST 8

using namespace std;


// Function for i) assigment using LSH.
void time_series_LSH(string inputFile, string queryFile, string outputFile, int k, int L)
{
    cout << "Executing LSH for vectors." << endl;
    if (inputFile == "0")
    {
        cout << "Give path to input file: ";
        cin >> inputFile;
    }

    Vector_of_points inputData;
    inputData = parsing(inputFile);
    
    int vectorsNumber = inputData.points.size();
    int dimension = inputData.points[0].vpoint.size();
    int bucketsNumber = vectorsNumber/8;
    LSH_hash_info hInfo(k, dimension, L);

    vector<HashTable> hashTables;
    for (int i = 0; i < L; i++)
    {
        HashTable ht(bucketsNumber);
        hashTables.push_back(ht);
    }
    
    for (int i = 0; i < L; i++)
    {
        hashTables[i].v = compute_v(k,dimension);
        hashTables[i].t = compute_t(k);
        hashTables[i].r = compute_r(k);
        for (int j = 0; j < vectorsNumber; j++)
        {
            hashTables[i].HTinsert(&inputData.points[j], &hInfo);
        }
    }

    if (queryFile == "0")
    {
        cout << "Give path to query file: ";
        cin >> queryFile;
    }

    if (outputFile == "0")
    {
        cout << "Give path to output file: ";
        cin >> outputFile;
    }

    Vector_of_points queryData;
    queryData = parsing(queryFile);

    ofstream out (outputFile);

    std::cout << "Writing to output file..." << endl;

    double lshSumTime = 0.0;
    double trueSumTime = 0.0;
    double lshTime,trueTime;
    double MAF = 0.0;
    for (int i = 0; i < queryData.points.size(); i++)
    {
        out << "Query: " << queryData.points[i].itemID << endl;
        out << "Algorithm: " << "LSH_Vector" << endl;
        pair<ClassPoint,double> lshBestCurveDist;
        pair<ClassPoint,double> trueBestCurveDist;
        lshTime = 0.0;
        trueTime = 0.0;

        lshBestCurveDist = lsh_approximate_NN(queryData.points[i], hashTables, &hInfo, lshTime);
        lshSumTime += lshTime;
        // cout << "Lsh time " << lshTime << ". Now total lsh time is " << lshSumTime << "." << endl;
        out << "Approximate Nearest neighbor: " << lshBestCurveDist.first.itemID << endl;
        
            trueBestCurveDist = true_NN(queryData.points[i],inputData,trueTime);
            trueSumTime += trueTime;
            // cout << "True time " << trueTime << ". Now total true time is " << trueSumTime << "." << endl;
            out << "True Nearest neighbor: " << trueBestCurveDist.first.itemID << endl;
        
        out << "distanceApproximate: " << lshBestCurveDist.second << endl;
        out << "distanceTrue: " << trueBestCurveDist.second << endl;
        
        out << endl;
        
        double AF = lshBestCurveDist.second/trueBestCurveDist.second;
        if (MAF <= AF)
        {
            MAF = AF;
        }
    }
    
    out << endl;
    out << "tApproximateAverage: " << lshSumTime/queryData.points.size() << endl;
    out << "tTrueAverage: " << trueSumTime/queryData.points.size() << endl;
    out << "MAF: " << MAF << endl;

    out << endl;
    out.close();

    std::cout << "Operation completed successfully." << endl << "Exiting." << endl;
    
}

// Function for i) assigment using Hypercube.
void time_series_Hypercube(string inputFile, string queryFile, string outputFile, int k, int M, int probes)
{
    cout << "Executing Hypercube." << endl;
    if (inputFile == "0")
    {
        cout << "Give path to input file: ";
        cin >> inputFile;
    }
    
    Vector_of_points inputData;
    inputData = parsing(inputFile);
    
    int vectorsNumber = inputData.points.size();
    int dimension = inputData.points[0].vpoint.size();
    int bucketsNumber = pow(2,k);
    CUBE_hash_info hInfo(k, dimension, M, probes, MAX_HAMMING_DIST);

    CubeTable cubeTable(bucketsNumber);
    cubeTable.v = compute_v(k,dimension);
    cubeTable.t = compute_t(k);
    for (int i = 0; i < vectorsNumber; i++)
    {
        cubeTable.CTinsert(&inputData.points[i], &hInfo);
    }

    if (queryFile == "0")
    {
        cout << "Give path to query file: ";
        cin >> queryFile;
    }

    if (outputFile == "0")
    {
        cout << "Give path to output file: ";
        cin >> outputFile;
    }
    
    Vector_of_points queryData;
    queryData = parsing(queryFile);

    ofstream out (outputFile);

    std::cout << "Writing to output file..." << endl;

    double cubeSumTime = 0.0;
    double trueSumTime = 0.0;
    double cubeTime,trueTime;
    double MAF = 0.0;
    for (int i = 0; i < queryData.points.size(); i++)
    {
        out << "Query: " << queryData.points[i].itemID << endl;
        out << "Algorithm: " << "Hypercube" << endl;
        set<pair<ClassPoint,double>, CompDist> cubeBestPointsDists;
        set<pair<ClassPoint,double>, CompDist> trueBestPointsDists;
        cubeTime = 0.0;
        trueTime = 0.0;

        cubeBestPointsDists = cube_approximate_nNN(queryData.points[i], 1, cubeTable, &hInfo, cubeTime);
        cubeSumTime += cubeTime;
        trueBestPointsDists = true_nNN(queryData.points[i], 1, inputData, trueTime);
        trueSumTime += trueTime;
        
        auto it1 = cubeBestPointsDists.begin();
        auto it2 = trueBestPointsDists.begin();
        out << "Approximate Nearest neighbor: " << it1->first.itemID << endl;
        out << "True Nearest neighbor: " << it2->first.itemID << endl;
        out << "distanceApproximate: " << it1->second << endl;
        out << "distanceTrue: " << it2->second << endl;

        out << endl;
        
        double AF = it1->second/it2->second;
        if (MAF <= AF)
        {
            MAF = AF;
        }
    }

    out << endl;
    out << "tApproximateAverage: " << cubeSumTime/queryData.points.size() << endl;
    out << "tTrueAverage: " << trueSumTime/queryData.points.size() << endl;
    out << "MAF: " << MAF << endl;

    out << endl;
    out.close();

    std::cout << "Operation completed successfully." << endl << "Exiting." << endl;

}

// Function for ii) assigment using Discrete Ferchet.
void time_series_DiscreteFrechet(string inputFile, string queryFile, string outputFile, int k, int L, double delta, bool FrechetBrute)
{
    cout << "Executing Discrete Frechet for curves." << endl;
    if (inputFile == "0")
    {
        cout << "Give path to input file: ";
        cin >> inputFile;
    }

    std::cout << "Parsing data file..." << endl;
    Vector_of_curves inputData;
    inputData = curve_parsing(inputFile, 2);

    int curveDim = inputData.curves[0].cpoints.size(); // number of points in curves
    int pointDim = inputData.curves[0].cpoints[0].vpoint.size(); // dimension of points in curve
    int curvesNumber = inputData.curves.size();
    int bucketsNumber = curvesNumber/8;
    
    LSH_hash_info hInfo(k, curveDim, L);

    vector<GridTable> gridTables;
    for (int i = 0; i < L; i++)
    {
        GridTable gt(bucketsNumber, delta, curveDim, pointDim);
        gridTables.push_back(gt);
    }
    
    std::cout << "Storing data in hash tables..." << endl;
    for (int i = 0; i < L; i++)
    {
        gridTables[i].v = compute_v(k,2*curveDim);
        gridTables[i].t = compute_t(k);
        gridTables[i].r = compute_r(k);
        for (int j = 0; j < curvesNumber; j++)
        {
            gridTables[i].GridInsert(&inputData.curves[j], &hInfo);
        }
    }

    if (queryFile == "0")
    {
        cout << "Give path to query file: ";
        cin >> queryFile;
    }

    if (outputFile == "0")
    {
        cout << "Give path to output file: ";
        cin >> outputFile;
    }

    Vector_of_curves queryData;
    queryData = curve_parsing(queryFile, pointDim);

    ofstream out (outputFile);

    std::cout << "Writing to output file..." << endl;

    double lshSumTime = 0.0;
    double trueSumTime = 0.0;
    double lshTime,trueTime;
    double MAF = 0.0;
    for (int i = 0; i < queryData.curves.size(); i++)
    {
        out << "Query: " << queryData.curves[i].curveID << endl;
        out << "Algorithm: " << "LSH_Frechet_Discrete" << endl;
        pair<ClassCurve,double> lshBestCurveDist;
        pair<ClassCurve,double> trueBestCurveDist;
        lshTime = 0.0;
        trueTime = 0.0;

        lshBestCurveDist = lsh_approximate_NN(queryData.curves[i], gridTables, &hInfo, lshTime);
        lshSumTime += lshTime;
        // cout << "Lsh time " << lshTime << ". Now total lsh time is " << lshSumTime << "." << endl;
        out << "Approximate Nearest neighbor: " << lshBestCurveDist.first.curveID << endl;
        
        if (FrechetBrute == 0)
        {
            trueBestCurveDist = true_NN(queryData.curves[i],inputData,trueTime);
            trueSumTime += trueTime;
            // cout << "True time " << trueTime << ". Now total true time is " << trueSumTime << "." << endl;
            out << "True Nearest neighbor: " << trueBestCurveDist.first.curveID << endl;
        }
        
        out << "distanceApproximate: " << lshBestCurveDist.second << endl;
        if (FrechetBrute == 0)
        {
            out << "distanceTrue: " << trueBestCurveDist.second << endl;
        }
        out << endl;
        
        if (FrechetBrute == 0)
        {
            double AF = lshBestCurveDist.second/trueBestCurveDist.second;
            if (MAF <= AF)
            {
                MAF = AF;
            }
        }
    }
    
    out << endl;
    out << "tApproximateAverage: " << lshSumTime/queryData.curves.size() << endl;
    if (FrechetBrute == 0)
    {
        out << "tTrueAverage: " << trueSumTime/queryData.curves.size() << endl;
        out << "MAF: " << MAF << endl;
    }
    out << endl;
    out.close();

    std::cout << "Operation completed successfully." << endl << "Exiting." << endl;

}

// Function for iii) assigment using Continuous Frechet.
void time_series_ContinuousFrechet(string inputFile, string queryFile, string outputFile, int k, int L, double delta, bool FrechetBrute)
{
    cout << "Executing Continuous Frechet for curves." << endl;
    if (inputFile == "0")
    {
        cout << "Give path to input file: ";
        cin >> inputFile;
    }

    std::cout << "Parsing data file..." << endl;
    Vector_of_curves inputData;
    inputData = curve_parsing(inputFile, 1);

    int curveDim = inputData.curves[0].cpoints.size(); // number of points in curves
    int pointDim = inputData.curves[0].cpoints[0].vpoint.size(); // dimension of points in curve
    int curvesNumber = inputData.curves.size();
    int bucketsNumber = curvesNumber/8;
    
    LSH_hash_info hInfo(k, curveDim, L);

    vector<GridTable> gridTables;
    for (int i = 0; i < L; i++)
    {
        GridTable gt(bucketsNumber, delta, curveDim, pointDim);
        gridTables.push_back(gt);
    }
    
    std::cout << "Storing data in hash tables..." << endl;
    for (int i = 0; i < L; i++)
    {
        gridTables[i].v = compute_v(k,curveDim);
        gridTables[i].t = compute_t(k);
        gridTables[i].r = compute_r(k);
        for (int j = 0; j < curvesNumber; j++)
        {
            gridTables[i].GridInsert(&inputData.curves[j], &hInfo);
        }
    }

    if (queryFile == "0")
    {
        cout << "Give path to query file: ";
        cin >> queryFile;
    }

    if (outputFile == "0")
    {
        cout << "Give path to output file: ";
        cin >> outputFile;
    }

    Vector_of_curves queryData;
    queryData = curve_parsing(queryFile, pointDim);

    ofstream out (outputFile);

    double lshSumTime = 0.0;
    double trueSumTime = 0.0;
    double lshTime,trueTime;
    double MAF = 0.0;
    for (int i = 0; i < queryData.curves.size(); i++)
    {
        out << "Query: " << queryData.curves[i].curveID << endl;
        out << "Algorithm: " << "LSH_Frechet_Continuous" << endl;
        pair<ClassCurve,double> lshBestCurveDist;
        pair<ClassCurve,double> trueBestCurveDist;
        lshTime = 0.0;
        trueTime = 0.0;

        lshBestCurveDist = lsh_approximate_NN(queryData.curves[i], gridTables, &hInfo, lshTime);
        lshSumTime += lshTime;
        // cout << "Lsh time " << lshTime << ". Now total lsh time is " << lshSumTime << "." << endl;
        out << "Approximate Nearest neighbor: " << lshBestCurveDist.first.curveID << endl;
        
        if (FrechetBrute == 0)
        {
            trueBestCurveDist = true_NN(queryData.curves[i],inputData,trueTime);
            trueSumTime += trueTime;
            // cout << "True time " << trueTime << ". Now total true time is " << trueSumTime << "." << endl;
            out << "True Nearest neighbor: " << trueBestCurveDist.first.curveID << endl;
        }
        
        out << "distanceApproximate: " << lshBestCurveDist.second << endl;
        if (FrechetBrute == 0)
        {
            out << "distanceTrue: " << trueBestCurveDist.second << endl;
        }
        out << endl;
        
        if (FrechetBrute == 0)
        {
            double AF = lshBestCurveDist.second/trueBestCurveDist.second;
            if (MAF <= AF)
            {
                MAF = AF;
            }
        }
    }
    
    out << endl;
    out << "tApproximateAverage: " << lshSumTime/queryData.curves.size() << endl;
    if (FrechetBrute == 0)
    {
        out << "tTrueAverage: " << trueSumTime/queryData.curves.size() << endl;
        out << "MAF: " << MAF << endl;
    }
    out << endl;
    out.close();

    std::cout << "Operation completed successfully." << endl << "Exiting." << endl;

}