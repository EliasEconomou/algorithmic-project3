#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include "../../include/search_methods.hpp"


using namespace std;


int main(int argc, char** argv) {

    string inputFile="0", queryFile="0", outputFile="0", algorithm="NULL", metric="NULL";
    int k=-1; //k will get a default value later according to algorithm given
    int L=5, M=10, probes=2;
    double delta=-1;
    bool FrechetBrute=0; //0 computes true distance in frechet algorithms, 1 does not. (default 0 - compute)
    string word;

    // for (int i = 1; i < argc; i+=2) //argument options always start with '-'
    // {
    //     if ((argv[i][0]!='-') || (argv[i+1]==NULL))
    //     {
    //         cout << "Error: Wrong argument." << endl;
    //         return -1;
    //     }
    // }
    
    //---Parse arguments---
    for (int i = 1; i < argc; i++){
        word  = argv[i];
        if(word == "-i"){
            if(i+1 < argc){
                inputFile = argv[++i];
            }else{
                cout << "No argument for input file." << endl;
                return -1;
            }
        }
        else if (word == "-q"){
            if(i+1 <= argc){
                queryFile = argv[++i];
            }else{
                cout << "No argument for query file." << endl;
                return -1;
            }
        }
        else if(word == "-k"){
            if(i+1 <= argc){
                k = stoi(argv[++i]);
            }else{
                cout << "Error: No argument for k value." << endl;
                return -1;
            }
        }
        else if(word == "-L"){
            if(i+1 <= argc){
                L = stoi(argv[++i]);
            }else{
                cout << "Error: No argument for L value." << endl;
                return -1;
            }
        }
        else if(word == "-M"){
            if(i+1 <= argc){
                M = stoi(argv[++i]);
            }else{
                cout << "Error: No argument for M value." << endl;
                return -1;
            }
        }
        else if(word == "-probes"){
            if(i+1 <= argc){
                probes = stoi(argv[++i]);
            }else{
                cout << "Error: No argument for probes value." << endl;
                return -1;
            }
        }
        else if (word == "-o"){
            if(i+1 <= argc){
                outputFile = argv[++i];
            }else{
                cout << "Error: No argument for output file." << endl;
                return -1;
            }
        }
        else if (word == "-algorithm"){
            if(i+1 <= argc){
                algorithm = argv[++i];
                if ((algorithm!="LSH") && (algorithm!="Hypercube") && (algorithm!="Frechet"))
                {
                    cout << "Error: Wrong algorithm argument given." << endl;
                    return -1;
                }
                
            }else{
                cout << "Error: No argument for algorithm method." << endl;
                return -1;
            }
        }
        else if(word == "-metric"){
            if(i+1 <= argc){
                metric = argv[++i];
                if ((metric!="discrete") && (metric!="continuous"))
                {
                    cout << "Error: Wrong metric argument given." << endl;
                    return -1;
                }
            }
        }
        else if(word == "-delta"){
            if(i+1 <= argc){
                delta = stod(argv[++i]);
            }else{
                cout << "Error: No argument for delta value." << endl;
                return -1;
            }
        }
        else if (word == "-nobrute")
        {
            FrechetBrute = 1;
        }
        
        else{

            cout << "Error: Wrong argument: " << word << endl;
            return -1;
        }
    }
    if (algorithm == "NULL")
    {
        cout << "Error: Algorithm is required." << endl;
        return -1;
    }
    
    if (algorithm == "Frechet") //Frechet must be discrete/continuous and have a delta value
    {   if (metric == "NULL")
        {
            cout << "Error: Metric is required for Frechet." << endl;
            return -1;
        }
        if (delta == -1)
        {
            cout << "Error: Delta is required for Frechet." << endl;
            return -1;
        }
    }
    if (((algorithm == "LSH") || (algorithm == "Hypercube")) && (metric != "NULL"))
    {
        cout << "Error: LSH/Hypercube can't have Frechet metric." << endl;
        return -1;
    }
    if (k==-1) //k not given - use default
    {
        if ((algorithm=="LSH") || (algorithm=="Frechet"))
            k=4;
        else if (algorithm=="Hypercube")
            k=14;
    }
    
    // Print arguments
    cout << "inputFile = " << inputFile << ". queryFile = " << queryFile << ". k = " << k << ". L = " << L;
    cout << ". M = " << M << ". probes = " << probes << ". outputFile = " << outputFile << ". algorithm = " << algorithm;
    cout << ". metric = " << metric << ". delta = " << delta << "." << endl;


    // Call appropriate function based on algorithm given.
    if (algorithm == "LSH")
    {
        time_series_LSH(inputFile, queryFile, outputFile, k, L);
    }
    else if (algorithm == "Hypercube")
    {
        time_series_Hypercube(inputFile, queryFile, outputFile, k, M, probes);
    }
    else if (algorithm == "Frechet")
    {
        if (metric == "discrete")
        {
            time_series_DiscreteFrechet(inputFile, queryFile, outputFile, k, L, delta, FrechetBrute);
        }
        else if (metric == "continuous")
        {
            time_series_ContinuousFrechet(inputFile, queryFile, outputFile, k, L, delta, FrechetBrute);
        }
    }
    
    

}