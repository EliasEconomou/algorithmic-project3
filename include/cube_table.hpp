#ifndef CUBE_TABLE_H
#define CUBE_TABLE_H

#include <iostream>
#include <list>
#include <vector>
#include "hash_functions.hpp"
#include "point_functions.hpp"


class Vertice
{
public:
    ClassPoint *point;
    Vertice(ClassPoint *p);
};

class CubeTable
{
    int bucketsNumber;
    std::vector<std::list<Vertice>> lists; //a list for every bucket
    
public:
    std::vector<double> t;
    std::vector<std::vector<double> > v; //k vectors to use to compute every h

    CubeTable(int bucketsNumber);
    
    // Insert item in cube table
    void CTinsert(ClassPoint *p, CUBE_hash_info *hInfo);

    int get_bucketsNumber();

    std::list<Vertice> get_bucketList(int g);
};


#endif