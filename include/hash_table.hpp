#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <iostream>
#include <list>
#include <vector>
#include "hash_functions.hpp"
#include "point_functions.hpp"

class ClassPoint;

class HashNode
{
public:
    long int ID;
    ClassPoint *point;
    HashNode(ClassPoint *p, long int ID);
};

class HashTable
{
    int bucketsNumber;
    std::vector<std::list<HashNode>> lists; //a list for every bucket
    
public:
    std::vector<double> t;
    std::vector<std::vector<double> > v; //k vectors to use to compute every h
    std::vector<int> r;

    HashTable(int bucketsNumber);
    
    // Insert item in hash table
    void HTinsert(ClassPoint *p, LSH_hash_info *hInfo);

    // Display hash table (debug)
    void HTdisplay();

    int get_bucketsNumber();

    std::list<HashNode> get_bucketList(int g);
};


#endif