#include "../include/cube_table.hpp"

using namespace std;


Vertice::Vertice(ClassPoint *p)
{
  this->point = p;
}


CubeTable::CubeTable(int bucketsNumber)
{
  this->bucketsNumber = bucketsNumber;
  this->lists.resize(this->bucketsNumber);
}


void CubeTable::CTinsert(ClassPoint *p, CUBE_hash_info *hInfo)
{
  vector<int> hValues;
  hInfo->update_v(this->v);
  hInfo->update_t(this->t);
  int k = hInfo->get_k();
  vector<double> vp = p->vpoint;
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
  lists[g].push_back(Vertice(p));
}


int CubeTable::get_bucketsNumber()
{
  return this->bucketsNumber;
}


list<Vertice> CubeTable::get_bucketList(int g)
{
  return this->lists[g];
}