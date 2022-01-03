#include "../include/hash_table.hpp"

using namespace std;


HashNode::HashNode(ClassPoint *p, long int ID)
{
  this->ID = ID;
  this->point = p;
}


HashTable::HashTable(int bucketsNumber)
{
  this->bucketsNumber = bucketsNumber;
  this->lists.resize(this->bucketsNumber);
}


void HashTable::HTinsert(ClassPoint *p, LSH_hash_info *hInfo)
{
  vector<int> hValues;
  hInfo->update_v(this->v);
  hInfo->update_t(this->t);
  hInfo->update_r(this->r);
  int k = hInfo->get_k();
  vector<double> vp = p->vpoint;
  // cout << p->itemID << "  ";
  for (int i = 0; i < k; i++)
  {
    hValues.push_back(compute_hValue(i, vp, hInfo));
  }

  long int ID = compute_IDvalue(hValues, hInfo);
  int g = compute_gValue(ID, this->bucketsNumber);
  lists[g].push_back(HashNode(p, ID));
}


void HashTable::HTdisplay() 
{
  for (int k=0 ; k < bucketsNumber ; k++){
    cout << "In bucket #" << k << " of hashtable: \n" ;
    typename list<HashNode>::iterator current;
      for (current = lists[k].begin() ; current != lists[k].end() ; ++current ){
        // std::cout << current->point->size() << endl;
        for (auto j = current->point->vpoint.begin() ; j != current->point->vpoint.end() ; ++j){
          cout << *j << " ";
        }
        cout << endl;
      }
    cout << endl;
  }
}


int HashTable::get_bucketsNumber()
{
  return this->bucketsNumber;
}


list<HashNode> HashTable::get_bucketList(int g)
{
  return this->lists[g];
}