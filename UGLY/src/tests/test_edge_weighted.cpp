
#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <set>
#include <unordered_map>

#include "../../include/ugly/edge_weighted.hpp"

using namespace ugly;
using namespace std;

int main(){

    cout << "Beginning Test" << endl;
    cout << "Testing: Constructor" << endl;
    { 
        EdgeWeighted edgeweighted(2,3);
    }

    cout << "Testing: getClassType" << endl;
    {
        assert(EdgeWeighted::getClassType()==constants::EdgeType::weighted);
    }

    cout << "Testing: getEdgeType" << endl;
    {
        EdgeWeighted edgeweighted;
        assert(edgeweighted.getEdgeType()==constants::EdgeType::weighted);
    }

    cout << "Testing: directional" << endl;
    {
      EdgeWeighted edgeweighted;
      assert(!edgeweighted.directional());
    }

    cout << "Testing: operator << " << endl;
    {
        EdgeWeighted edgeweighted(2,3);
        cout << edgeweighted << endl; 
    }

    cout << "Testing: ==" << endl;
    {
        EdgeWeighted edgeweighted(2,3);
        EdgeWeighted edgeweighted2(2,4);
        EdgeWeighted edgeweighted3(3,2);
        assert(edgeweighted==edgeweighted);
        assert((edgeweighted==edgeweighted2)==0);
        assert(edgeweighted==edgeweighted3);
    }

    cout << "Testing: !=" << endl;
    {
        EdgeWeighted edgeweighted(2,3);
        EdgeWeighted edgeweighted2(2,4);
        assert(edgeweighted!=edgeweighted2);
        assert((edgeweighted!=edgeweighted)==0);
    }

    cout << "Testing: getVertex1 and getVertex2" << endl;
    {
        EdgeWeighted edgeweighted(2,3);
        assert(edgeweighted.getVertex1()==2);
        assert(edgeweighted.getVertex2()==3);
    }

    cout << "Testing: < " << endl;
    {
        EdgeWeighted edgeweighted1(1,2);
        EdgeWeighted edgeweighted2(2,1);
        assert((edgeweighted1<edgeweighted2)==0);
        EdgeWeighted edgeweighted3(3,0);
        assert((edgeweighted1<edgeweighted3)==0);
        EdgeWeighted edgeweighted4(1,5);
        assert(edgeweighted1<edgeweighted4);
        EdgeWeighted edgeweighted5(2,2);
        assert(edgeweighted4<edgeweighted5);
    }

    cout << "Testing: > " << endl;
    {
       
        EdgeWeighted edgeweighted1(1,2);
        EdgeWeighted edgeweighted2(2,1);
        assert((edgeweighted2>edgeweighted1)==0);
        EdgeWeighted edgeweighted3(3,0);
        assert((edgeweighted3>edgeweighted1)==0);
    }

    cout << "Testing: EdgeWeighted in set" << endl;
    {
        set<EdgeWeighted> e_set;
        EdgeWeighted edgeweighted(23,43);
        e_set.insert(edgeweighted);
    }

    cout << "Testing: EdgeWeighted in unordered_map" << endl;
    {

        unordered_map<int,EdgeWeighted> e_map;
        EdgeWeighted edgeweighted(1,32);
        e_map[2] = edgeweighted;
    }
  return 0;
}
