
#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <set>
#include <unordered_map>

#include "../../include/ugly/edge_directed_weighted.hpp"

using namespace ugly;
using namespace std;

int main(){

    cout << "Beginning Test" << endl;
    cout << "Testing: Constructor" << endl;
    { 
        EdgeDirectedWeighted edgedirectedweighted(2,3);
    }

    cout << "Testing: getClassType" << endl;
    {
        assert(EdgeDirectedWeighted::getClassType()==constants::EdgeType::directed_weighted);
    }

    cout << "Testing: getEdgeType" << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted;
        assert(edgedirectedweighted.getEdgeType()==constants::EdgeType::directed_weighted);
    }

    cout << "Testing: directional" << endl;
    {
      EdgeDirectedWeighted EdgeDirectedWeighted;
      assert(EdgeDirectedWeighted.directional());
    }

    cout << "Testing: operator << " << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted(2,3);
        cout << edgedirectedweighted << endl; 
    }

    cout << "Testing: ==" << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted(2,3);
        EdgeDirectedWeighted edgedirectedweighted2(2,4);
        EdgeDirectedWeighted edgedirectedweighted3(3,2);
        assert(edgedirectedweighted==edgedirectedweighted);
        assert((edgedirectedweighted==edgedirectedweighted2)==0);
        assert((edgedirectedweighted==edgedirectedweighted3)==0);
    }

    cout << "Testing: !=" << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted(2,3);
        EdgeDirectedWeighted edgedirectedweighted2(2,4);
        assert(edgedirectedweighted!=edgedirectedweighted2);
        assert((edgedirectedweighted!=edgedirectedweighted)==0);
    }

    cout << "Testing: getVertex1 and getVertex2" << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted(2,3);
        assert(edgedirectedweighted.getVertex1()==2);
        assert(edgedirectedweighted.getVertex2()==3);
    }

    cout << "Testing: < " << endl;
    {
        EdgeDirectedWeighted edgedirectedweighted1(1,2);
        EdgeDirectedWeighted edgedirectedweighted2(2,1);
        assert(edgedirectedweighted1<edgedirectedweighted2);
        EdgeDirectedWeighted edgedirectedweighted3(3,0);
        assert(edgedirectedweighted1<edgedirectedweighted3);
        EdgeDirectedWeighted edgedirectedweighted4(1,5);
        assert(edgedirectedweighted1<edgedirectedweighted4);
        EdgeDirectedWeighted edgedirectedweighted5(2,2);
        assert(edgedirectedweighted4<edgedirectedweighted5);
    }

    cout << "Testing: > " << endl;
    {
       
        EdgeDirectedWeighted edgedirectedweighted1(1,2);
        EdgeDirectedWeighted edgedirectedweighted2(2,1);
        assert(edgedirectedweighted2>edgedirectedweighted1);
        EdgeDirectedWeighted edgedirectedweighted3(3,0);
        assert(edgedirectedweighted3>edgedirectedweighted1);
    }

    cout << "Testing: EdgeDirectedWeighted in set" << endl;
    {
        set<EdgeDirectedWeighted> e_set;
        EdgeDirectedWeighted edgedirectedweighted(23,43);
        e_set.insert(edgedirectedweighted);
    }

    cout << "Testing: EdgeDirectedWeighted in unordered_map" << endl;
    {

        unordered_map<int,EdgeDirectedWeighted> e_map;
        EdgeDirectedWeighted edgedirectedweighted(1,32);
        e_map[2] = edgedirectedweighted;
    }
  return 0;
}
