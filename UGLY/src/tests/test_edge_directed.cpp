
#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <set>
#include <unordered_map>

#include "../../include/ugly/edge_directed.hpp"

using namespace ugly;
using namespace std;

int main(){

    cout << "Beginning Test" << endl;
    cout << "Testing: Constructor" << endl;
    { 
        EdgeDirected edgedirected(2,3);
    }

    cout << "Testing: getClassType" << endl;
    {
        assert(EdgeDirected::getClassType()==constants::EdgeType::directed);
    }

    cout << "Testing: getEdgeType" << endl;
    {
        EdgeDirected edgedirected;
        assert(edgedirected.getEdgeType()==constants::EdgeType::directed);
    }
    
    cout << "Testing: directional" << endl;
    {
        EdgeDirected edgedirected;
        assert(edgedirected.directional());
    }

    cout << "Testing: operator << " << endl;
    {
        EdgeDirected edgedirected(2,3);
        cout << edgedirected << endl; 
    }

    cout << "Testing: ==" << endl;
    {
        EdgeDirected edgedirected(2,3);
        EdgeDirected edgedirected2(2,4);
        EdgeDirected edgedirected3(3,2);
        assert(edgedirected==edgedirected);
        assert(!(edgedirected==edgedirected3));
        assert(!(edgedirected==edgedirected2));
    }

    cout << "Testing: !=" << endl;
    {
        EdgeDirected edgedirected(2,3);
        EdgeDirected edgedirected2(2,4);
        assert(edgedirected!=edgedirected2);
        assert((edgedirected!=edgedirected)==0);
    }

    cout << "Testing: getVertex1 and getVertex2" << endl;
    {
        EdgeDirected edgedirected(2,3);
        assert(edgedirected.getVertex1()==2);
        assert(edgedirected.getVertex2()==3);

        // Conserves the location of the vertices stored in the edge
        EdgeDirected edgedirected2(3,2);
        assert(edgedirected2.getVertex1()==3);
        assert(edgedirected2.getVertex2()==2);

    }

    cout << "Testing: < " << endl;
    {
        EdgeDirected edgedirected1(1,2);
        EdgeDirected edgedirected2(2,1);
        assert(edgedirected1<edgedirected2);
        EdgeDirected edgedirected3(3,0);
        assert(edgedirected1<edgedirected3);
        EdgeDirected edgedirected4(1,5);
        assert(edgedirected1<edgedirected4);
        EdgeDirected edgedirected5(2,2);
        assert(edgedirected4<edgedirected5);
    }

    cout << "Testing: > " << endl;
    {
       
        EdgeDirected edgedirected1(1,2);
        EdgeDirected edgedirected2(2,1);
        assert(!(edgedirected1>edgedirected2));
        EdgeDirected edgedirected3(3,0);
        assert(edgedirected3>edgedirected1);
    }

    cout << "Testing: EdgeDirected in set" << endl;
    {
        set<EdgeDirected> e_set;
        EdgeDirected edgedirected(23,43);
        e_set.insert(edgedirected);
    }

    cout << "Testing: EdgeDirected in unordered_map" << endl;
    {

        unordered_map<int,EdgeDirected> e_map;
        EdgeDirected edgedirected(1,32);
        e_map[2] = edgedirected;
    }
  return 0;
}
