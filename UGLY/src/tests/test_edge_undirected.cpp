
#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <set>
#include <unordered_map>

#include "../../include/ugly/edge_undirected.hpp"

using namespace ugly;
using namespace std;

int main(){

    cout << "Beginning Test" << endl;
    cout << "Testing: Constructor" << endl;
    { 
        EdgeUndirected edgeundirected(2,3);
    }

    cout << "Testing: getClassType" << endl;
    {
        assert(EdgeUndirected::getClassType()==constants::EdgeType::undirected);
    }

    cout << "Testing: getEdgeType" << endl;
    {
        EdgeUndirected edgeundirected;
        assert(edgeundirected.getEdgeType()==constants::EdgeType::undirected);
    }

    cout << "Testing: directional" << endl;
    {
      EdgeUndirected EdgeUndirected;
      assert(!EdgeUndirected.directional());
    }

    cout << "Testing: operator << " << endl;
    {
        EdgeUndirected edgeundirected(2,3);
        cout << edgeundirected << endl; 
    }

    cout << "Testing: ==" << endl;
    {
        EdgeUndirected edgeundirected(2,3);
        EdgeUndirected edgeundirected2(2,4);
        EdgeUndirected edgeundirected3(3,2);
        assert(edgeundirected==edgeundirected);
        assert(edgeundirected==edgeundirected3);
        assert((edgeundirected==edgeundirected2)==0);
    }

    cout << "Testing: !=" << endl;
    {
        EdgeUndirected edgeundirected(2,3);
        EdgeUndirected edgeundirected2(2,4);
        assert(edgeundirected!=edgeundirected2);
        assert((edgeundirected!=edgeundirected)==0);
    }

    cout << "Testing: getVertex1 and getVertex2" << endl;
    {
        EdgeUndirected edgeundirected(2,3);
        assert(edgeundirected.getVertex1()==2);
        assert(edgeundirected.getVertex2()==3);
        // Does not conserve the location, always places the smaller vertex id
        // at vertex1 
        EdgeUndirected edgeundirected2(2,3);
        assert(edgeundirected2.getVertex1()==2);
        assert(edgeundirected2.getVertex2()==3);
    }

    cout << "Testing: < " << endl;
    {
        EdgeUndirected edgeundirected1(1,2);
        EdgeUndirected edgeundirected2(2,1);
        assert(!(edgeundirected1<edgeundirected2));
        EdgeUndirected edgeundirected3(3,0);
        assert(edgeundirected3<edgeundirected1);
        EdgeUndirected edgeundirected4(1,5);
        assert(edgeundirected1<edgeundirected4);
        EdgeUndirected edgeundirected5(2,2);
        assert(edgeundirected4<edgeundirected5);
    }

    cout << "Testing: > " << endl;
    {
       
        EdgeUndirected edgeundirected1(1,2);
        EdgeUndirected edgeundirected2(2,1);
        assert(!(edgeundirected1>edgeundirected2));
        EdgeUndirected edgeundirected3(3,0);
        assert(!(edgeundirected3>edgeundirected1));
    }

    cout << "Testing: EdgeUndirected in set" << endl;
    {
        set<EdgeUndirected> e_set;
        EdgeUndirected edgeundirected(23,43);
        e_set.insert(edgeundirected);
    }

    cout << "Testing: EdgeUndirected in unordered_map" << endl;
    {

        unordered_map<int,EdgeUndirected> e_map;
        EdgeUndirected edgeundirected(1,32);
        e_map[2] = edgeundirected;
    }
  return 0;
}
