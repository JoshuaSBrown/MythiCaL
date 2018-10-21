
#include <iostream>
#include <string>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "../../include/ugly/graph.hpp"
#include "../../include/ugly/graph_node.hpp"
#include "../../include/ugly/graph_algorithms.hpp"
#include "../../include/ugly/edge_weighted.hpp"
#include "../../include/ugly/edge_directed_weighted.hpp"

using namespace ugly;
using namespace std;
using namespace graphalgorithms;

int main(void){

  cout << "Testing: dijkstraGoingFrom" << endl;
  cout << "Test: with weighted edges 1" << endl;
  {
    // Let's begin by first creating a pentagon
    //
    //      -1c-
    //     |    |
    //    5c    2c
    //     |    |
    //    4c -- 3c
    //     
    GraphNode<string> GN0("C");
    GraphNode<string> GN1("C");
    GraphNode<string> GN2("C");
    GraphNode<string> GN3("C");
    GraphNode<string> GN4("C");

    // Edge 1-2
    shared_ptr<Edge> ed1( new EdgeWeighted(1,2,1.0));
    // Edge 2-3
    shared_ptr<Edge> ed2( new EdgeWeighted(2,3,1.0));
    // Edge 3-4
    shared_ptr<Edge> ed3( new EdgeWeighted(3,4,2.0));
    // Edge 4-5
    shared_ptr<Edge> ed4( new EdgeWeighted(4,5,1.0));
    // Edge 5-1
    shared_ptr<Edge> ed5( new EdgeWeighted(5,1,1.0));

    list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

    map<int,GraphNode<string>> nds;
    nds[1]= GN0;
    nds[2]= GN1;
    nds[3]= GN2;
    nds[4]= GN3;
    nds[5]= GN4;

    Graph<string> gc(eds,nds);

    double shortest_distance = dijkstraGoingFrom<string>(2,4,gc);
    assert(static_cast<int>(round(shortest_distance))==3);

    shortest_distance = dijkstraGoingFrom<string>(3,4,gc);
    assert(static_cast<int>(round(shortest_distance))==2);
  }

  cout << "Test: with weighted edges 2" << endl;
  {
    // Let's begin by first creating a pentagon
    //
    //      -1c-     7c - 8c
    //     |    |     |   |
    //    5c    2c - 9c - 6c
    //     |    |
    //    4c -- 3c
    //     
    GraphNode<string> GN0("C");
    GraphNode<string> GN1("C");
    GraphNode<string> GN2("C");
    GraphNode<string> GN3("C");
    GraphNode<string> GN4("C");

    // Edge 1-2
    shared_ptr<Edge> ed1( new EdgeWeighted(1,2,1.0));
    // Edge 2-3
    shared_ptr<Edge> ed2( new EdgeWeighted(2,3,1.0));
    // Edge 3-4
    shared_ptr<Edge> ed3( new EdgeWeighted(3,4,1.0));
    // Edge 4-5
    shared_ptr<Edge> ed4( new EdgeWeighted(4,5,1.0));
    // Edge 5-1
    shared_ptr<Edge> ed5( new EdgeWeighted(5,1,1.0));
    // Edge 9-2
    shared_ptr<Edge> ed6( new EdgeWeighted(9,2,1.0));
    // Edge 9-7
    shared_ptr<Edge> ed7( new EdgeWeighted(9,7,1.0));
    // Edge 7-8
    shared_ptr<Edge> ed8( new EdgeWeighted(7,8,1.0));
    // Edge 8-6
    shared_ptr<Edge> ed9( new EdgeWeighted(8,6,1.0));
    // Edge 9-6
    shared_ptr<Edge> ed10( new EdgeWeighted(9,6,1.0));

    list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9, ed10 };        

    map<int,GraphNode<string>> nds;
    nds[1]= GN0;
    nds[2]= GN1;
    nds[3]= GN2;
    nds[4]= GN3;
    nds[5]= GN4;

    Graph<string> gc(eds,nds);

    double shortest_distance = dijkstraGoingFrom<string>(4,8,gc);
    assert(static_cast<int>(round(shortest_distance))==5);

    shortest_distance = dijkstraGoingFrom<string>(1,6,gc);
    assert(static_cast<int>(round(shortest_distance))==3);
  }

  cout << "Test: with directed and weighted edges 3" << endl;
  {
    // Let's begin by first creating a pentagon
    //
    //      -1c-
    //     |    |
    //    5c    2c
    //     |    |
    //    4c -- 3c
    //     
    GraphNode<string> GN0("C");
    GraphNode<string> GN1("C");
    GraphNode<string> GN2("C");
    GraphNode<string> GN3("C");
    GraphNode<string> GN4("C");

    // Edge 1->2
    shared_ptr<Edge> ed1( new EdgeDirectedWeighted(1,2,1.0));
    // Edge 2->3
    shared_ptr<Edge> ed2( new EdgeDirectedWeighted(2,3,1.0));
    // Edge 3->4
    shared_ptr<Edge> ed3( new EdgeDirectedWeighted(3,4,2.0));
    // Edge 4->5
    shared_ptr<Edge> ed4( new EdgeDirectedWeighted(4,5,1.0));
    // Edge 5->1
    shared_ptr<Edge> ed5( new EdgeDirectedWeighted(5,1,1.0));

    list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

    map<int,GraphNode<string>> nds;
    nds[1]= GN0;
    nds[2]= GN1;
    nds[3]= GN2;
    nds[4]= GN3;
    nds[5]= GN4;

    Graph<string> gc(eds,nds);

    double shortest_distance = dijkstraGoingFrom<string>(2,3,gc);
    assert(static_cast<int>(round(shortest_distance))==1);

    shortest_distance = dijkstraGoingFrom<string>(3,2,gc);
    assert(static_cast<int>(round(shortest_distance))==5);
  }

  cout << "Test: with directed and weighted edges 4" << endl;
  {
    // Let's begin by first creating a pentagon
    //
    //      -1c-
    //     |    |
    //    5c    2c
    //     |    |
    //    4c -- 3c
    //     
    GraphNode<string> GN0("C");
    GraphNode<string> GN1("C");
    GraphNode<string> GN2("C");
    GraphNode<string> GN3("C");
    GraphNode<string> GN4("C");

    // Switched the order so it should not be possible to find a shortest
    // distance between the vertices
    // Edge 2->1
    shared_ptr<Edge> ed1( new EdgeDirectedWeighted(2,1,1.0));
    // Edge 2->3
    shared_ptr<Edge> ed2( new EdgeDirectedWeighted(2,3,1.0));
    // Edge 3->4
    shared_ptr<Edge> ed3( new EdgeDirectedWeighted(3,4,2.0));
    // Edge 4->5
    shared_ptr<Edge> ed4( new EdgeDirectedWeighted(4,5,1.0));
    // Edge 5->1
    shared_ptr<Edge> ed5( new EdgeDirectedWeighted(5,1,1.0));

    list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

    map<int,GraphNode<string>> nds;
    nds[1]= GN0;
    nds[2]= GN1;
    nds[3]= GN2;
    nds[4]= GN3;
    nds[5]= GN4;

    Graph<string> gc(eds,nds);

    double shortest_distance = dijkstraGoingFrom<string>(2,3,gc);
    assert(static_cast<int>(round(shortest_distance))==1);

    bool throwError=false;
    try{
      shortest_distance = dijkstraGoingFrom<string>(3,2,gc);
    }catch(...){
      throwError = true;
    }
    assert(throwError);
  }

  cout << "Testing: maxMinimumDistanceBetweenEveryVertex" << endl;
  {

    // Let's begin by first creating a pentagon
    //
    //      -1c-
    //     |    |
    //    5c    2c
    //     |    |
    //    4c -- 3c
    //     
    GraphNode<string> GN0("C");
    GraphNode<string> GN1("C");
    GraphNode<string> GN2("C");
    GraphNode<string> GN3("C");
    GraphNode<string> GN4("C");

    // Switched the order so it should not be possible to find a shortest
    // distance between the vertices
    // Edge 1->2
    shared_ptr<Edge> ed1( new EdgeDirectedWeighted(1,2,2.0));
    // Edge 2->3
    shared_ptr<Edge> ed2( new EdgeDirectedWeighted(2,3,1.0));
    // Edge 3->4
    shared_ptr<Edge> ed3( new EdgeDirectedWeighted(3,4,1.0));
    // Edge 4->5
    shared_ptr<Edge> ed4( new EdgeDirectedWeighted(4,5,1.0));
    // Edge 5->1
    shared_ptr<Edge> ed5( new EdgeDirectedWeighted(5,1,1.0));

    list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

    map<int,GraphNode<string>> nds;
    nds[1]= GN0;
    nds[2]= GN1;
    nds[3]= GN2;
    nds[4]= GN3;
    nds[5]= GN4;

    Graph<string> gc(eds,nds);

    auto maxMinDistances = maxMinimumDistanceBetweenEveryVertex<string>(gc);
    assert(maxMinDistances.size()==5);
    
    // Max distance starting at vertex 
    // 1 should be a distance of 5
    // 2 should be a distance of 4
    // 3 should be a distance of 5
    // 4 should be a distance of 5
    // 5 should be a distance of 5

    int count4s = 0;
    int count5s = 0;
    for(auto map_iterator : maxMinDistances){
      if(static_cast<int>(round(map_iterator.second))==4) ++count4s;
      if(static_cast<int>(round(map_iterator.second))==5) ++count5s;
    }

    assert(count4s==1);
    assert(count5s==4);
  }

  return 0;
}
