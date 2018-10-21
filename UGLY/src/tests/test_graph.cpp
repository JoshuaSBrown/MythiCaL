
#include <iostream>
#include <string>
#include <cassert>
#include <list>
#include <vector>
#include <map>

#include "../../include/ugly/graph.hpp"
#include "../../include/ugly/graph_node.hpp"
#include "../edge/edge.hpp"

using namespace ugly;
using namespace std;

int main(void){

    cout << "Testing: Constructor " << endl;
    {
        GraphNode<string> GN0("B");
        GraphNode<string> GN1("C");
        GraphNode<string> GN2("A");

        shared_ptr<Edge> ed1( new Edge(1,2));
        shared_ptr<Edge> ed2( new Edge(2,3));

        list<weak_ptr<Edge>> eds = { ed1, ed2 };        
    
        map<int,GraphNode<string>> nds;
        nds[1]=GN0;
        nds[2]=GN1;
        nds[3]=GN2;
        Graph<string> gc(eds,nds);
        assert(gc.getLabel()=="A,B,C,");
        
    }

    cout << "Testing: directional " << endl;
    {
        GraphNode<string> GN0("B");
        GraphNode<string> GN1("C");
        GraphNode<string> GN2("A");

        shared_ptr<Edge> ed1( new Edge(1,2));
        shared_ptr<Edge> ed2( new Edge(2,3));

        list<weak_ptr<Edge>> eds = { ed1, ed2 };        
    
        map<int,GraphNode<string>> nds;
        nds[1]=GN0;
        nds[2]=GN1;
        nds[3]=GN2;
        Graph<string> gc(eds,nds);
        // Edges are by default directional so this should evaluate to true
        assert(gc.directional()); 
    }

    cout << "Testing: == & != " << endl;
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

        shared_ptr<Edge> ed1( new Edge(1,2));
        shared_ptr<Edge> ed2( new Edge(2,3));
        shared_ptr<Edge> ed3( new Edge(3,4));
        shared_ptr<Edge> ed4( new Edge(4,5));
        shared_ptr<Edge> ed5( new Edge(5,1));

        list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

        map<int,GraphNode<string>> nds;
        nds[1]= GN0;
        nds[2]= GN1;
        nds[3]= GN2;
        nds[4]= GN3;
        nds[5]= GN4;

        Graph<string> gc(eds,nds);
        assert( gc.getLabel()=="C,C,C,C,C,");

        Graph<string> gc2(eds,nds);

        assert(gc==gc2);
        //
        //      -1c-
        //     |    |
        //    5H    2c
        //     |    |
        //    4c -- 3c
        //     
        GraphNode<string> GN5("H");
        map<int,GraphNode<string>> nds2;
        nds2[1]= GN0;
        nds2[2]= GN1;
        nds2[3]= GN2;
        nds2[4]= GN3;
        nds2[5]= GN5;


        Graph<string> gc3(eds,nds2);
        assert( gc3.getLabel()=="C,C,C,C,H,");

        assert(gc3!=gc2);
    }

    cout << "Testing: getVertices" << endl;
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

      shared_ptr<Edge> ed1( new Edge(1,2));
      shared_ptr<Edge> ed2( new Edge(2,3));
      shared_ptr<Edge> ed3( new Edge(3,4));
      shared_ptr<Edge> ed4( new Edge(4,5));
      shared_ptr<Edge> ed5( new Edge(5,1));

      list<weak_ptr<Edge>> eds = { ed1, ed2, ed3, ed4, ed5 };        

      map<int,GraphNode<string>> nds;
      nds[1]= GN0;
      nds[2]= GN1;
      nds[3]= GN2;
      nds[4]= GN3;
      nds[5]= GN4;

      Graph<string> gc(eds,nds);

      auto vertices = gc.getVertices();

      assert(vertices.size()==5);
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      bool found5 = false;
      for(auto vertex : vertices ){
        if(vertex==1) found1=true;
        if(vertex==2) found2=true;
        if(vertex==3) found3=true;
        if(vertex==4) found4=true;
        if(vertex==5) found5=true;
      }
      assert(found1);
      assert(found2);
      assert(found3);
      assert(found4);
      assert(found5);
    }
    return 0;
}
