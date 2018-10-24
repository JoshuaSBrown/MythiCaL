
#include <iostream>
#include <string>
#include <cassert>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "../libugly/graphvisitor/graphvisitor_depth_first.hpp"
#include "../../include/ugly/edge_weighted.hpp"
#include "../../include/ugly/edge_undirected.hpp"
#include "../../include/ugly/edge_directed.hpp"
#include "../libugly/weak_pointer_supplement.hpp"
using namespace ugly;
using namespace std;

int main(void){

    cout << "Testing: Constructor " << endl;
    {
      GraphVisitorDepthFirst graphvisitordepthfirst;
    }

    cout << "Testing: addEdge" << endl;
    {
      shared_ptr<EdgeWeighted> ed1( new EdgeWeighted(1,2));
      shared_ptr<EdgeWeighted> ed2( new EdgeWeighted(2,3));
      shared_ptr<EdgeWeighted> ed3( new EdgeWeighted(2,1));

      GraphVisitorDepthFirst graphvisitordepthfirst;
      
      bool throwError = false;
      try{
        graphvisitordepthfirst.addEdge(ed1);
      }catch(...){
        throwError=true;
      }
      assert(throwError);

      graphvisitordepthfirst.setStartingVertex(1);
      graphvisitordepthfirst.addEdge(ed1);
     
      // Cannot add the same edge twice 
      throwError=false;
      try{
        graphvisitordepthfirst.addEdge(ed1);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

      // Cannot add an edge that is not associated with a vertex that has been 
      // explored, currently the starting vertex is the only one that satisfies
      // that criteria
      throwError=false;
      try{
        graphvisitordepthfirst.addEdge(ed2);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

      // The depth first graph visitor only allows weighted edges 
      shared_ptr<Edge> ed4( new EdgeUndirected(4,1));

      throwError=false;
      try{
        graphvisitordepthfirst.addEdge(ed4);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

    }

    cout << "Testing: addEdges" << endl;
    {

      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };        

      GraphVisitorDepthFirst graphvisitordepthfirst;
      bool throwError=false;
      try{
        graphvisitordepthfirst.addEdges(eds);
      }catch(...){
        throwError=true;
      }
      assert(throwError);

      graphvisitordepthfirst.setStartingVertex(1);

      throwError=false;
      try{
        graphvisitordepthfirst.addEdges(eds);
      }catch(...){
        throwError=true;
      }
      assert(throwError);

     
      shared_ptr<Edge> ed11( new EdgeUndirected(1,2));
      shared_ptr<Edge> ed22( new EdgeDirected(2,3));
      vector<weak_ptr<Edge>> eds2 = { ed11, ed22 };        

      cout << "ed11 directional? " << ed11->directional() << endl;
      GraphVisitorDepthFirst graphvisitordepthfirst2;
      graphvisitordepthfirst2.setStartingVertex(2);
      
      // Only allows weighted edges
      throwError=false;
      try{
        graphvisitordepthfirst2.addEdges(eds2);
      }catch(...){
        throwError=true;
      }
      assert(throwError);

    }

    cout << "Testing: exploreEdge" << endl;
    {
      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };        

      GraphVisitorDepthFirst graphvisitordepthfirst;

      bool throwError = false;
      try{
        graphvisitordepthfirst.exploreEdge(ed1);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      graphvisitordepthfirst.exploreEdge(ed1);

      // Should now throw an error because ed1 has now been explored
      throwError = false;
      try{
        graphvisitordepthfirst.exploreEdge(ed1);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

      graphvisitordepthfirst.exploreEdge(ed2);
    }


    cout << "Testing: exploreEdge" << endl;
    {
      shared_ptr<Edge> ed1( new Edge(1,2));
      shared_ptr<Edge> ed2( new Edge(2,3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };

      GraphVisitor graphvisitor;

      bool throwError = false;
      try{
        graphvisitor.exploreEdge(ed1);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

    }

    cout << "Testing: getUnexploredVertex" << endl;
    {

      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3));
      shared_ptr<Edge> ed3( new EdgeWeighted(2,4));
      shared_ptr<Edge> ed4( new EdgeWeighted(3,4));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };

      GraphVisitorDepthFirst graphvisitordepthfirst;
      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      auto unexploredvertex = graphvisitordepthfirst.getUnexploredVertex(ed1);
      assert(unexploredvertex==1);

      unexploredvertex = graphvisitordepthfirst.getUnexploredVertex(ed3);
      assert(unexploredvertex==4);

      // Neither vertex of edge 4 is listed as explored so it will
      // throw an error
      bool throwError = false;
      try{
        graphvisitordepthfirst.getUnexploredVertex(ed4);
      }catch(...){
        throwError = true;
      }
      assert(throwError);
    }

    cout << "Testing: getExploredVertex" << endl;
    {
      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3));
      shared_ptr<Edge> ed3( new EdgeWeighted(4,3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };

      GraphVisitorDepthFirst graphvisitordepthfirst;
      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      bool throwError = false;
      try{
        graphvisitordepthfirst.getExploredVertex(ed3);
      }catch(...){
        throwError = true;
      }
      assert(throwError);

      auto exploredvertex = graphvisitordepthfirst.getExploredVertex(ed1);
      assert(exploredvertex==2);
    }

    cout << "Testing: allEdgesExplored" << endl;
    {
      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };        
      GraphVisitorDepthFirst graphvisitordepthfirst;
      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      bool complete = graphvisitordepthfirst.allEdgesExplored();
      assert(complete==false);
      graphvisitordepthfirst.exploreEdge(ed1);
      complete = graphvisitordepthfirst.allEdgesExplored();
      assert(complete==false);
      graphvisitordepthfirst.exploreEdge(ed2);
      complete = graphvisitordepthfirst.allEdgesExplored();
      assert(complete==true);
    }

    cout << "Testing: getDistanceOfVertex" << endl;
    {
      shared_ptr<Edge> ed1( new EdgeWeighted(1,2));
      shared_ptr<Edge> ed2( new EdgeWeighted(2,3,32));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };        
      GraphVisitorDepthFirst graphvisitordepthfirst;
      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      graphvisitordepthfirst.exploreEdge(ed1);
      graphvisitordepthfirst.exploreEdge(ed2);

      auto dist = graphvisitordepthfirst.getDistanceOfVertex(2);
      assert(static_cast<int>(round(dist))==0);
      dist = graphvisitordepthfirst.getDistanceOfVertex(1);
      cout << "distance " << dist << endl;
      //assert(static_cast<int>(round(dist))==1);
      dist = graphvisitordepthfirst.getDistanceOfVertex(3);
      cout << "distance " << dist << endl;
      assert(static_cast<int>(round(dist))==32);
    }

    cout << "Testing: getNextEdge" << endl;
    {
      shared_ptr<EdgeWeighted> ed1( new EdgeWeighted(1,2,0.5));
      shared_ptr<EdgeWeighted> ed2( new EdgeWeighted(2,3,2.3));
      vector<weak_ptr<Edge>> eds = { ed1, ed2 };        
      GraphVisitorDepthFirst graphvisitordepthfirst;
      graphvisitordepthfirst.setStartingVertex(2);
      graphvisitordepthfirst.addEdges(eds);

      auto ed = graphvisitordepthfirst.getNextEdge<EdgeWeighted>();
      assert(ed==ed1);

      graphvisitordepthfirst.exploreEdge(ed);
      ed = graphvisitordepthfirst.getNextEdge<EdgeWeighted>();
      assert(ed==ed2);
    }
    return 0;
}
