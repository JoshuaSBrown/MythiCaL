#include <map>
#include <iostream>
#include <string>
#include <cassert>

#include "../../include/ugly/graph_node.hpp"

using namespace ugly;
using namespace std;

int main(void){
    
    cout << "Testing: Contstructor" << endl;
    {
        GraphNode<string,double> GN2("One",1.1);
        GraphNode<double,int,double,string> GN3(3.4,1,45,"two");
    }
    
    cout << "Testing; == && !=" << endl;
    {
        GraphNode<string,double> GN2("One",1.1);
        GraphNode<string,double> GN4("One",1.1);
        assert(GN2==GN2);
        assert(GN2==GN4);

        GraphNode<string,double> GN3("Two",-1.1);
        assert(GN2!=GN3);
    }

    cout << "Testing: < && >" << endl;
    {
        GraphNode<string,double,int> GN0("One",1.1,4);
        GraphNode<string,double,int> GN1("Two",1.1,4);
        GraphNode<string,double,int> GN2("One",2.3,4);

        assert(GN0<GN1);
        assert(GN1>GN0);
        assert(GN0<GN2);
        assert(GN2>GN0);
    }

    cout << "Testing: str()"<< endl;
    {
        GraphNode<string,double,int> GN0("One",1.3,4);
        string var = GN0.getLabel();
        cout << var << endl; 
    }

    cout << "Testing: assignment operator = " << endl;
    {
      map<int, GraphNode<string> > test_map;
      GraphNode<string> GN("One");
      GraphNode<string> GN1;
      GN1 = GN;
      test_map[1]=GN;
    }

    return 0;
}
