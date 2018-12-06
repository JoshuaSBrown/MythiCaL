
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include "../../libkmccoarsegrain/kmc_basin_explorer.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

  cout << "Testing: findBasin" << endl;
  cout << "Test 1" << endl;
  {

    // 
    // site1 -> site2 -> site3 -> site4
    //       <-       <-       <-
    //
    // Site2 and Site3 should be identified as clusters if exploration is 
    // started on either of them. 


    KMC_Site site1;
    KMC_Site site2;
    KMC_Site site3;
    KMC_Site site4;

    site1.setId(1);
    site2.setId(2);
    site3.setId(3);
    site4.setId(4);

    double fast = 100;
    double slow = 1;

    double rate1_2 = fast;
    double rate2_1 = slow;
    double rate2_3 = fast;
    double rate3_2 = fast;
    double rate3_4 = slow;
    double rate4_3 = fast;

    unordered_map<int, double *> neigh_rates_site1;
    neigh_rates_site1[2] = &rate1_2;
    site1.setRatesToNeighbors(neigh_rates_site1);

    unordered_map<int, double *> neigh_rates_site2;
    neigh_rates_site2[1] = &rate2_1;
    neigh_rates_site2[3] = &rate2_3;
    site2.setRatesToNeighbors(neigh_rates_site2);
     
    unordered_map<int, double *> neigh_rates_site3;
    neigh_rates_site3[2] = &rate3_2;
    neigh_rates_site3[4] = &rate3_4;
    site3.setRatesToNeighbors(neigh_rates_site3);

    unordered_map<int, double *> neigh_rates_site4;
    neigh_rates_site4[3] = &rate4_3;
    site4.setRatesToNeighbors(neigh_rates_site4);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site1);
    site_container.addKMC_Site(site2);
    site_container.addKMC_Site(site3);
    site_container.addKMC_Site(site4);

    // Need at least an empty Cluster container
    KMC_Cluster_Container clusters;

    BasinExplorer basin_explorer; 
    auto vertices = basin_explorer.findBasin(site_container,clusters,2);
    assert(vertices.size()==2);

    bool found2 = false;
    bool found3 = false;
    for(auto siteId : vertices){
      if(siteId==2) found2 = true;
      if(siteId==3) found3 = true;
    }
    assert(found2);
    assert(found3);

  }

  cout << "Test 2" << endl;
  {

    // 
    // site1 -> site2 -> site3 -> site4
    //   |^  <-  |^   <-  |^   <-  ^|
    //   V|      v|       v|       |v
    // site5 -> site6 -> site7 -> site8 
    //   |^  <-  |^   <-  |^   <-  ^|
    //   V|      v|       v|       |v
    // site9 -> site10-> site11-> site12
    //   |^  <-  |^   <-  |^   <-  ^|
    //   V|      v|       v|       |v
    // site13-> site14-> site15-> site16
    //
    // Here we will assume hopping between sites 6, 7 and 11 is very fast

    double fast = 10;
    double slow = 1;

    double rate1_2 = slow*0.94;
    double rate1_5 = slow*0.3;
    unordered_map<int,double *> site_1_rates;
    site_1_rates[2] = &rate1_2;
    site_1_rates[5] = &rate1_5;

    double rate2_1 = slow*0.47;
    double rate2_6 = fast*1.1;
    double rate2_3 = fast*1.4;
    unordered_map<int,double *> site_2_rates;
    site_2_rates[1] = &rate2_1;
    site_2_rates[6] = &rate2_6;
    site_2_rates[3] = &rate2_3;

    double rate3_2 = slow*0.8;
    double rate3_4 = fast*1.3;
    double rate3_7 = fast*0.8;
    unordered_map<int,double *> site_3_rates;
    site_3_rates[2] = &rate3_2;
    site_3_rates[4] = &rate3_4;
    site_3_rates[7] = &rate3_7;

    double rate4_3 = slow*1.0;
    double rate4_8 = fast*1.0;
    unordered_map<int,double *> site_4_rates;
    site_4_rates[3] = &rate4_3;
    site_4_rates[8] = &rate4_8;

    double rate5_1 = slow*1.3;
    double rate5_6 = fast*0.4;
    double rate5_9 = fast*1.1;
    unordered_map<int,double *> site_5_rates;
    site_5_rates[1] = &rate5_1;
    site_5_rates[6] = &rate5_6;
    site_5_rates[9] = &rate5_9;
    
    double rate6_2 = slow*0.01;
    double rate6_5 = slow*0.8;
    double rate6_7 = fast*1.2;
    double rate6_10 = slow*1.1;
    unordered_map<int,double *> site_6_rates;
    site_6_rates[2] = &rate6_2;
    site_6_rates[5] = &rate6_5;
    site_6_rates[7] = &rate6_7;
    site_6_rates[10] = &rate6_10;

    double rate7_3 = slow*1.0;
    double rate7_6 = fast*1.5;
    double rate7_8 = slow*1.2;
    double rate7_11 = fast*0.9;
    unordered_map<int,double *> site_7_rates;
    site_7_rates[3] = &rate7_3;
    site_7_rates[6] = &rate7_6;
    site_7_rates[8] = &rate7_8;
    site_7_rates[11] = &rate7_11;

    double rate8_4 = slow*1.2;
    double rate8_7 = fast*1.9;
    double rate8_12 = slow*1.0;
    unordered_map<int,double *> site_8_rates;
    site_8_rates[4] = &rate8_4;
    site_8_rates[7] = &rate8_7;
    site_8_rates[12] = &rate8_12;

    double rate9_5 = slow*0.4;
    double rate9_10 = fast*1.4;
    double rate9_13 = fast*0.8;
    unordered_map<int,double *> site_9_rates;
    site_9_rates[5] = &rate9_5;
    site_9_rates[10] = &rate9_10;
    site_9_rates[13] = &rate9_13;

    double rate10_6 = fast*2.1;
    double rate10_9 = fast*0.8;
    double rate10_11 = fast*1.2;
    double rate10_14 = slow*1.3;
    unordered_map<int,double *> site_10_rates;
    site_10_rates[6] = &rate10_6;
    site_10_rates[9] = &rate10_9;
    site_10_rates[11] = &rate10_11;
    site_10_rates[14] = &rate10_14;

    double rate11_7 = fast*0.8;
    double rate11_10 = slow*1.0;
    double rate11_12 = slow*1.3;
    double rate11_15 = slow*0.3;
    unordered_map<int,double *> site_11_rates;
    site_11_rates[7] = &rate11_7;
    site_11_rates[10] = &rate11_10;
    site_11_rates[12] = &rate11_12;
    site_11_rates[15] = &rate11_15;

    double rate12_8 = slow*1.2;
    double rate12_11 = fast*1.1;
    double rate12_16 = fast*0.8;
    unordered_map<int,double *> site_12_rates;
    site_12_rates[8] = &rate12_8;
    site_12_rates[11] = &rate12_11;
    site_12_rates[16] = &rate12_16;
    
    double rate13_9 = slow*0.8;
    double rate13_14 = slow*0.3;
    unordered_map<int,double *> site_13_rates;
    site_13_rates[9] = &rate13_9;
    site_13_rates[14] = &rate13_14;

    double rate14_10 = slow*0.75;
    double rate14_13 = fast*1.0;
    double rate14_15 = slow*1.1;
    unordered_map<int,double *> site_14_rates;
    site_14_rates[10] = &rate14_10;
    site_14_rates[13] = &rate14_13;
    site_14_rates[15] = &rate14_15;

    double rate15_11 = fast*1.5;
    double rate15_14 = slow*0.8;
    double rate15_16 = slow*1.0;
    unordered_map<int,double *> site_15_rates;
    site_15_rates[11] = &rate15_11;
    site_15_rates[14] = &rate15_14;
    site_15_rates[16] = &rate15_16;

    double rate16_12 = fast*1.1;
    double rate16_15 = slow*0.9;
    unordered_map<int,double *> site_16_rates;
    site_16_rates[12] = &rate16_12;
    site_16_rates[15] = &rate16_15;

    KMC_Site site1;
    KMC_Site site2;
    KMC_Site site3;
    KMC_Site site4;
    KMC_Site site5;
    KMC_Site site6;
    KMC_Site site7;
    KMC_Site site8;
    KMC_Site site9;
    KMC_Site site10;
    KMC_Site site11;
    KMC_Site site12;
    KMC_Site site13;
    KMC_Site site14;
    KMC_Site site15;
    KMC_Site site16;
    
    site1.setRatesToNeighbors(site_1_rates);
    site2.setRatesToNeighbors(site_2_rates);
    site3.setRatesToNeighbors(site_3_rates);
    site4.setRatesToNeighbors(site_4_rates);
    site5.setRatesToNeighbors(site_5_rates);
    site6.setRatesToNeighbors(site_6_rates);
    site7.setRatesToNeighbors(site_7_rates);
    site8.setRatesToNeighbors(site_8_rates);
    site9.setRatesToNeighbors(site_9_rates);
    site10.setRatesToNeighbors(site_10_rates);
    site11.setRatesToNeighbors(site_11_rates);
    site12.setRatesToNeighbors(site_12_rates);
    site13.setRatesToNeighbors(site_13_rates);
    site14.setRatesToNeighbors(site_14_rates);
    site15.setRatesToNeighbors(site_15_rates);
    site16.setRatesToNeighbors(site_16_rates);

    site1.setId(1);
    site2.setId(2);
    site3.setId(3);
    site4.setId(4);
    site5.setId(5);
    site6.setId(6);
    site7.setId(7);
    site8.setId(8);
    site9.setId(9);
    site10.setId(10);
    site11.setId(11);
    site12.setId(12);
    site13.setId(13);
    site14.setId(14);
    site15.setId(15);
    site16.setId(16);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site1);
    site_container.addKMC_Site(site2);
    site_container.addKMC_Site(site3);
    site_container.addKMC_Site(site4);
    site_container.addKMC_Site(site5);
    site_container.addKMC_Site(site6);
    site_container.addKMC_Site(site7);
    site_container.addKMC_Site(site8);
    site_container.addKMC_Site(site9);
    site_container.addKMC_Site(site10);
    site_container.addKMC_Site(site11);
    site_container.addKMC_Site(site12);
    site_container.addKMC_Site(site13);
    site_container.addKMC_Site(site14);
    site_container.addKMC_Site(site15);
    site_container.addKMC_Site(site16);
  
    // Need at least an empty Cluster container
    KMC_Cluster_Container clusters;

    BasinExplorer basin_explorer; 
    auto vertices = basin_explorer.findBasin(site_container,clusters,6);

    for(auto siteId : vertices){
      cout << "site id " << siteId << endl;
    }

    vertices = basin_explorer.findBasin(site_container,clusters,7);

    for(auto siteId : vertices){
      cout << "site id " << siteId << endl;
    }

    vertices = basin_explorer.findBasin(site_container,clusters,11);

    for(auto siteId : vertices){
      cout << "site id " << siteId << endl;
    }

  }

  return 0;
}
