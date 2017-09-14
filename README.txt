OPV Coarse Grain
Riley Hadjis


Note use:
	make VAR="-D _E_" 
to enalbe error outut to the error stream. 

Note use:
	make &>Output
to output all warnings and errors to text file Output

___________________Update Log_________________________
9/13
Additions:
	Prob hop function
	Prob hop off function
	started escape time
To do:
	consolidate functions that are actually needed
	refresh the data struct and class go over old code and notes of how to
	fix shit ton of bugs

7/25
Additions:
	Finished the cluster class
	Wrote the dwell time function
	began Pval function
To do:
	test code
	squash bugs
	make clusterId apart of site
	add to cluster function
	functions have the following dependencies

Name:         | Dwell time | Pvals for neighbors | P hop to neigh      | 
Dependencies: | none       | none                | Dwell time and Pval |

Name:         | Prob off cluster    | Convergence     | Escape Time | escape cluster |
Dependencies: | Dwell time and Pval | Pvals and Phops | Phops       |escapeTime     |  
	

7/24
Additions:
	reworked the cluster class and implemented vector use
To do:
	squash bugs (pushback -> push_back, remove cluster::)
	write test code
	start on real functions

6/12
Additions:
	started to fix some issues in the data struct, it still needs work
	wrote print info function to use in testind
To do:
	clean up testing flags in cluster contructer
	clean up constructer, both site and cluster
	find way around siteInCluster capped at 100 or a fixed value linked list?
	write a more fluid test function
	potentialCluster doesnt work
	figure EXACTLY which functions need to be public and which could be private

6/18
Additions:
	fixed the readability of some of the fuctions
	started to tix the potential cluster and cluster or site functions
	cleaned up error handeling a little
	thought and diagramed what the overall hierarchy of the code should look like:
				
				User---package opv class---clusters using vector wrapper--- site struct

To do:
	Implement the above heirarchy, starting with reconrstructin the site struct
	replace all arrays with vectors
	clean up test code

