#include <iostream>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <vector>

#define MAX_ALLOW_GEN 1000
#define CHROMO_LENGTH 5
#define POOL_SIZE 100
#define MAX_NEIGHBOURS 4 //Set equal to 2*number of parents in edge recombination operator, default = 4
#define RECOM_CHANCE 0.2
#define RANDOM_FLOAT static_cast<float>(rand())/static_cast<float>(RAND_MAX)

using namespace std;

class Node{
	private:
		int x = 0, y = 0;
	public:
		Node(int param_x, int param_y){
			x = param_x;
			y = param_y;
			//DEBUG : Displays node coordinates on creation
			//cout << "Made node : " << x << " " << y << endl;
		}
		
		friend float const calcDistance(Node &n1, Node &n2);
};

class Chromosome{
	private:
		int path[CHROMO_LENGTH];
		
	public:
		Chromosome(){
			bool nodeUsed[CHROMO_LENGTH];
			
			for (int i=0; i<CHROMO_LENGTH; i++){
				nodeUsed[i] = false;
			}
			
			for (int i = 0; i<CHROMO_LENGTH; i++){
				int randInt;
				do{
					randInt = rand()%CHROMO_LENGTH;
				}while (nodeUsed[randInt]);
				nodeUsed[randInt] = true;
				path[i] = randInt;
			}
		}
		
		Chromosome(const Chromosome &origin){
			for (int i = 0; i<CHROMO_LENGTH; i++){
				path[i] = origin.path[i];
			}
		}
		
		Chromosome(int givenPath[]){
			for (int i = 0; i<CHROMO_LENGTH; i++){
				path[i] = givenPath[i];
			}
		}
		
		float evaluate(vector<Node> &nodeList){
			float pathLength = 0.0;
			for (int i=1; i<CHROMO_LENGTH; i++){
				pathLength += calcDistance(nodeList[path[i]], nodeList[path[i-1]]);
			}
			pathLength += calcDistance(nodeList[path[0]], nodeList[path[CHROMO_LENGTH-1]]);
			
			return pathLength;
		}
		
		friend ostream &operator<<(ostream& os, const Chromosome &obj){
			for (int i = 0; i<CHROMO_LENGTH; i++){
				cout << setw(2) << obj.path[i] << " ";
			}
			cout << endl;
			return os;
		}
		
		friend Chromosome* edge_recom(Chromosome&, Chromosome&);
		
		friend void updateAdjMatrix(Chromosome&, int[][MAX_NEIGHBOURS]);
};

//FX to generate Node Map
vector<Node> genNodes(int numNodes){

	vector<Node> temp;
	vector<Node>::iterator it = temp.begin();

	for (int i = 0; i<numNodes; i++){
		it = temp.insert(it, Node(rand()%100,rand()%100));
	}

	return temp;
}

//FX to calculate Distace between 2 nodes
float const calcDistance(Node &n1, Node &n2){
	return sqrt(static_cast<float>(pow((n1.x-n2.x),2)) + static_cast<float>(pow((n1.y-n2.y),2)));
}

//FX to create Adjacency Matrix, given 2 chromosomes
void updateAdjMatrix(Chromosome &cr, int adj_matrix[][MAX_NEIGHBOURS]){
	for (int targetGene=0; targetGene<CHROMO_LENGTH; targetGene++){
		for (int i=0; i<MAX_NEIGHBOURS; i++){
			if (adj_matrix[cr.path[targetGene]][i] == cr.path[targetGene==0?CHROMO_LENGTH-1:targetGene-1]){
				break;
			}
			else if ( adj_matrix[cr.path[targetGene]][i]==-1){
				adj_matrix[cr.path[targetGene]][i] = cr.path[targetGene==0?CHROMO_LENGTH-1:targetGene-1];
				break;
			}
		}

		for (int i=0; i<MAX_NEIGHBOURS; i++){
			if (adj_matrix[cr.path[targetGene]][i] == cr.path[targetGene==CHROMO_LENGTH-1?0:targetGene+1]){
				break;
			}
			else if (adj_matrix[cr.path[targetGene]][i]==-1){
				adj_matrix[cr.path[targetGene]][i] = cr.path[targetGene==CHROMO_LENGTH-1?0:targetGene+1];
				break;
			}
		}
	}
}

//FX to create 1 offspring from 2 parents
Chromosome* edge_recom(Chromosome &cr1, Chromosome &cr2, int adj_matrix[][MAX_NEIGHBOURS]){

	//DEBUG: Display cr1, cr2 and Adjacency Matrix for cr1
	/*
	cout << cr1 << endl;
	cout << cr2 << endl;
	
	for (int i=0; i<CHROMO_LENGTH; i++){
		cout << "[" << i << "] : ";
		for (int j=0; j<MAX_NEIGHBOURS; j++){
			cout << adj_matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/

	int newPath[CHROMO_LENGTH];
	
	//DEBUG: Test code for passing of Chromosome new Path
	for (int i = 0; i<CHROMO_LENGTH; i++){
		newPath[i] = i;
	}
	
	int neighbour_count[CHROMO_LENGTH];
	int least_neighbours_count = 99;
	int least_neighbours_node = -1;

	for (int i = 0; i<CHROMO_LENGTH; i++){
		neighbour_count[i] = 0;
	}
	
	for (int node = 0; node<CHROMO_LENGTH; node++){
		for (int neighbour = 0; neighbour<MAX_NEIGHBOURS; neighbour++){
			if (adj_matrix[node][neighbour] == -1) continue;
			else neighbour_count[node]++;
		}
		if (neighbour_count[node] < least_neighbours_count){
			least_neighbours_node = node;
		}
	}
	
	
	return new Chromosome(newPath);
}

int main(){
	// Seed RNG
	srand(time(NULL));
	
	// Generate node map
	vector<Node> nodeList = genNodes(15);
	
	// Generate chromosome pool
	vector<Chromosome> chromoPool;

	vector<Chromosome>::iterator poolIt = chromoPool.begin();
	
	for (int i=0; i<POOL_SIZE; i++){
		Chromosome cr;
		poolIt = chromoPool.insert(poolIt, cr);
		// DEBUG: Displays Chromosome path after insertion into pool
		// cout << cr;
	}
	
	for (int iter = 1; iter <= MAX_ALLOW_GEN; iter++){
		//EVO CODE HERE
	}
	
	//LOOP FOR SELECTION HERE

	//Parent Chromosome selection here
	Chromosome cr1, cr2;
	
	int adj_matrix[CHROMO_LENGTH][MAX_NEIGHBOURS];

	for (int i=0; i<CHROMO_LENGTH; i++){
		for (int j=0; j<MAX_NEIGHBOURS; j++){
			adj_matrix[i][j] = -1;
		}
	}
	
	//Edge Recombination Operator
	if (RANDOM_FLOAT < RECOM_CHANCE){
		updateAdjMatrix(cr1, adj_matrix);
		updateAdjMatrix(cr2, adj_matrix);
		
		Chromosome *temp1 = edge_recom(cr1, cr2, adj_matrix);
		Chromosome *temp2 = edge_recom(cr1, cr2, adj_matrix);
		
		cr1 = *temp1;
		cr2 = *temp2;
	}
	
	cout << cr1 << cr2;

	return 0;
}
