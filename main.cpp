#include <iostream>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <vector>

#define MAX_ALLOW_GEN 1000
#define CHROMO_LENGTH 15
#define POOL_SIZE 100
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
		short int path[CHROMO_LENGTH];
		
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
};

float const calcDistance(Node &n1, Node &n2){
	return sqrt(static_cast<float>(pow((n1.x-n2.x),2)) + static_cast<float>(pow((n1.y-n2.y),2)));
}

vector<Node> genNodes(int numNodes){

	vector<Node> temp;
	vector<Node>::iterator it = temp.begin();

	for (int i = 0; i<numNodes; i++){
		it = temp.insert(it, Node(rand()%100,rand()%100));
	}
	
	return temp;
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

	}

	return 0;
}
