#include <iostream>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <vector>
#include <chrono>

#include <windows.h>    // Win32API Header File
#include <cstring>
#include <cstdio>
#include <unistd.h>

int iter = 0;

#define CHROMO_LENGTH 50
#define MAX_ALLOW_GEN CHROMO_LENGTH*5000
#define POOL_SIZE 100
#define MAX_NEIGHBOURS 4 //Set equal to 2*number of parents in edge recombination operator, default = 4
#define RECOM_CHANCE 0.2
#define MUT_CHANCE 0.05
#define UPDATE_INTERVAL 20 //The number of generations evaluated before a screen redraw of the best chromosomes, lower will be more frequent updates, at cost of performance
#define REINSERT_CHANCE 0.005
#define REINSERT_GEN static_cast<int>(min(1000,(MAX_ALLOW_GEN/50))) //The number of generations required to pass before algorithm tries to reinsert best chromosome into the pool
#define SELECT_PRESSURE (iter<MAX_ALLOW_GEN/2)?0.5:1.5 //Default to 1, Higher numbers bias selection towards better fitness
#define DIFF_FRACTION_BEST 0.25 //The threshold for random reinsertion of best chromosome, triggers when difference of gen best with world best > fraction of best
#define RANDOM_FLOAT static_cast<float>(rand())/static_cast<float>(RAND_MAX)
#define HARD_THRESHOLD MAX_ALLOW_GEN/10 // Forces exit if HARD_THRESHOLD generations have passed without improvement of global best. set to MAX_ALLOW_GEN for no hard_threshold

#define Red  RGB (255,0,0)
#define Lime RGB (206,255,0)
#define Blue RGB (0,0,255)

#define XScalar 5
#define YScalar 5
#define Xoffset1 100
#define Xoffset2 Xoffset1+600
#define Yoffset 225
#define CircleRadius static_cast<int>(XScalar+1/2)

using namespace std;

HWND    GetConsoleWndHandle (void);
static HWND    hConWnd = GetConsoleWndHandle();
int     BCX_Line (HWND,int,int,int,int,int=0,HDC=0);
int     BCX_Circle (HWND,int,int,int,int=0,int=0,HDC=0);


class Node{
	private:
		
	public:
		int x = 0, y = 0;
		Node(int param_x, int param_y){
			x = param_x;
			y = param_y;
			//DEBUG : Displays node coordinates on creation
			//cout << "Made node : " << x << " " << y << endl;
		}
		
		friend float const calcDistance(Node &n1, Node &n2);
		
		friend ostream &operator<<(ostream &os, const vector<Node> &obj){
			for (int i = 0; i<CHROMO_LENGTH; i++){
				cout << "[" << setw(2) << i << "] = (" << setw(2) << obj[i].x << "," << setw(2) << obj[i].y << ")" << endl;
			}
			
			return os;
		}
};

class Chromosome{
	private:
		
		
	public:
		int path[CHROMO_LENGTH];
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
		
		friend Chromosome* edge_recom(Chromosome&, Chromosome&,int adj_matrix[][MAX_NEIGHBOURS]);
		
		friend void updateAdjMatrix(Chromosome&, int[][MAX_NEIGHBOURS]);
		
		void mutate(){
			int rand1 = rand()%CHROMO_LENGTH;
			int rand2 = rand()%CHROMO_LENGTH;
			
			int temp = path[rand1];
			path[rand1] = path[rand2];
			path[rand2] = temp;
		}
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

//FX to calculate Distance between 2 nodes
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

	if (RANDOM_FLOAT<0.5) newPath[0] = cr1.path[0];
	else newPath[0] = cr2.path[0];

	int addedNode = newPath[0];

	for (int newPathLength = 1; newPathLength<CHROMO_LENGTH; newPathLength++){
		
		int neighbour_count[CHROMO_LENGTH];
		int least_neighbours_count = 99;
		vector<int> least_neighbours_node;
	
		//Reset neighbour count
		for (int i = 0; i<CHROMO_LENGTH; i++){
			neighbour_count[i] = 0;
		}
		
		//Void out added node
		for (int i = 0; i<MAX_NEIGHBOURS; i++){
			adj_matrix[addedNode][i] = 99;
		}

		//Find out least node
		for (int node = 0; node<CHROMO_LENGTH; node++){
			for (int neighbour = 0; neighbour<MAX_NEIGHBOURS; neighbour++){
				if (adj_matrix[node][neighbour] == addedNode)
					adj_matrix[node][neighbour] = -1;
				
				if (adj_matrix[node][neighbour] == -1)
					continue;
				else
					neighbour_count[node]++;
			}
			if (neighbour_count[node] < least_neighbours_count){
				least_neighbours_count = neighbour_count[node];
				least_neighbours_node.clear();
				least_neighbours_node.push_back(node);
			}
			else if (neighbour_count[node] == least_neighbours_count){
				least_neighbours_node.push_back(node);
			}
		}
		
		float nodeProbSelect = static_cast<float>(1.0/least_neighbours_node.size());
		float randomNodeProb = RANDOM_FLOAT;
		
		int i = 1;
		
		while(randomNodeProb > i*nodeProbSelect){
			i++;
		}
		newPath[newPathLength] = least_neighbours_node[i-1];
		
		addedNode = newPath[newPathLength];
	}
	
	return new Chromosome(newPath);
}

std::ostream&
display(std::ostream& os, std::chrono::nanoseconds ns)
{
    using namespace std;
    using namespace std::chrono;
    typedef duration<int, ratio<86400>> days;
    char fill = os.fill();
    os.fill('0');
    auto d = duration_cast<days>(ns);
    ns -= d;
    auto h = duration_cast<hours>(ns);
    ns -= h;
    auto m = duration_cast<minutes>(ns);
    ns -= m;
    auto s = duration_cast<seconds>(ns);
    os << setw(2) << d.count() << "d:"
       << setw(2) << h.count() << "h:"
       << setw(2) << m.count() << "m:"
       << setw(2) << s.count() << 's';
    os.fill(fill);
    return os;
};

int BCX_Line (HWND Wnd,int x1,int y1,int x2,int y2,int Pen,HDC DrawHDC)
{
  int a,b=0;
  HPEN hOPen;
  // penstyle, width, color
  HPEN hNPen = CreatePen(PS_SOLID, 2, Pen);
  if (!DrawHDC) DrawHDC = GetDC(Wnd), b = 1;
  hOPen = (HPEN)SelectObject(DrawHDC, hNPen);
  // starting point of line
  MoveToEx(DrawHDC, x1, y1, NULL);
  // ending point of line
  a = LineTo(DrawHDC, x2, y2);
  DeleteObject(SelectObject(DrawHDC, hOPen));
  if (b) ReleaseDC(Wnd, DrawHDC);
  return a;
}
// converts circle(centerX,centerY,radius,pen) to WinApi function
// ellipse inside box with upper left and lower right coordinates
int BCX_Circle(HWND Wnd,int X,int Y,int R,int Pen,int Fill,HDC DrawHDC)
{
  int a, b = 0;
  if (!DrawHDC) DrawHDC = GetDC(Wnd), b = 1;
  // penstyle, width, color
  HPEN   hNPen = CreatePen(PS_SOLID, 2, Pen);
  HPEN   hOPen = (HPEN)SelectObject(DrawHDC, hNPen);
  HBRUSH hOldBrush;
  HBRUSH hNewBrush;
  // if true will fill circle with pencolor
  if (Fill)
  {
    hNewBrush = CreateSolidBrush(Pen);
    hOldBrush = (HBRUSH)SelectObject(DrawHDC, hNewBrush);
  }
  else
  {
    hNewBrush = (HBRUSH)GetStockObject(NULL_BRUSH);
    hOldBrush = (HBRUSH)SelectObject(DrawHDC, hNewBrush);
  }
  a = Ellipse(DrawHDC, X-R, Y+R, X+R, Y-R);
  DeleteObject(SelectObject(DrawHDC, hOPen));
  DeleteObject(SelectObject(DrawHDC, hOldBrush));
  if (b) ReleaseDC(Wnd, DrawHDC);
  return a;
}
// the hoop ...
HWND GetConsoleWndHandle(void)
{
  HWND hConWnd;
  OSVERSIONINFO os;
  char szTempTitle[64], szClassName[128], szOriginalTitle[1024];
  os.dwOSVersionInfoSize = sizeof( OSVERSIONINFO );
  GetVersionEx( &os );
  // may not work on WIN9x
  if ( os.dwPlatformId == VER_PLATFORM_WIN32s ) return 0;
  GetConsoleTitle( szOriginalTitle, sizeof( szOriginalTitle ) );
  sprintf( szTempTitle,"%u - %u", GetTickCount(), GetCurrentProcessId() );
  SetConsoleTitle( szTempTitle );
  Sleep( 40 );
  // handle for NT
  hConWnd = FindWindow( NULL, szTempTitle );
  SetConsoleTitle( szOriginalTitle );
  // may not work on WIN9x
  if ( os.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS )
  {
    hConWnd = GetWindow( hConWnd, GW_CHILD );
    if ( hConWnd == NULL )  return 0;
    GetClassName( hConWnd, szClassName, sizeof ( szClassName ) );
    while ( strcmp( szClassName, "ttyGrab" ) != 0 )
    {
      hConWnd = GetNextWindow( hConWnd, GW_HWNDNEXT );
      if ( hConWnd == NULL )  return 0;
      GetClassName( hConWnd, szClassName, sizeof( szClassName ) );
    }
  }
  return hConWnd;
}

int main(){
	
	std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    // Seed RNG
	srand(time(NULL));
	
	// Generate node map
	vector<Node> nodeList = genNodes(CHROMO_LENGTH);
	
	// Generate chromosome pool
	vector<Chromosome> chromoPool;

	for (int i=0; i<POOL_SIZE; i++){
		Chromosome cr;
		chromoPool.push_back(cr);
		// DEBUG: Displays Chromosome path after insertion into pool
		// cout << cr;
	}
	
	Chromosome bestEver;
	Chromosome bestThisGen;
	int bestEverAtGen = 0;
	vector<int> bestTracker;
	bestTracker.reserve(MAX_ALLOW_GEN/50);
	
	bestTracker.push_back(0);
	
	cout << MAX_ALLOW_GEN << endl; cin.get();
	
	for (iter = 1; (iter <= MAX_ALLOW_GEN); iter++){
		if ((iter-bestTracker.back()) > HARD_THRESHOLD) break;
		//DEBUG Display Chromosome Pool
		/*
		cout << " ==== POOL " << iter << " ==== " << endl;
		for (int i = 0; i<POOL_SIZE; i++){
			cout << chromoPool[i] << "\t" << chromoPool[i].evaluate(nodeList) << endl;
		}
		cin.get();
		*/

		bestThisGen = chromoPool[0];
		
		//EVAL CODE
		float fitness[POOL_SIZE];
		float cumulativeProb[POOL_SIZE];
		float curCumProb = 0.0;
		float curTotalFitness = 0.0;
		float curWorstFitness = 0.0;
		int trackWorst = -1;
		
		// Eval all chromosomes, find worst fitness
		for (int chromoNum = 0; chromoNum<POOL_SIZE; chromoNum++){
			fitness[chromoNum] = chromoPool[chromoNum].evaluate(nodeList);
			
			if (fitness[chromoNum] > curWorstFitness){
				curWorstFitness = fitness[chromoNum];
				trackWorst = chromoNum;
				//DEBUG : Show worst ftiness this gen
				//cout << " Worst" <<  curWorstFitness << endl;
			}
			
			if (fitness[chromoNum] < bestEver.evaluate(nodeList)){
				bestEver = Chromosome(chromoPool[chromoNum]);
				bestEverAtGen = iter;
				//DEBUG : For tracking best chromosome updates occuring at what gens
				bestTracker.push_back(iter);
			}
			
			if (fitness[chromoNum] < bestThisGen.evaluate(nodeList)){
				bestThisGen = Chromosome(chromoPool[chromoNum]);
			}
		}
		
		// Randomly reinsert best ever chromosome into the pool, replacing the worst this generation
		// If iterations since last improvment > REINSERT_GEN
		if ((bestThisGen.evaluate(nodeList)!=bestEver.evaluate(nodeList))&&(((iter - bestTracker.back()) > REINSERT_GEN)||((bestThisGen.evaluate(nodeList) - bestEver.evaluate(nodeList)) >= (DIFF_FRACTION_BEST*bestEver.evaluate(nodeList)))) && RANDOM_FLOAT < REINSERT_CHANCE){
			chromoPool.erase(chromoPool.begin() + trackWorst);
			chromoPool.insert(chromoPool.begin() + trackWorst, bestEver);
			
			fitness[trackWorst] = bestEver.evaluate(nodeList);
		}
		
		// Adjust all fitness, sum total adjusted fitness
		for (int chromoNum = 0; chromoNum<POOL_SIZE; chromoNum++){
			//DEBUG : Show Pre-adjust fitness
			//cout << " Fitness " <<  fitness[chromoNum] << endl;
			fitness[chromoNum] = pow((-fitness[chromoNum] + curWorstFitness),SELECT_PRESSURE);
			
			curTotalFitness += fitness[chromoNum];
		}
		
		// Calculate cumulative probability
		for (int chromoNum = 0; chromoNum<POOL_SIZE; chromoNum++){
			curCumProb += fitness[chromoNum] / curTotalFitness;
			cumulativeProb[chromoNum] = curCumProb;
		}
		//DEBUG : Show Cumulative Prob
		//for (int chromoNum = 0; chromoNum<POOL_SIZE; chromoNum++)			cout << cumulativeProb[chromoNum] << endl;
		
		//DEBUG : Show best and worst fitness every gen
		if (iter%UPDATE_INTERVAL==0){
			system("CLS");
		}
		else{
			HANDLE hOutput = ::GetStdHandle(STD_OUTPUT_HANDLE);

		   	COORD coord = {0,0};
		   	::SetConsoleCursorPosition(hOutput, coord);

		   	char buff[] = "\n\n\n\n\n\n";
		   	::WriteConsoleA(hOutput, buff, 6, NULL, NULL);
		   	
		   	::SetConsoleCursorPosition(hOutput, coord);
		}

    	chrono::duration<double> elapsed_seconds = chrono::system_clock::now()-start;
	
        cout << "   Elapsed time : " << setw(6) << setprecision(5) << fixed << elapsed_seconds.count() << "s\n";
        cout << "Avg time per gen: " << setw(6) << setprecision(5) << fixed << elapsed_seconds.count()/iter << "s" << endl;
		
		cout << "   Worst for gen " << setw(4) << iter << " : " << curWorstFitness << endl;
		cout << "    Best for gen " << setw(4) << iter << " : " << bestThisGen.evaluate(nodeList)<< endl;
		cout << "Best ever at gen " << setw(4) << bestEverAtGen<< " : " << bestEver.evaluate(nodeList)<< endl;
		cout << "Best updated at : ";
		
		for (int i = 0; i<bestTracker.size(); i++){
			if (i%22==0) cout << endl;
			cout << setw(6) << bestTracker[i] << ' ';
		}
		
		if (hConWnd && iter%UPDATE_INTERVAL==1)
		{
			for (int i = 0; i<CHROMO_LENGTH; i++){
				BCX_Line(hConWnd, nodeList[bestThisGen.path[i]].x*XScalar+Xoffset1, nodeList[bestThisGen.path[i]].y*YScalar+Yoffset, nodeList[bestThisGen.path[i==CHROMO_LENGTH-1?0:i+1]].x*XScalar+Xoffset1, nodeList[bestThisGen.path[i==CHROMO_LENGTH-1?0:i+1]].y*YScalar+Yoffset, Red);
				BCX_Circle(hConWnd, nodeList[i].x*XScalar+Xoffset1, nodeList[i].y*YScalar+Yoffset, CircleRadius, Lime);
				
				BCX_Line(hConWnd, nodeList[bestEver.path[i]].x*XScalar+Xoffset2, nodeList[bestEver.path[i]].y*YScalar+Yoffset, nodeList[bestEver.path[i==CHROMO_LENGTH-1?0:i+1]].x*XScalar+Xoffset2, nodeList[bestEver.path[i==CHROMO_LENGTH-1?0:i+1]].y*YScalar+Yoffset, Blue);
				BCX_Circle(hConWnd, nodeList[i].x*XScalar+Xoffset2, nodeList[i].y*YScalar+Yoffset, CircleRadius, Lime);
		  	}
		}

		
		//LOOP FOR SELECTION
		vector<Chromosome> newPool;
		
		for (int newGenPop = 0; newGenPop<POOL_SIZE; newGenPop+=2){
			float rand1 = RANDOM_FLOAT;
			float rand2 = RANDOM_FLOAT;
			
			//Parent Chromosome selection
			Chromosome cr1;
			Chromosome cr2;
			
			for (int searchChromo = 0; searchChromo < POOL_SIZE; searchChromo++){
				if (rand1 < cumulativeProb[searchChromo]){
					cr1 = Chromosome(chromoPool[searchChromo]);
					rand1 = 999;
				}
				if (rand2 < cumulativeProb[searchChromo]){
					cr2 = Chromosome(chromoPool[searchChromo]);
					rand2 = 999;
				}
				if (rand1 == 999 &&  rand2 == 999){
					break;
				}
			}

			//Edge Recombination Operator
			if (RANDOM_FLOAT < RECOM_CHANCE){
				
				//Reset Adjacency Matrix
				int adj_matrix[CHROMO_LENGTH][MAX_NEIGHBOURS];

				for (int i=0; i<CHROMO_LENGTH; i++){
					for (int j=0; j<MAX_NEIGHBOURS; j++){
						adj_matrix[i][j] = -1;
					}
				}
				
				updateAdjMatrix(cr1, adj_matrix);
				updateAdjMatrix(cr2, adj_matrix);

				//Copy adjacency matrix
				int adj_matrix_copy[CHROMO_LENGTH][MAX_NEIGHBOURS];
				for (int i = 0; i<CHROMO_LENGTH; i++){
					for (int j = 0; j<MAX_NEIGHBOURS; j++){
						adj_matrix_copy[i][j] = adj_matrix[i][j];
					}
				}

				Chromosome *temp1 = edge_recom(cr1, cr2, adj_matrix);
				Chromosome *temp2 = edge_recom(cr1, cr2, adj_matrix_copy);

				cr1 = *temp1;
				cr2 = *temp2;
			}
			
			newPool.push_back(cr1);
			newPool.push_back(cr2);
		}
		chromoPool.clear();
		chromoPool.assign(newPool.begin(), newPool.end());
  		newPool.clear();
  		
  		for (int i = 0; i<POOL_SIZE; i++){
  			if (RANDOM_FLOAT < MUT_CHANCE){
  				chromoPool[i].mutate();
  			}
  		}
	}
	
	//End of computation, display results
	system("CLS");

	chrono::duration<double> elapsed_seconds = chrono::system_clock::now()-start;

    cout << "   Elapsed time : " << setw(6) << setprecision(5) << fixed << elapsed_seconds.count() << "s\n";
    cout << "Avg time per gen: " << setw(6) << setprecision(5) << fixed << elapsed_seconds.count()/iter << "s" << endl;

	cout << "Best updated at : ";

	for (int i = 0; i<bestTracker.size(); i++){
		cout << setw(4) << bestTracker[i] << " ";
  	}
	
	cout <<endl;
	
	if (iter-bestTracker.back() <= HARD_THRESHOLD){
		cout << " Terminated at " << iter << endl;
	}
	
	std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	cout << " Started computation at " << std::ctime(&start_time) << endl;
	cout << "Finished computation at " <<std::ctime(&end_time) << endl;

	cout << "Best path found : " << endl << bestEver <<endl;
	
	if (hConWnd)
	{
		for (int i = 0; i<CHROMO_LENGTH; i++){
   			BCX_Line(hConWnd, nodeList[bestEver.path[i]].x*XScalar+(Xoffset1 +Xoffset2)/2, nodeList[bestEver.path[i]].y*YScalar+Yoffset, nodeList[bestEver.path[i==CHROMO_LENGTH-1?0:i+1]].x*XScalar+(Xoffset1 +Xoffset2)/2, nodeList[bestEver.path[i==CHROMO_LENGTH-1?0:i+1]].y*YScalar+Yoffset, Blue);
			BCX_Circle(hConWnd, nodeList[i].x*XScalar+(Xoffset1 +Xoffset2)/2, nodeList[i].y*YScalar+Yoffset, CircleRadius, Lime);
	  	}
	}
	
	cin.get();

 	return 0;
}
