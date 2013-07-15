#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
using namespace std;
/*
class ProblemInstance
{
	int ncustomers;
	int ndepots;
	int nnodes;
	int nperiods;
	int nvehicles;
	int* vehicleAllocation; // HOW MANY VEHICLES PER DEPOT
	int** cost;

public:
	ProblemInstance(int ncustomers, int ndepots, int nperiods, int nvehicles,
			int* vehicleAllocation, int** cost)
	{
		this->ncustomers = ncustomers;
		this->ndepots = ndepots;
		this->nperiods = nperiods;
		this->nvehicles = nvehicles;
		this->nnodes = ndepots + ncustomers;

		this->vehicleAllocation = new int[ndepots];
		for (int i = 0; i < ndepots; i++)
		{
			this->vehicleAllocation[i] = vehicleAllocation[i];
		}

		this->cost = new int*[nnodes];

		for (int i = 0; i < nnodes; i++)
		{
			this->cost[i] = new int[nnodes];
		}

		for (int i = 0; i < nnodes; i++)
		{
			for (int j = 0; j < nnodes; j++)
			{
				this->cost[i][j] = cost[i][j];
			}
		}
	}

};

*/

#define POPULATION_SIZE 1
#define NUMBER_OF_GENERATION 5

int testCases ;
ifstream in ("input.txt");


void escapeComment()
{
	char comment[255];
	in.getline(comment,255);

}

class ProblemInstance
{
public:
	int depotCount,customerCount,periodCount;
	int* vehicleAllocation; // kon depot er koyta kore vehicle
	double **costMatrix;
	int nodeCount;
	int *frequencyAllocation;
	int vehicleCount;
	int* depotAllocation; // kon vehicle kon depot er under a


void print()
{
	cout<<"Period : "<<periodCount<<endl;
	cout<< "Number of depots : " <<depotCount <<endl;


	cout<<"Number of vehicles : "<< vehicleCount <<" Allocation : ";
	for(int j=0;j<depotCount;j++) cout<<vehicleAllocation[j]<<" ";
	cout<<endl;

	cout<<"Vehicle to Depot Mapping : ";
	for(int j=0;j<vehicleCount;j++) cout << depotAllocation[j]<<" ";
	cout << endl;

	cout<<"Clients : "<<customerCount<<endl;

	cout << "Frequency allocation : ";
	for(int i =0;i<customerCount ;i++) cout << frequencyAllocation[i] << " ";
	cout <<endl;

	cout<< "Printing cost matrix : \n";

	for(int row=0;row<nodeCount;row++)
	{
		for(int col=0;col<nodeCount;col++)
			cout<<costMatrix[row][col]<<" ";
		cout<<endl;

	}
	cout<<endl;
}

};

class Individual
{
public:
	bool** periodAssignment;
	int** permutation;
	int* routePartition;
	double fitness;
	bool isFeasible;
	bool feasibilitySet;

	ProblemInstance problemInstance;


	Individual()
	{
		periodAssignment = NULL;
		permutation = NULL;
		routePartition = NULL;
		fitness = -1;
		feasibilitySet = false;
		isFeasible = false;
	}

	//calculate and return fitness of individual
	double calculateFitness()
	{
		double cost = 0;
		for(int i=0;i<problemInstance.periodCount;i++)
		{
			for(int j=0;j<problemInstance.vehicleCount;j++)
			{
				cost += calculateFitness(i,j);
			}
		}

		fitness = cost;
		return fitness;
	}

	//calcuate fitness for each period for each vehicle
	// route for vehicle i is  [ routePartition[i-1]+1 , routePartition[i] ]
	// given that routePartition[i-1]+1 <= routePartition[i]

	double calculateFitness(int period,int vehicle)
	{
		int assignedDepot;
		assignedDepot = problemInstance.depotAllocation[vehicle];
		double cost = 0;
		int start,end; // marks the first and last position of corresponding route for the array permutation

		if(vehicle == 0) start = 0;
		else start = routePartition[vehicle-1]+1;

		end = routePartition[vehicle];

		if(end<start) return 0;

		int activeStart=-1,activeEnd,previous=-1,clientNode;


		for(int i=start;i<=end;i++)
		{
			clientNode = permutation[period][i];
			if(periodAssignment[period][clientNode]==false) continue;

			if(activeStart == -1) activeStart = clientNode;
			activeEnd = clientNode;

			if(previous == -1)
			{
				previous = clientNode;
				continue;
			}

			cost +=	problemInstance.costMatrix[previous+problemInstance.depotCount][clientNode+problemInstance.depotCount];

			previous = clientNode;

		}

		cost += problemInstance.costMatrix[assignedDepot][activeStart+problemInstance.depotCount];
		cost += problemInstance.costMatrix[activeEnd+problemInstance.depotCount][assignedDepot];


		return cost;

	}

	void initialise(ProblemInstance problemInstance)
	{
		this->problemInstance = problemInstance;

		// ALLOCATING periodCount * customerCount Matrix for Period Assignment
		periodAssignment = new bool*[problemInstance.periodCount];
		for(int i=0; i< problemInstance.periodCount; i++)
		{
			periodAssignment[i] = new bool[problemInstance.customerCount];

			//initialise every period assignment with false
			for(int j=0;j<problemInstance.customerCount;j++)
				periodAssignment[i][j] = false;
		}


		//ALLOCATING permutation map matrix -> period * customer
		permutation = new int*[problemInstance.periodCount];
		for(int i=0;i<problemInstance.periodCount;i++)
		{
			permutation[i] = new int[problemInstance.customerCount];

			// initially every permutation is identity permutation
			for(int j=0;j<problemInstance.customerCount;j++)
			{
				permutation[i][j] = j;
			}
		}

		//allocating routeAllocation
		routePartition = new int[problemInstance.vehicleCount];


		// NOW INITIALISE WITH VALUES

		//initialize period assignment

		int freq,allocated,random,tmp,j;

		//Randomly allocate period to clients equal to their frequencies
		for(int client=0; client < problemInstance.customerCount; client++)
		{
			freq = problemInstance.frequencyAllocation[client];
			allocated=0;

			while(allocated!=freq)
			{
				random = rand() % problemInstance.periodCount;
				if(periodAssignment[random][client]==false)
				{
					periodAssignment[random][client]=true;
					allocated++;
				}
			}
		}

		//initialize permutation map - KNUTH SHUFFLE
		for(int period=0; period < problemInstance.periodCount;period++)
		{
			//apply knuths shuffle
			for(int i = problemInstance.customerCount -1 ;i>0 ;i-- )
			{
				j = rand() % i ;

				//if(i == j) continue;

				tmp = permutation[period][i];
				permutation[period][i] = permutation[period][j];
				permutation[period][j] = tmp;
			}
		}


		//NEED TO GENERATE #vehicle-1 (not distinct - distinct) random numbers in increasing order from [0,#customer - 1]
		// DEVICE some faster and smarter algorithm

		// route for vehicle i is  [ routePartition[i-1]+1 , routePartition[i] ]
		// given that routePartition[i-1]+1 <= routePartition[i]

		//bool found;
		allocated = 0;
		while(allocated != problemInstance.vehicleCount-1)
		{
			random = rand() % problemInstance.customerCount;

			//found = binarySearch(random,allocated-1);

			//if(found) continue;
			//else
			//{
			routePartition[allocated]=random;
			sort(random,allocated);
			allocated++;

			//}
		}
		routePartition[problemInstance.vehicleCount-1] = problemInstance.customerCount-1;

		calculateFitness();
	}

	// sorts the array routePartition in increasing order
	// input -> routePartition array [0, upperbound ], with,n inserted at the last in the array
	// output -> sorted array [0, upperbound]
	void sort(int n,int upperbound)
	{
		int tmp;
		for(int i = upperbound-1;i>=0;i--)
		{
			if(routePartition[i]>routePartition[i+1])
			{
				tmp = routePartition[i];
				routePartition[i] = routePartition[i+1];
				routePartition[i+1] = tmp;
			}
			else
				break;
		}
	}
	//Searches if the number already exists or not, within the given range [0,upperBound] of routePartition
	bool binarySearch(int n, int upperBound)
	{
		if(upperBound<0) return false;

		int low=0;
		int high=upperBound;
		int mid;

		while(low<=high)
		{
			mid = (high+low)/2;

			if(routePartition[mid]==n) return true;
			else if(n < routePartition[mid])
				high = mid-1;
			else
				low = mid+1;
		}

		return false;
	}

	void print()
	{
		cout << "PERIOD ASSIGMENT : \n";
		for(int i=0;i<problemInstance.periodCount;i++)
		{
			for(int j=0;j<problemInstance.customerCount;j++)
			{
				cout<<periodAssignment[i][j]<<" ";
			}
			cout << endl;
		}

		cout << "Permutation : \n";
		for(int i=0; i<problemInstance.periodCount;i++)
		{
			for(int j=0;j<problemInstance.customerCount;j++)
			{
				cout << permutation[i][j]<<" ";
			}
			cout << endl;
		}

		cout<< "Route partition : ";
		for(int i=0;i<problemInstance.vehicleCount;i++)cout<< routePartition[i] <<" ";
		cout << endl;
		cout << "Fitness/Cost : " << fitness << endl <<endl;
	}

	void mutatePermutation(int period)
	{
		int first = rand() % problemInstance.customerCount;

		int second;
		do
		{
			second = rand() % problemInstance.customerCount;
		}
		while(periodAssignment[period][second]==false || second == first);

		int temp = permutation[period][first];
		permutation[period][first] = permutation[period][second];
		permutation[period][second] = temp;

		// FITNESS CAN BE UPDATED HERE
	}

	//moves some red line
	//no effect if only one vehicle
	void mutateRoutePartition()
	{
		//pick a red line/seperator
		//generate random number in [0,vehicleCount-1)

		if(problemInstance.vehicleCount == 1) return ;
		int seperatorIndex = rand() % (problemInstance.vehicleCount-1);
		int dir = rand() %2; // 0-> left , 1-> right

		//if(seperatorIndex == 0) dir = 1; // 0th seperator will always go right

		int difference,increment;
		if(dir==0)//move the seperator left
		{
			if(seperatorIndex==0) difference = routePartition[0] ;
			else difference = routePartition[seperatorIndex] - routePartition[seperatorIndex-1] ;
			increment = 1 + ( rand() % difference );
			routePartition[seperatorIndex] -= increment;
		}
		else	//move the seperator right
		{
			difference = routePartition[seperatorIndex+1] - routePartition[seperatorIndex] ;
			increment = 1 + (rand() % difference );
			routePartition[seperatorIndex] += increment;
		}


	}

};

class GeneticAlgo
{
	ProblemInstance problemInstance;
	Individual population[POPULATION_SIZE];

public:
	GeneticAlgo(ProblemInstance problemInstance)
	{
		this -> problemInstance = problemInstance;
	}
	void run()
	{
		//problemInstance.print();

		// INITIALISE POPULATION
		initialisePopulation();

		for(int generation=0;generation<NUMBER_OF_GENERATION;generation++)
		{

			cout << "GENERATION : "<<generation<<" \n";
			population[0].mutateRoutePartition();
			population[0].calculateFitness();
			population[0].print();

		}
	}

	void initialisePopulation()
	{
		for(int i=0; i<POPULATION_SIZE; i++)
		{
			population[i].initialise(problemInstance);
			cout<<"Printing individual "<<i <<" : "<<endl;
			population[i].print();
		}
	}
};
ProblemInstance problemInstance;

void parseInputFile()
{

	cout << "OPENING INPUT FILE" << endl;
	if(!in)
	{
		cout << "ERROR OPENING FILE!! EXIT !!"<< endl;
		exit(1);
	}


	in>>testCases;
	cout<<"TEST CASES : "<<testCases<<endl;
	escapeComment();

	for(int t=0 ; t<testCases;t++)
	{
		cout<< "TEST CASE "<<t+1<<" : "<<endl<<endl;

		in>>problemInstance.periodCount;
		escapeComment();

		in>>problemInstance.depotCount;
		escapeComment();

		in>>problemInstance.vehicleCount;
		escapeComment();


		problemInstance.vehicleAllocation = new int[problemInstance.depotCount];
		problemInstance.depotAllocation = new int[problemInstance.vehicleCount];
		int vehicleCursor = 0;

		for(int j=0;j<problemInstance.depotCount;j++)
		{
			in>>problemInstance.vehicleAllocation[j];

			for(int i=0;i<problemInstance.vehicleAllocation[j];i++)
			{
				problemInstance.depotAllocation[vehicleCursor]=j;
				vehicleCursor++;
			}
		}
		escapeComment();

		in>>problemInstance.customerCount;
		escapeComment();

		problemInstance.frequencyAllocation = new int[problemInstance.customerCount];

		for(int i=0 ; i<problemInstance.customerCount; i++) in>> problemInstance.frequencyAllocation[i];
		escapeComment();

		problemInstance.nodeCount = problemInstance.customerCount+problemInstance.depotCount;

		problemInstance.costMatrix = new double*[problemInstance.nodeCount];
		for(int i=0;i<problemInstance.nodeCount;i++)problemInstance.costMatrix[i] = new double[problemInstance.nodeCount];


		for(int row=0;row<problemInstance.nodeCount;row++)
			for(int col=0;col<problemInstance.nodeCount;col++)
				in>>problemInstance.costMatrix[row][col];

		problemInstance.print();

		GeneticAlgo ga(problemInstance);
		ga.run();
	}
	in.close();
}





int main()
{
	srand (time(NULL));
	cout << "VRP CPP V1" << endl;
	parseInputFile();

	return 0;
}
