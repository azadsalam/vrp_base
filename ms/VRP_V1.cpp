#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
using namespace std;

#define POPULATION_SIZE 20
#define NUMBER_OF_OFFSPRING 15

#define NUMBER_OF_GENERATION 20

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
	double* loadCapacity; // kon vehicle max koto load nite parbe
	double* serviceTime;  // kon client kototuk time lage service pete
	double* demand; 	  // kon client koto demand
	double** timeConstraintsOfVehicles; // periodCount * vehicleCount
void print()
{
    int i,j;
	cout<<"Period : "<<periodCount<<endl;
	cout<< "Number of depots : " <<depotCount <<endl;


	cout<<"Number of vehicles : "<< vehicleCount <<" Allocation : ";
	for( j=0;j<depotCount;j++) cout<<vehicleAllocation[j]<<" ";
	cout<<endl;

	cout<<"Vehicle to Depot Mapping : ";
	for( j=0;j<vehicleCount;j++) cout << depotAllocation[j]<<" ";
	cout << endl;

	cout<<"Load Capacity/ Max allowable load for each vehicle : ";
	for( i=0;i<vehicleCount;i++) cout << loadCapacity[i] << " ";
	cout<< endl;


    cout<<"Time constraints for vehicles : "<<endl;
    for( i=0;i<periodCount;i++)
    {
        for( j=0;j<vehicleCount;j++)
        {
            cout<<timeConstraintsOfVehicles[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    

	cout<<"Clients : "<<customerCount<<endl;

	cout << "Frequency allocation : ";
	for( i =0;i<customerCount ;i++) cout << frequencyAllocation[i] << " ";
	cout <<endl;

	cout << "Service Time : ";
	for( i =0;i<customerCount ;i++) cout << serviceTime[i] << " ";
	cout <<endl;

	cout << "Demand (load) : ";
	for( i =0;i<customerCount ;i++) cout << demand[i] << " ";
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
	double cost;
	bool isFeasible;
	bool feasibilitySet;

	double** loadViolation;
	double totalLoadViolation;

	ProblemInstance problemInstance;


	Individual()
	{
		periodAssignment = NULL;
		permutation = NULL;
		routePartition = NULL;
		cost = -1;
		feasibilitySet = false;
		isFeasible = false;

		loadViolation = NULL;

	}

	// make a copy cat individual
	//copies problem instance, periodAssignment, permutation, routePartition
	void makeCopy(Individual original)
	{
	    int i,j;
		this->problemInstance = original.problemInstance;


		periodAssignment = new bool*[problemInstance.periodCount];
		for( i=0;i<problemInstance.periodCount;i++)
		{
			periodAssignment[i] = new bool[problemInstance.customerCount];

			for( j=0;j<problemInstance.customerCount;j++)
			{
				periodAssignment[i][j] = original.periodAssignment[i][j];
			}
		}



		permutation = new int*[problemInstance.periodCount];
		for( i=0;i<problemInstance.periodCount;i++)
		{
			permutation[i] = new int[problemInstance.customerCount];

			for( j=0;j<problemInstance.customerCount;j++)
			{
				permutation[i][j] = original.permutation[i][j];
			}
		}


		routePartition = new int[problemInstance.vehicleCount];
		for( i=0;i<problemInstance.vehicleCount;i++)
		{
			routePartition[i]=original.routePartition[i];
		}

		feasibilitySet = original.feasibilitySet;
		cost = original.cost;
		isFeasible = original.isFeasible;

		//allocate demanViolationMatrix

        loadViolation = new double*[problemInstance.periodCount];
		for(int period=0; period < problemInstance.periodCount;period++)
        {
            loadViolation[period] = new double[problemInstance.vehicleCount];
        }

	}
	//calculate and return fitness of individual
	double calculateFitness()
	{
		double tempCost = 0;

		//totalLoadViolation = 0;
        double temlLoad;
		for(int i=0;i<problemInstance.periodCount;i++)
		{
			for(int j=0;j<problemInstance.vehicleCount;j++)
			{
				tempCost += calculateFitness(i,j);
                //calculate the total load violation
                //Add only when actually the load is violated i.e. violation is positive
                if(loadViolation[i][j]>0) totalLoadViolation += loadViolation[i][j];
			}
		}

		cost = tempCost;
		return cost;
	}

	//calcuate fitness for each period for each vehicle
	// route for vehicle i is  [ routePartition[i-1]+1 , routePartition[i] ]
	// given that routePartition[i-1]+1 <= routePartition[i]

	double calculateFitness(int period,int vehicle)
	{
		int assignedDepot;
		assignedDepot = problemInstance.depotAllocation[vehicle];
		double costForPV = 0;
		int start,end; // marks the first and last position of corresponding route for the array permutation

		if(vehicle == 0) start = 0;
		else start = routePartition[vehicle-1]+1;

		end = routePartition[vehicle];

		if(end<start) return 0;

		int activeStart=-1,activeEnd=-1,previous=-1,clientNode;

        double clientDemand=0;
		for(int i=start;i<=end;i++)
		{
			clientNode = permutation[period][i];
			if(periodAssignment[period][clientNode]==false) continue;

			if(activeStart == -1) activeStart = clientNode;
			activeEnd = clientNode;

            //Caluculate total client demand for corresponding period,vehicle
            clientDemand += problemInstance.demand[clientNode];

			if(previous == -1)
			{
				previous = clientNode;
				continue;
			}

			costForPV +=	problemInstance.costMatrix[previous+problemInstance.depotCount][clientNode+problemInstance.depotCount];

			previous = clientNode;

		}

        if(activeStart!=-1 && activeEnd != -1)
        {
            costForPV += problemInstance.costMatrix[assignedDepot][activeStart+problemInstance.depotCount];
            costForPV += problemInstance.costMatrix[activeEnd+problemInstance.depotCount][assignedDepot];
        }
        loadViolation[period][vehicle] = clientDemand - problemInstance.loadCapacity[vehicle];

		return costForPV;

	}

	void initialise(ProblemInstance problemInstance)
	{
	    int i,j;
		this->problemInstance = problemInstance;

		// ALLOCATING periodCount * customerCount Matrix for Period Assignment
		periodAssignment = new bool*[problemInstance.periodCount];
		for( i=0; i< problemInstance.periodCount; i++)
		{
			periodAssignment[i] = new bool[problemInstance.customerCount];

			//initialise every period assignment with false
			for( j=0;j<problemInstance.customerCount;j++)
				periodAssignment[i][j] = false;
		}


		//ALLOCATING permutation map matrix -> period * customer
		permutation = new int*[problemInstance.periodCount];
		for( i=0;i<problemInstance.periodCount;i++)
		{
			permutation[i] = new int[problemInstance.customerCount];

			// initially every permutation is identity permutation
			for( j=0;j<problemInstance.customerCount;j++)
			{
				permutation[i][j] = j;
			}
		}

		//allocating routeAllocation
		routePartition = new int[problemInstance.vehicleCount];


		// NOW INITIALISE WITH VALUES

		//initialize period assignment

		int freq,allocated,random,tmp;

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
			for( i = problemInstance.customerCount -1 ;i>0 ;i-- )
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



		loadViolation = new double*[problemInstance.periodCount];
		for( period=0; period < problemInstance.periodCount;period++)
        {
            loadViolation[period] = new double[problemInstance.vehicleCount];
        }

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
		int i,j;
		cout << "PERIOD ASSIGMENT : \n";
		for( i=0;i<problemInstance.periodCount;i++)
		{
			for( j=0;j<problemInstance.customerCount;j++)
			{
				cout<<periodAssignment[i][j]<<" ";
			}
			cout << endl;
		}

		cout << "Permutation : \n";
		for( i=0; i<problemInstance.periodCount;i++)
		{
			for( j=0;j<problemInstance.customerCount;j++)
			{
				cout << permutation[i][j]<<" ";
			}
			cout << endl;
		}

		cout<< "Route partition : ";
		for( i=0;i<problemInstance.vehicleCount;i++)cout<< routePartition[i] <<" ";
		cout << endl;

        // print load violation

		cout<<endl<<endl<<"LOAD VIOLATION MATRIX : \n";
        for( i=0;i<problemInstance.periodCount;i++)
        {
            for( j=0;j<problemInstance.vehicleCount;j++)
            {
                cout << loadViolation[i][j]<<" ";
            }
            cout<<endl;
        }
        cout<<"Total Load Violation : "<<totalLoadViolation<<endl;


		cout << "\nFitness/Cost : " << cost << endl <<endl;
	}

	void mutatePermutation(int period)
	{
		int first = rand() % problemInstance.customerCount;

		int second;
		int count=0;
		do
		{
			second = rand() % problemInstance.customerCount;
			count++;
			if(count==problemInstance.customerCount)break;
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
		//nothing to do if only one vehicle
		if(problemInstance.vehicleCount == 1) return ;

		//pick a red line/seperator
		//generate random number in [0,vehicleCount-1)


		int distance,increment;

		while(true)
		{
			int seperatorIndex = rand() % (problemInstance.vehicleCount-1);
			int dir = rand() %2; // 0-> left , 1-> right
			if(dir==0)//move the seperator left
			{
				if(seperatorIndex==0) distance = routePartition[0] ;
				else distance = routePartition[seperatorIndex] - routePartition[seperatorIndex-1];
				// if the line can not merge with the previous one ,
				// difference = routePartition[seperatorIndex] - 1 - routePartition[seperatorIndex-1]

				// increment should be in [1,distance]
				if(distance==0)continue;
				increment = 1 + ( rand() % distance );
				routePartition[seperatorIndex] -= increment;
				return;
			}
			else	//move the seperator right
			{
				distance = routePartition[seperatorIndex+1] - routePartition[seperatorIndex] ;
				if(distance==0)continue;
				increment = 1 + (rand() % distance );
				routePartition[seperatorIndex] += increment;
				return;
			}
		}

	}


	//returns 0 if it couldnt mutate as period == freq
	int mutatePeriodAssignment(int clientNo)
	{
		//no way to mutate per. ass. as freq. == period
		if(problemInstance.frequencyAllocation[clientNo] == problemInstance.periodCount) return 0;

		int previouslyAssigned; // one period that was assigned to client
		do
		{
			previouslyAssigned = rand() % problemInstance.periodCount;
		} while (periodAssignment[previouslyAssigned][clientNo]==false);

		int previouslyUnassigned;//one period that was NOT assigned to client
		do
		{
			previouslyUnassigned = rand() % problemInstance.periodCount;
		} while (periodAssignment[previouslyUnassigned][clientNo]==true);

		periodAssignment[previouslyAssigned][clientNo] = false;
		periodAssignment[previouslyUnassigned][clientNo]= true;

		return 1;
	}
};

class GeneticAlgo
{
	ProblemInstance problemInstance;
	Individual population[POPULATION_SIZE];

	//for storing new generated offsprings
	Individual offspringPopulation[NUMBER_OF_OFFSPRING];

	//for temporary storing
	Individual temporaryPopulation[POPULATION_SIZE];

	// for selection - roulette wheel
	double fitness[POPULATION_SIZE];
	double cdf[POPULATION_SIZE];


public:
	GeneticAlgo(ProblemInstance problemInstance)
	{
		this -> problemInstance = problemInstance;
	}
	void run()
	{
		int selectedParent;

		Individual parent,offspring;

		//problemInstance.print();

		// INITIALISE POPULATION
		initialisePopulation();

		sort(population,POPULATION_SIZE);

    int i;
		for(int generation=0;generation<NUMBER_OF_GENERATION;generation++)
		{
			//sort function uses selection sort, replace with some O(n lg n) sort algthm

			cout << "--------------------------\nGENERATION : "<<generation<<"\n\n";

			//Select a parent and apply genetic operator
			for( i=0;i<NUMBER_OF_OFFSPRING;i++)
			{

					selectedParent=rouletteWheelSelection();

					parent = population[selectedParent];
					offspring.makeCopy(parent);


					applyMutation(offspring);
					cout << "Selected Parent : " << selectedParent <<endl;
					parent.print();
					offspring.calculateFitness();
					cout << "Offspring : \n";
					offspring.print();

					offspringPopulation[i] = offspring;

			}


			cout <<"\n\n\n\n\n---TESTING (lamba + miu)"<<endl<<endl;

			//////////////////////////////////////////////////////////////////////////
			cout <<"PRINTING PARENT POPULATION\n";
			for( i=0;i<POPULATION_SIZE;i++)
			{
				cout << "parent " << i << " :\n";
				population[i].print();
			}
			cout << endl<<endl;
			//////////////////////////////////////////////////////////////////////////////

			//TAKE THE BEST "POPULATION_SIZE" individuals from the set of all parents and children
			sort(offspringPopulation,NUMBER_OF_OFFSPRING);


			//////////////////////////////////////////////////////////////////////////////
			cout <<"PRINTING OFFSPRING POPULATION\n";
			for( i=0;i<NUMBER_OF_OFFSPRING;i++)
			{
				cout << "offspring " << i << " :\n";
				offspringPopulation[i].print();
			}
			cout << endl<<endl;
			//////////////////////////////////////////////////////////////////////////////

			//first select best indivdls in the temporary array
			//afterwards replace population with it
			i = 0;
			int j = 0;
			int cursor = 0;

			while(cursor<POPULATION_SIZE)
			{
				if(i == POPULATION_SIZE)
				{
					temporaryPopulation[cursor] = offspringPopulation[j];
					j++;
				}
				else if(j== NUMBER_OF_OFFSPRING)
				{
					temporaryPopulation[cursor] = population[i];
					i++;
				}
				else if(population[i].cost <= offspringPopulation[j].cost)
				{
					temporaryPopulation[cursor] = population[i];
					i++;
				}
				else
				{
					temporaryPopulation[cursor] = offspringPopulation[j];
					j++;
				}
				cursor++;
			}

			//replace population with temporary array
			for(i=0;i<POPULATION_SIZE;i++)
			{
				population[i] = temporaryPopulation[i];
			}


			//////////////////////////////////////////////////////////////////////////
			cout <<"PRINTING NEW GENERATION\n";
			for( i=0;i<POPULATION_SIZE;i++)
			{
				cout << "parent " << i << " :\n";
				population[i].print();
			}
			cout << endl<<endl;
			//////////////////////////////////////////////////////////////////////////////



		}


		cout<<"\n\n\n\n\n--------------------------------------------------\n";
		cout<<"FINAL POPULATION"<<endl;
		for( i=0;i<POPULATION_SIZE;i++)
		{
			cout<<"Individual : "<<i<<endl;
			population[i].print();
		}

	}


	// for now not applying periodAssignment Mutation operator
	// for now working with only MDVRP ->  period = 1
	void applyMutation(Individual offspring)
	{


		int selectedMutationOperator = selectMutationOperator();
		if(selectedMutationOperator==0)offspring.mutateRoutePartition();
		else if (selectedMutationOperator == 1)offspring.mutatePermutation(0);//for now single period
		else
		{
			offspring.mutateRoutePartition();
			offspring.mutatePermutation(0);
		}

	}

	//0 -> route partition
	//1 ->	permutation
	//2 -> both
	int selectMutationOperator()
	{
		return rand()%3;
	}

	//SORT THE INDIVIDUALS ON ASCENDING ORDER OF COST
	//BETTER INDIVIDUALS HAVE LOWER INDEX
	//COST LESS, INDEX LESS ;-)
	void sort(Individual* array,int length)
	{
		Individual temp;
		//FOR NOW DONE SELECTION SORT
		//AFTERWARDS REPLACE IT WITH QUICK SORT OR SOME OTHER O(n logn) sort
		for(int i=0;i<length;i++)
		{
			for(int j=i+1;j<length;j++)
			{
				if(array[i].cost >array[j].cost)
				{
					temp = array[i];
					array[i] =array[j];
					array[j] = temp;
				}
			}
		}

	}


	// it also calculates cost of every individual
	int rouletteWheelSelection()
	{
    int i,j;
		//SELECTION -> Roulette wheel
		double sumOfFitness = 0,sumOfCost=0;
		double sumOfProability = 0;

		cout<< "SELECTION\nCost : ";
		for( i=0;i<POPULATION_SIZE;i++)
		{
			population[i].calculateFitness();
			fitness[i] = population[i].cost;
			sumOfCost += fitness[i];
			cout << " "<<fitness[i];
		}
		cout <<"   Total cost : "<<sumOfCost<<endl;

//		cout<< "Fitness : ";
		for( i=0;i<POPULATION_SIZE;i++)
		{
			fitness[i] = sumOfCost / fitness[i]; // the original fitness
			// incorporate penalty
			sumOfFitness += fitness[i];
		//	cout << " "<< fitness[i];
		}
		//cout <<"    Total fitness : "<<sumOfFitness<<endl;

		for( i=0;i<POPULATION_SIZE;i++)
		{
			sumOfProability = cdf[i] = sumOfProability + ((double)fitness[i]/sumOfFitness);
		}

		double num = rand()%101; // generate random number from [0,100]
		double indicator = num/100;

		//find the smallest index i, with cdf[i] greater than indicator

		int par =  findParent(indicator);
		cout <<"Selected Parent : "<< par<<endl;
		return par;

	}

	//binary search for smallest index i, having cdf[i] greater than indicator
	int findParent(double indicator)
	{
		//for now linear search, do binary search later
		for(int i=0;i<POPULATION_SIZE;i++)
			if(cdf[i]>=indicator)
				return i;
		return POPULATION_SIZE-1;
	}

	void initialisePopulation()
	{
		cout << "Initial population : \n";
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
    int i,j;
	for(int t=0 ; t<testCases;t++)
	{
		cout<< "TEST CASE "<<t+1<<" : "<<endl<<endl;



		in>>problemInstance.periodCount;
		escapeComment();

		in>>problemInstance.depotCount;
		escapeComment();

		in>>problemInstance.vehicleCount;
		escapeComment();


		//vehicle per depot
		problemInstance.vehicleAllocation = new int[problemInstance.depotCount];
		problemInstance.depotAllocation = new int[problemInstance.vehicleCount];
		int vehicleCursor = 0;

		for( j=0;j<problemInstance.depotCount;j++)
		{
			in>>problemInstance.vehicleAllocation[j];

			for( i=0;i<problemInstance.vehicleAllocation[j];i++)
			{
				problemInstance.depotAllocation[vehicleCursor]=j;
				vehicleCursor++;
			}
		}
		escapeComment();

		//capacity of vehicle
		problemInstance.loadCapacity = new double[problemInstance.vehicleCount];
		for( i=0;i<problemInstance.vehicleCount;i++)
		{
			in >> problemInstance.loadCapacity[i];
		}
		escapeComment();

        
        //time constraints
        escapeComment(); // escape the line "; t(total period) lines containg  v (total vehicle)
                         //values each referring maximum time limit for that day for that vehicle (NEW)"

        //read periodCount lines
        problemInstance.timeConstraintsOfVehicles = new double* [problemInstance.periodCount];
        for(int i=0;i<problemInstance.periodCount;i++)
        {
            problemInstance.timeConstraintsOfVehicles[i] = new double[problemInstance.vehicleCount];
            //read vehicleCOunt values
            for(int j=0;j<problemInstance.vehicleCount;j++)
            {
                in>>problemInstance.timeConstraintsOfVehicles[i][j];
            }
        }
            
		//CLIENT COUNT
		in>>problemInstance.customerCount;
		escapeComment();

		//frequency
		problemInstance.frequencyAllocation = new int[problemInstance.customerCount];
		for( i=0 ; i<problemInstance.customerCount; i++) in>> problemInstance.frequencyAllocation[i];
		escapeComment();


		//service time
		problemInstance.serviceTime = new double[problemInstance.customerCount];
		for( i=0; i<problemInstance.customerCount; i++)
		{
			in >> problemInstance.serviceTime[i];
		}
		escapeComment();

		//demand
		problemInstance.demand = new double[problemInstance.customerCount];
		for( i=0 ; i<problemInstance.customerCount; i++)
		{
			in >> problemInstance.demand[i];
		}
		escapeComment();

		// cost matrix
		escapeComment(); // escapes the line ";escape comment"
		problemInstance.nodeCount = problemInstance.customerCount+problemInstance.depotCount;
		problemInstance.costMatrix = new double*[problemInstance.nodeCount];
		for( i=0;i<problemInstance.nodeCount;i++)problemInstance.costMatrix[i] = new double[problemInstance.nodeCount];


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
	//srand(1);
	ofstream file;
	file.open ("out.txt");
	streambuf* sbuf = cout.rdbuf();
	cout.rdbuf(file.rdbuf());

	cout << "VRP CPP V1" << endl;



	 	parseInputFile();

    cout<<"END\n";
    cout.rdbuf(sbuf);

	return 0;
}

