#include <iostream>
#include <numeric>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

class Graph
{
private:
	int size, rando, maxWeight;
	int ctr, currentRow, sumEdges, totalCost = 0;
	float gen, prob, density;
	vector <int> sourceRow, MST, Weights, closedSet, values, node;
	vector < vector <int> > openSet, graph, weightedGraph;
	int dest;
	bool isDestination = false; //Initially we assume the next node placed in the closed set is not the destination

public:
	Graph(int nodes, float density) //Default constructor for calling Graph class
	{
		this->size = nodes; //Set size of adjacency matrix to the size passed in from file - 20 in the case of the example
		graph.resize(size); //Resize the graph & weighted graph to the same size
		weightedGraph.resize(size);
		values.resize(3);

		for (int i = 0; i < size; ++i)
		{
			graph[i] = vector <int>(size); //Set up graph vectors. After fully populating all the vectors
			weightedGraph[i] = vector <int>(size);
		}

		for (int i = 0; i < size; ++i) //This nested loop will set up an adjacency matrix representing a random graph
		{
			for (int j = 0; j < size; ++j)
			{
				graph[i][j] = graph[j][i] = false; //Create a zero by default
				weightedGraph[i][j] = weightedGraph[j][i] = false;

				if (i == j) //Create zeros along diagonal
				{
					graph[i][j] = graph[j][i] = false;
					weightedGraph[i][j] = weightedGraph[j][i] = false;
				}

				else
				{
					prob = ((float)rand() / (RAND_MAX)); // Define probability function using RAND_MAX
					graph[i][j] = graph[j][i] = (prob < density); //Create a 1 if prob is less than density
				}
			}
		}
	};

	void addWeights(int maxWeight) //This function takes the random graph and updates the existing edges to have weights 
	{                                     //The maximum weight will be passed in by the calling function
		this->maxWeight = maxWeight;

		for (int i = 0; i < size; ++i)
		{
			for (int j = i; j < size; ++j)
			{
				gen = (float)rand() / (RAND_MAX)* maxWeight + 1; //Create random weight between 1 and maximum weight
				rando = (int)gen; // Cast the random weight to an integer
				weightedGraph[i][j] = weightedGraph[j][i] = graph[i][j] * rando; // Multiply the edge of the graph by the weight of the edge
			}
		}
	};

	void addToOpenSet(vector <int> node) //Created this function just so when it's called later its purpose is more clear
	{
		openSet.push_back(node);
	};

	void analyzeEdges(vector <int> sourceRow) //This function analyzes the possible paths from the current node
	{
		node.resize(2);
		for (int i = 0; i < size; ++i) //This loop will find the non-zero edges and add the corresponding node/edge pairs to the open set
		{
			if (sourceRow[i] != 0)
			{
				node = { i, sourceRow[i] }; //Create node/edge pair
				if (find(closedSet.begin(), closedSet.end(), i) == closedSet.end())//If i is not in closed set:
				{
					if (find(openSet.begin(), openSet.end(), node) == openSet.end())//If node is not in open set:
						node[1] = node[1] + totalCost; //Add distance traveled so far to edge cost
						addToOpenSet(node); //Add node/edge pair to open set
				}
			}
		}
	};

	vector<vector<int>> evalOpenSet(vector < vector < int > > openSet) //This function will search the open set and replace a node if a lower cost is found for the same node
	{
		vector <int> Nodes(openSet.size()); //Create node list from open set
		vector <int> foundNode(2);// Create a node placeholder in case a current node needs to be replaced

		for (int i = 0; i < openSet.size(); ++i)
		{
			for (int i = 0; i < openSet.size(); ++i)// List nodes currently in open set without weights
			{
				Nodes[i] = openSet[i][0]; //Create node list before each iteration
			};

			int nodeName = Nodes[i];
			auto it = find(Nodes.begin() + (i + 1), Nodes.end(), nodeName);//Look for nodeName in remainder of open set
			if (it != Nodes.end()) //If the node is present in the open set:
			{
				ptrdiff_t pos = distance(Nodes.begin(), it);

				if (pos >= openSet.size()) //If the position of the found node is the last node there are no duplicates
				{
					break;
				}

				node = openSet[i]; //Pull out node at current first location                          
				foundNode = openSet[pos]; //Pull out node at second location
				if (node[1] < foundNode[1]) //If original node cost is smaller than the newer node:
				{
					openSet.erase(openSet.begin() + pos); //Remove newer node from open set
				}
				else
				{
					openSet.erase(openSet.begin() + i); //Otherwise remove original node from open set and replace with newer node of equal or lower cost
				}

			}

			this->openSet = openSet; //Update the open set for the next pass
		}

		return openSet;
	};

	//Next steps: Populate closed set with lowest cost node and check if that node is the destination
	bool selectNode(vector <vector <int>> openSet)
	{
		if (ctr == 0) //Put the source node in closed set to begin
		{
			closedSet.push_back(0);
		};

		openSet = evalOpenSet(openSet); //Evaluate the open set for redundant nodes and update

		vector <int> edgeCosts(openSet.size()); //Convert the current open set to a vector 

		for (int i = 0; i < openSet.size(); ++i)
		{
			edgeCosts[i] = openSet[i][1]; //Create a vector to keep track of the edge costs to each node in the open set
		};

		if (openSet.size() != 0) //If open set is empty, there are no possible edges and the code below will crash
		{
			int indexOfLeast = distance(edgeCosts.begin(), min_element(edgeCosts.begin(), edgeCosts.end())); // Returns the index of the smallest edge cost
			vector <int> nextNode = openSet[indexOfLeast]; //Choose the next node out of the open set, located at the index of the smallest next edge cost

			closedSet.push_back(nextNode[0]); //Add node to MST and therefore the closed set
			Weights.push_back(nextNode[1]);
			totalCost = nextNode[1];
			openSet.erase(openSet.begin() + indexOfLeast); //Remove current node from open set
			this->openSet = openSet; //Update global open set with the current removed
			currentRow = nextNode[0]; //Update value of current row to represent traveling to the next node

		    //Need to check if node moved to closed set is the destination. If so, set isDestination to TRUE
			if (nextNode[0] == dest)
			{
				isDestination = true;
			}
		}
		else
		{

		}

		ctr++;
		return isDestination;
	}

	bool isConnected() //This function will check if the graph is connected. It will detect any grouping of nodes that are an island.
	{
		bool connected = 1; //Assume graph is connected. If algorithm detects graph is disconnected the value of connected will be set to zero.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				sumEdges += graph[i][j]; //Add up all of the non-weighted edges in the graph. These will show up as 1's in the adjacency matrix.
			}
		}

		sumEdges = sumEdges / 2; //Each edge is counted twice in the adjacency matrix, so we divide by 2 to find the actual number of edges

		if (sumEdges >= size - 1) //For n bodes, if there are less than n-1 edges the graph is disconnected and we can set connected equal to 0
		{
			sourceRow.resize(size); //Allocate memory for sourceRow
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < size; ++j)
				{
					sourceRow[j] = graph[i][j];
				}
				if (accumulate(sourceRow.begin(), sourceRow.end(), 0) == 0) //If there are no edges, the sum will be zero and the node is an island.
				{
					connected = 0; //Set the value of connected to zero
					break; //Break out of isConnected() function
				}
			};

			for (int i = 1; i < size - 1; ++i)
			{
				//The following code checks for islands. An island will appear is there is one and only one link between a pair of nodes.
				vector <int> checkRow = graph[0];
				vector <int> reverseRow = graph[graph.size() - i]; //Start at last row and work backwards
				reverse(reverseRow.begin(), reverseRow.end());
				if (reverseRow == sourceRow)
				{
					connected = 0;
					break;
				}
			}
		}


		else
		{
			connected = 0;
		}
		return connected;

	};

	void printGraph()
	{
		cout << "Adjacency Matrix Representation:\n" << endl;
		for (int x = 0; x < size; ++x)
		{
			for (int y = 0; y < size; ++y)
			{
				cout << graph[x][y];
				if (y == size - 1)
					cout << endl;
			}
		}
	};

	void printWeightedGraph()
	{
		cout << "\nMatrix with weights added to edges:\n" << endl;
		for (int x = 0; x < size; ++x)
		{
			for (int y = 0; y < size; ++y)
			{
				cout << weightedGraph[x][y];
				if (y == size - 1)
					cout << endl;
			}
		}
	};

	void runDijktra(int dest)
	{
		this->dest = dest;

		if (isConnected()) //First make sure graph is connected
		{
			while (isDestination == false)
			{
				sourceRow = weightedGraph[currentRow]; //Extract row of edges to next nodes to evaluate
				analyzeEdges(sourceRow); //Populate open set with all of the next possible nodes and their respective edge costs
				selectNode(openSet); //Add next node to open set
			}

			//Print results
			cout << "\n\nDestination: " << dest << "  Least cost to destination: " << totalCost; //Sum up all of edge costs
		}

		else
		{
			cout << "\n\n Graph is not connected. Not every destination can be reached.\n\n" << endl;
		}
	};
};


int main()
{
	Graph myGraph(8, 0.4);
	myGraph.printGraph();
	cout << endl;
	myGraph.addWeights(7);
	myGraph.printWeightedGraph();
	myGraph.runDijktra(7);
	getchar();
}