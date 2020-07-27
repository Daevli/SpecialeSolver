// ---------------------------------
// Author: Johan Arendal J�rgensen
// Title:  Tournament Planning Tool
// Version: 1.2.4
// ---------------------------------

#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include <conio.h>
#include <fstream>
#include <algorithm>
#include <functional>
#include <queue>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <windows.h>
#include <mmsystem.h>
ILOSTLBEGIN

using namespace std;

int n; // Number of teams.
int m; // Number of rounds in the first half of the tournament.
int topTeamsCaseSix = 0; // Number of top teams declared in the constraint case 6 in phase 1.
int topTeamsCaseSeven = 0; // Same for case 7.
vector<int> bannedRoundsForTopMatches; // List containing the round where a top match is not allowed.
int topTeamsCaseEight = 0; // Same for case 8.
int roundsBetweenTopMatch = 0; // Minimum number of rounds without a top match between two top matches.
vector<int> latinSquare;		// Latin square used for solving the edge coloring problem in phase 1.
vector<long> M1;
vector<long> M2;
vector<long> M3;
vector<vector<int>> costMatrix;

int main();

// Declare refactoring methods
vector<vector<long>> convertVectorOneToTwoDimensions(vector<long> v, int rows, int cols);
void printMat(vector<long> &arr, int n, int m);
int numberOfBreaks(vector<long> plan);
void swapRounds(vector<long> &mat, int k, int l);
void swapRows(vector<long> &mat, int k, int l);
void swapNumbers(vector<long> &mat, int k, int l);
void swapTeams(vector<long>& mat, int k, int l);
void variableNeighborhoodSearch(vector<long> &plan);
void shakeSolution(vector<long>& plan);
void tabuSearch(vector<long>& plan, int tolerence);
long cost(vector<long>& plan);
int modMod(int a, int b);
int flipOneAndZero(int t);
bool firstAcceptVariableNeighborhoodDescent(vector<long>& plan);
bool isNegative(long t);
bool columnsOk(vector<long>& plan);
bool rowsOk(vector<long>& plan);
bool breaksOk(vector<long>& plan);
bool topTeamsOk(vector<long> plan);
bool isFeasible(vector<long>& plan);

int main() {
	// Greeting/Dialog: number of teams?
	cout << "Program Start...\n" << "How many teams? ";

	// Save number of teams, n, and number of rounds, m
	// Add a dummy team to n, if the number given is odd
	cin >> n;
	if (n % 2 != 0) { n++; }
	m = n - 1;

	// Dialog: Any hard constraints? 
	cout << "------\n" << "Number of teams set to: " << n 
		 << "\nNumber of rounds: " << 2 * m << "\nTotal number of matches: " << n * m << "\n------\n";
	cout << "Are there any hard problem-specific constraints on matchups? (1/0)\n"; bool doPhaseOne; 
	cin >> doPhaseOne;
	if (doPhaseOne) {
		cout << "Will search for hard constriants in hardConstraints.txt file.\n";
	}
	else { cout << "No hard constraints. Using circle method...\n"; }
	cout << "Are there any hard problem-specific constraints on location? (1/0)\n"; bool doPhaseTwo; 
	cin >> doPhaseTwo;
	if (doPhaseTwo) {
		cout << "Will search for hard constriants in hardConstraints.txt file.\n";
	}
	else { cout << "No hard constraints. Using modified canonical pattern...\n"; }

	auto start = chrono::steady_clock::now();
	/*****************************************************************************************************************************/
	/*******************************************| Phase 1: Edge-coloring/Latin Square |*******************************************/
	/**********************************************/ cout << "\n\n--- Phase 1:\n";/***********************************************/
	M1.resize(n * m);

	if (doPhaseOne) {
	// Edge-coloring/latin square and solve it with Cplex
		latinSquare.resize(n * n);

		// Build Cplex model
		IloEnv env1;

		// try solving model
		try {
			IloModel model1(env1);
			env1 = model1.getEnv();
			IloCplex cplex1(model1);
			IloNumVarArray latinSquareCplexTemp(env1, n * n);
			IloExpr e1_1(env1);

			cplex1.extract(model1);

			// Add entry variables
			for (int i = 0; i < n * n; i++) {
				latinSquareCplexTemp[i] = IloNumVar(env1, 1, n, ILOINT);
			}

			// Create objective expression
			IloExpr obj1(env1);
			for (int i = 0; i < n * n; i++) {
				obj1 += latinSquareCplexTemp[i];
			}

			// Add objective to the model
			model1.add(IloMinimize(env1, obj1));

			// Constraint such that the LS is symmetric
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++)	{
					e1_1.clear();
					e1_1 += (latinSquareCplexTemp[i * n + j] == latinSquareCplexTemp[j * n + i]);
					model1.add(e1_1 == 1);
				}
			}

			// Constraint such that each team plays "itself" in the first round
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++)	{
					e1_1.clear();
					if (i == j) {
						e1_1 += (latinSquareCplexTemp[i * n + j] == 1);
						model1.add(e1_1 == 1);
					}
				}
			}

			// Constraint such that rows contain all different numbers
			for (int i = 0; i < n; i++)	{
				for (int j = 0; j < n - 1; j++)	{
					e1_1.clear();
					for (int k = j + 1; k < n; k++){
						e1_1 += (latinSquareCplexTemp[i * n + j] == latinSquareCplexTemp[i * n + k]);
					}
					model1.add(e1_1 == 0);
				}
			}

			// Read file 
			string line;
			string H;
			int f;
			int constraintCode;
			int round;
			vector<int> v;
			ifstream myFile("hardConstraints.txt");

			// Problem-specific constraints
			if (myFile.is_open()) {
				std::cout << "\n\nFound constraints: ";
				while (getline(myFile, line)) {			// Read string 
					stringstream ss(line);				// Make stringstream from s 
					v.clear();							// Clear vector v
					while (ss >> f) {					// Read ints from ss into f 
						v.push_back(f);					// Add f to vector 
					}
					constraintCode = v.operator[](0);	// Read first int of v, which determines the type of constraint
					switch (constraintCode) {
					case 3:
						std::cout << "\nHard pairing constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " must NOT play in round " << v.operator[](3);
						round = v.operator[](3) % m;	// If the specified round exceeds number of rounds, uses mod.
						e1_1.clear();
						e1_1 += (latinSquareCplexTemp[(v.operator[](1) - 1) * n + (v.operator[](2) - 1)] != round);
						model1.add(e1_1 == 1);
						break;
					case 5:
						std::cout << "\nHard pairing constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " must play in round " << v.operator[](3);
						round = v.operator[](3) % m;	// If the specified round exceeds number of rounds, uses mod.
						e1_1.clear();
						e1_1 += (latinSquareCplexTemp[(v.operator[](1) - 1) * n + (v.operator[](2) - 1)] == round);
						model1.add(e1_1 == 1);
						break;
					case 6:
						// constraint h5
						std::cout << "\nTop opponent constraint : The first " << v.operator[](1) << " teams are set as top teams in this constraint.";
						if (v.operator[](1) > topTeamsCaseSix) {
							topTeamsCaseSix = v.operator[](1);
						}
						e1_1.clear();
						for (int i = 0; i < n; i++) {
							for (int j = 0; j < v.operator[](1) - 1; j++) {
								for (int k = j + 1; k < v.operator[](1); k++) {
									if (k != i && j != i) {
										e1_1 += (latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] <= 1 
												&& latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] >= -1 
											|| latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] >= n - 2   
											|| latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] <= 2 - n); 
									}
								}							
							}
							model1.add(e1_1 == 0);
							e1_1.clear();
						}
						break;
					case 7: 
						// constraint h6
						std::cout << "\nTop match constraint : No top match is allowed in round " << v.operator[](1) << ". The first " << v.operator[](2) << " teams are set as top teams in this constraint.";
						if (v.operator[](2) > topTeamsCaseSeven) {
							topTeamsCaseSeven = v.operator[](2);
						}
						bannedRoundsForTopMatches.push_back(v.operator[](1));
						e1_1.clear();
						for (int i = 0; i < v.operator[](2); i++) {
							for (int j = 0; j < v.operator[](2); j++) {
								e1_1 += (latinSquareCplexTemp[i * n + j] == v.operator[](1) + 1); // +1 because the latin square has +1. E.g. round 1 in tnmt is 2 in the LatSq
							}
						}
						model1.add(e1_1 == 0);
						break;
					case 8:
						// Constraint h7
						std::cout << "\nTop match constraint : There must be at least " << v.operator[](1) << " rounds between top matches. The first " << v.operator[](2) << " teams are set as top teams in this constraint.";
						if (v.operator[](2) > topTeamsCaseEight) {
							topTeamsCaseEight = v.operator[](2);
						}
						if (v.operator[](1) > roundsBetweenTopMatch) {
							roundsBetweenTopMatch = v.operator[](1);
						}
						e1_1.clear();
						for (int i = 0; i < v.operator[](2); i++) {
							for (int j = 0; j < v.operator[](2) - 1; j++) {
								for (int k = j + 1; k < v.operator[](2); k++) {
									if (k != i && j != i) {
										e1_1 += (latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] <= v.operator[](1) + 1
												&& latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] >= -v.operator[](1) - 1
											|| latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] >= n - v.operator[](1) - 1
											|| latinSquareCplexTemp[j * n + i] - latinSquareCplexTemp[k * n + i] <= v.operator[](1) + 1 - n);
									}
								}
							}
							model1.add(e1_1 == 0);
							e1_1.clear();
						}
						break;
					default:
						break;
					}
				}
				myFile.close();
			}
			else { std::cout << "Error while trying to open file!"; }

			// Extract model and set parameters for the solve
			cplex1.extract(model1);
			// cplex1.setOut(env1.getNullStream()); // <-- Gets rid of cplex output
			cplex1.setParam(IloCplex::Param::MIP::Limits::Solutions, 1); // <-- Only need one solution

			// Solve the model
			cplex1.solve();

			// Save solution to the array latinSquare
			IloNumArray vals1(env1);
			cplex1.getValues(vals1, latinSquareCplexTemp);
			for (int i = 0; i < n * n; i++) {
				latinSquare[i] = vals1[i];
			}

		}
		// end try
		catch (IloException& e) { cerr << "\nConcert exception caught: " << e << "\n"; }
		catch (...) { cerr << "\nUnknown exception caught" << "\n"; }
		env1.end();

		// Print the found latin square
		std::cout << "\nSolution found: \n";
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				std::cout << setw(3) << latinSquare[i * n + j];
			}
			std::cout << "\n";
		}
		// Convert latin square to tournament
		int row;
		int col;
		for (int i = 0; i < latinSquare.size(); i++)
		{
			latinSquare[i]--;
			if (latinSquare[i] != 0) {
				row = floor(i / n);
				col = latinSquare[i] - 1;
				M1[row * m + col] = i % n + 1;
			}
		}
	}
	else {
		// Use circle method (reformulated by Matsui & Miyashiro (2006))
		cout << "No problem-specific matchup constraints. \nUsing the circle method as formulated by Matsui & Miyashiro (2006)";
		for (int i = 0; i < n * m; i++) {
			int t = floor(i / m) + 1;
			int r = (i % m) + 1;
			if (t == n) {
				M1[i] = r;
			}
			else if (t == r) {
				M1[i] = n;
			}
			else {
				M1[i] = modMod((2 * r) - t, n - 1);
			}
		}
	}
	std::cout << "\n\nM1: ";
	printMat(M1, n, m);

	// Extends the plan to a full DRR tournament
	vector<long> M1_2(n * 2 * m);
	M1_2.resize(n * 2 * m);
	vector<long> temp (m);
	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < m; i++) {
			temp[i] = M1[i + m * k];
		}
		for (int j = 0; j < m; j++) {
			M1_2[j + m * (2 * k)] = temp[j];
			M1_2[j + m * (2 * k + 1)] = temp[j];
		}
	}

	// Print solution found in phase 1
	std::cout << "\nM1_2 post phase 1:";
	printMat(M1_2, n, 2 * m);

	// For some reason, this solves a bug concerning the construction of M3 in phase 3
	vector<long> M1_safety(n * 2 * m);
	M1_safety = M1_2;
		
	// Save time of completion
	auto postPhaseOne = chrono::steady_clock::now();

	/*****************************************************************************************************************************/
	/*******************************************************| Phase 2: CP |*******************************************************/
	/***********************************************/ cout << "\n\n--- Phase 2:\n";/**********************************************/
	M2.resize(n * m);	
	vector<long> M2_2(n * 2 * m);
	M2_2.resize(n * 2 * m);
	bool failedHeuristics = false;

	if (!doPhaseTwo) {
		if (!doPhaseOne) {
			cout << "No problem-specific hard constraints in phase 2. \nUsing modified canonical pattern by de Werra (1981) to find the canonical tournament plan\n";
			// Canonical pattern (modified)
			for (int i = 0; i < m; i++) {
				for (int k = 1; k < n - 1; k++) {
					if (k % 2 == 0) {
						M2[((i - k) % (n - 1)) * m + i] = 1;
						M2[((M1[((i - k) % (n - 1)) * m + i] - 1) % (n - 1)) * m + i] = -1;
					}
					else {
						M2[((i + k) % (n - 1)) * m + i] = 1;
						M2[((M1[((i + k) % (n - 1)) * m + i] - 1) % (n - 1)) * m + i] = -1;
					}
				}
				if (i % 2 == 0 && i <= m - 4) {
					M2[(n - 1) * m + i] = -1;
					M2[(M1[(n - 1) * m + i] - 1) * m + i] = 1;
				}
				if (i % 2 == 0 && i > m - 4) {
					M2[(n - 1) * m + i] = 1;
					M2[(M1[(n - 1) * m + i] - 1) * m + i] = -1;
				}
				if (i % 2 != 0 && i <= m - 4) {
					M2[(n - 1) * m + i] = 1;
					M2[(M1[(n - 1) * m + i] - 1) * m + i] = -1;
				}
				if (i % 2 != 0 && i > m - 4) {
					M2[(n - 1) * m + i] = -1;
					M2[(M1[(n - 1) * m + i] - 1) * m + i] = 1;
				}
			}
		}
		else {
			std::cout << "No problem-specific hard constraints concerning location. \nUsing the SD algorithm.";
			vector<vector<long>> tempM2 = convertVectorOneToTwoDimensions(M2_2, n, 2 * m);
			//---------------------------- SD algorithm ------------------------------- 
			// First step: Fill out the first column and the first column in the second half
			for (int i = 0; i < n; i++) {
				if (tempM2[i][0] == 0) {
					tempM2[i][0] = 1;
					tempM2[M1[i * m] - 1][0] = -1;
					tempM2[i][m] = -1;
					tempM2[M1[i * m] - 1][m] = 1;
				}
			}
			// Second step: Find the next empty entry and fill in the negative value of the one to the left...::: (and the opponent and the second half)
			for (int j = 1; j < m; j++) {
				for (int i = 0; i < n; i++) {
					if (tempM2[i][j] == 0) {
						if (tempM2[i][j - 1] == -1) {
							tempM2[i][j] = 1;
							tempM2[M1[i * m + j] - 1][j] = -1;
							tempM2[i][j + m] = -1;
							tempM2[M1[i * m + j] - 1][j + m] = 1;
						}
						if (tempM2[i][j - 1] == 1) {
							tempM2[i][j] = -1;
							tempM2[M1[i * m + j] - 1][j] = 1;
							tempM2[i][j + m] = 1;
							tempM2[M1[i * m + j] - 1][j + m] = -1;
						}
						// Third step: Search for fires
						for (int k = 0; k < n; k++) {
							for (int l = 0; l < 2 * m - 2; l++) {
								// Fourth step: Put out fire and reset indices
								if (tempM2[k][l] + tempM2[k][l + 1] + tempM2[k][l + 2] == 2) {
									for (int h = 0; h < 3; h++) {
										if (tempM2[k][l + h] == 0) {
											tempM2[k][l + h] = -1;
											tempM2[M1[k * m + (l + h) % m] - 1][l + h] = 1;
											tempM2[k][l + h + m] = 1;
											tempM2[M1[k * m + (l + h) % m] - 1][l + h + m] = -1;
										}
									}
									k = 0;
									l = 0;
								}
								if (tempM2[k][l] + tempM2[k][l + 1] + tempM2[k][l + 2] == -2) {
									for (int h = 0; h < 3; h++) {
										if (tempM2[k][l + h] == 0) {
											tempM2[k][l + h] = 1;
											tempM2[M1[k * m + (l + h) % m] - 1][l + h] = -1;
											tempM2[k][l + h + m] = -1;
											tempM2[M1[k * m + (l + h) % m] - 1][l + h + m] = 1;
										}
									}
									k = 0;
									l = 0;
								}
							}
						}
					}
				}
			}
			// Fill in M2_2
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < 2 * m; j++) {
					M2_2[i * 2 * m + j] = tempM2[i][j];
				}
			}
			// A fail-safe mechanism is added, 
			// such that if it produces an infeasible solution, Cplex used to solve instead.
			// Note, that this is unlikely to happen, but if it does, the program will
			// read from the hard constraints file and enforce them
			
			if (!breaksOk(M2_2)) {
				std::cout << "\n\nOops! It seems that the SDMC algorithm didn't work this time! :( \nTrying Mishas algorithm instead!\n";
				// Wipe M2_2 and tempM2
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < 2 * m; j++) {
						tempM2[i][j] = 0;
						M2_2[i * 2 * m + j] = 0;
					}
				}
				//----------------------------------------------- Misha -------------------------------------
				int t = 0;
				for (int j = 0; j < m; j += 2) {
					for (int i = 0; i < n; i++) {
						if (tempM2[i][j] == 0) {
							t = i;
							if (tempM2[i][j + 1] == 0) {
								while (tempM2[t][j] == 0) {
									tempM2[t][j] = 1;
									tempM2[t][j + 1] = -1;
									tempM2[t][j + m] = -1;
									tempM2[t][j + m + 1] = 1;
									t = M1[t * m + j + 1] - 1;
									tempM2[t][j + 1] = 1;
									tempM2[t][j] = -1;
									tempM2[t][j + m + 1] = -1;
									tempM2[t][j + m] = 1;
									t = M1[t * m + j] - 1;
								}
							}
							// the last column is a special case
							else if (isNegative(tempM2[t][j - 1]) && isNegative(tempM2[t][j + 1])) {
								tempM2[t][j] = 1;
								tempM2[t][j + m] = -1;
								t = M1[t * m + j] - 1;
								tempM2[t][j] = -1;
								tempM2[t][j + m] = 1;
							}
							else if (!isNegative(tempM2[t][j - 1]) && !isNegative(tempM2[t][j + 1])) {
								tempM2[t][j] = -1;
								tempM2[t][j + m] = 1;
								t = M1[t * m + j] - 1;
								tempM2[t][j] = 1;
								tempM2[t][j + m] = -1;
							}
							else if (isNegative(tempM2[M1[t * m + j] - 1][j - 1]) && isNegative(tempM2[M1[t * m + j] - 1][j + 1])) {
								tempM2[t][j] = -1;
								tempM2[t][j + m] = 1;
								t = M1[t * m + j] - 1;
								tempM2[t][j] = 1;
								tempM2[t][j + m] = -1;
							}
							else if (!isNegative(tempM2[M1[t * m + j] - 1][j - 1]) && !isNegative(tempM2[M1[t * m + j] - 1][j + 1])) {
								tempM2[t][j] = 1;
								tempM2[t][j + m] = -1;
								t = M1[t * m + j] - 1;
								tempM2[t][j] = -1;
								tempM2[t][j + m] = 1;
							}
							else {
								tempM2[t][j] = -1;
								tempM2[t][j + m] = 1;
								t = M1[t * m + j] - 1;
								tempM2[t][j] = 1;
								tempM2[t][j + m] = -1;
							}
						}
					}
				}
				//---------------------------------------------------------------------
				// Fill in M2_2
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < 2 * m; j++) {
						M2_2[i * 2 * m + j] = tempM2[i][j];
					}
				}
			}
			if (!breaksOk(M2_2)) {
				std::cout << "\nYikes! That didn't work either! \nLast resort: Cplex... \n(Note that the hardConstraints.txt file will be read and hard constraints may be enforced).\n";
				failedHeuristics = true;
			}
		}
		
	}
	if (doPhaseTwo || failedHeuristics) {
		std::cout << "Problem-specific constraints detected! \nUsing Cplex to solve the CP model...";

		// Build model (CP model used in phase 2)
		IloEnv env;

		try {
			IloModel model(env);
			env = model.getEnv();
			IloCplex cplex(model);
			IloNumVarArray h(env, n * m);
			IloNumVarArray b(env, n * m);
			IloExpr e1(env);

			cplex.extract(model);

			// Add variables h_ir and b_ir
			for (int i = 0; i < n * m; i++) {
				h[i] = IloNumVar(env, 0, 1, ILOINT);
				b[i] = IloNumVar(env, 0, 1, ILOINT);
			}

			// Create objective expression
			IloExpr obj(env);
			for (int i = 0; i < n * m; i += m) {
				obj += b[i];
			}
			for (int i = 0; i < n * m; i++) {
				if (i % m != 0) {
					obj += 2 * b[i];
				}
			}

			// Add objective to the model
			model.add(IloMinimize(env, obj));

			// Constraints 1, 2
			for (int i = 0; i < n * m - 2; i++) {
				e1.clear();
				e1 += h[i] + h[i + 1] + h[i + 2];
				model.add(e1 <= 2); // Constraint 1
				e1.clear();
				e1 += h[i] + h[i + 1] + h[i + 2];
				model.add(e1 >= 1); // Constraint 2
			}

			// Constraint ensuring that all rows in the tournament are different
			for (int i = 0; i < n - 1; i++) {
				for (int k = 1; k < n - i; k++) {
					e1.clear();
					for (int j = 0; j < m; j++) {
						e1 += (h[j + i * m] != h[j + (i + k) * m]);
					}
					model.add(e1 >= 1);
				}
			}

			// Constraint making sure that the two teams playing, are playing opposite H/A
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					model.add(h[j * m + i] != h[(M1[j * m + i] - 1) * m + i]);
				}
			}

			// Constraints 3, 4, 5, 6
			for (int i = 0; i < n; i++) {
				e1.clear();
				e1 += h[i * m] + h[i * m + 1] - h[i * m + m - 1];
				model.add(e1 <= 1);
				e1.clear();
				e1 += h[i * m] + h[i * m + 1] - h[i * m + m - 1];
				model.add(e1 >= 0);
				e1.clear();
				e1 += h[i * m] - h[i * m + m - 2] - h[i * m + m - 1];
				model.add(e1 <= 0);
				e1.clear();
				e1 += h[i * m] - h[i * m + m - 2] - h[i * m + m - 1];
				model.add(e1 >= -1);
			}

			// CP equivalent of constraints 7, 8, 9 and 10
			for (int i = 1; i < n * m; i++) {
				model.add(IloIfThen(env, h[i - 1] == h[i], b[i] >= 1));
			}

			// Constraint 11
			for (int i = 0; i < m; i++) {
				e1.clear();
				for (int j = 0; j < n; j++) {
					e1 += h[i + j * m];
				}
				model.add(e1 == n / 2);
			}

			// Constraints 12, 13
			for (int i = 0; i < n; i++) {
				e1.clear();
				for (int j = 0; j < m; j++) {
					e1 += h[i * m + j];
				}
				model.add(e1 <= m - floor(m / 3)); // Constraint 12
				model.add(e1 >= floor(m / 3)); // Constraint 13
			}

			string line;
			string H;
			int f;
			int constraintCode;
			vector<int> v;
			ifstream myFile("hardConstraints.txt");

			// Problem-specific constraints
			if (myFile.is_open()) { 
				cout << "\n\nFound constraints: ";
				while (getline(myFile, line)) {			// Read string 
					stringstream ss(line);				// Make stringstream from s 
					v.clear();							// Clear vector v
					while (ss >> f) {					// Read ints from ss into f 
						v.push_back(f);					// Add f to vector 
					}
					constraintCode = v.operator[](0);	// Read first int of v, which determines the type of constraint
					switch (constraintCode)	{
					case 1: 
						// add constraint to the model. Note, that this ensures that the two teams do not play home at the same time IN THE FIRST HALF.
						for (int i = 0; i < m; i++)	
						{
							e1.clear();
							e1 = h[(v.operator[](1) - 1) * m + i] + h[(v.operator[](2) - 1) * m + i];
							model.add(e1 <= 1);
						}
						cout << "\nHard complementary constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " share stadium.";
						break;
					case 2: 
						if (v.operator[](2) == 0) { H = "away"; } else { H = "home"; } // Make string for print
						e1.clear();
						e1 = h[(v.operator[](1) - 1) * m + v.operator[](3) - 1];
						if(v.operator[](3) <= m) {	
							model.add(e1 == v.operator[](2));	// If the round is in first half, add constraint
						}
						else {						
							model.add(e1 == flipOneAndZero(v.operator[](2)));	// Otherwise, add opposite constraint for the first half (remember: tnmt is mirrored)
						}
						cout << "\nHard availability constraint  : Team " << v.operator[](1) << " must play " << H << " in round " << v.operator[](3) << ".";
						break;
					case 4:
						cout << "\nTriple condition              : At most two of teams " << v.operator[](1) << ", " << v.operator[](2) << " and " << v.operator[](3) << " may play home in any round.";
						// add constraint to the model
						for (int i = 0; i < m; i++)	{
							e1.clear();
							e1 = h[(v.operator[](1) - 1) * m + i] + h[(v.operator[](2) - 1) * m + i] + h[(v.operator[](3) - 1) * m + i];
							model.add(e1 <= 2);
							model.add(e1 >= 1);
						}
					default:
						break;
					}
				}
				myFile.close();
			}
			else { std::cout << "Error while trying to open file!"; }
			std::cout << "\n---\n";
			// Extract model and set parameters for the solve
			cplex.extract(model);
			cplex.setOut(env.getNullStream()); // <-- Gets rid of cplex output
			cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 1); // <-- Only need one solution
			// Solve the model
			cplex.solve();


			// Save H/A pattern (Solution) to the array M2
			IloNumArray vals(env);
			cplex.getValues(vals, h);
			for (int i = 0; i < n * m; i++) {
				M2[i] = vals[i];
				if (M2[i] == 0) {
					M2[i] = -1;
				}
			}

		}
		// end try
		catch (IloException& e) { cerr << "\nConcert exception caught: " << e << "\n"; }
		catch (...) { cerr << "\nUnknown exception caught" << "\n"; }
		env.end();
	}

	if (doPhaseTwo || !doPhaseOne || failedHeuristics) {
		// Extends the plan to a full DRR tournament
		for (int k = 0; k < n; ++k) {
			for (int i = 0; i < m; i++) {
				temp[i] = M2[i + m * k];
			}
			for (int j = 0; j < m; j++) {
				M2_2[j + m * (2 * k)] = temp[j];
				M2_2[j + m * (2 * k + 1)] = -temp[j];
			}
		}
	}

	std::cout << "\nDone!\n\nM2_2 post phase 2:";
	printMat(M2_2, n, 2 * m);
	auto postPhaseTwo = chrono::steady_clock::now();
	/*****************************************************************************************************************************/
	/*************************************************| Phase 3: Metaheuristics |*************************************************/
	/***********************************************/ cout << "\n\n--- Phase 3:\n";/**********************************************/
	M3.resize(n * 2 * m);

	// Create M3 by multiplying entries from the solutions found in the previous phases
	for (int i = 0; i < n * 2 * m; i++) {
		M3[i] = M1_2[i] * M2_2[i];
	}

	// Print M3 before solving
	cout << "\nM3: feasible: " << isFeasible(M3);
	printMat(M3, n, 2 * m);

	// Read costs file 
	string line;
	int f;
	int constraintCode;
	int round;
	long team1;
	long team2;
	vector<int> v;
	ifstream myFile("costs.txt");
	if (myFile.is_open()) {
		std::cout << "\n\nFound costs: ";
		while (getline(myFile, line)) {			// Read string
			stringstream ss(line);				// Make stringstream from s
			v.clear();							// Clear vector v
			while (ss >> f) {					// Read ints from ss into f
				v.push_back(f);					// Add these to vector
			}
			switch (v[0]) {
			case 1: 
				std::cout << "\nSoft constraint s1: A cost of " << setw(9) << v[3] << " if team " << v[1] << " and " << v[2] << " play home at the same time.";
				break;
			case 2: 
				std::cout << "\nSoft constraint s2: A cost of " << setw(9) << v[3] << " if team " << v[1] << " plays home in round " << v[2];
				break;
			case 3: 
				std::cout << "\nSoft constraint s3: A cost of " << setw(9) << v[4] << " if team " << v[1] << " and team " << v[2] << " play in round " << v[3];
				break;
			default:
				break;
			}
			costMatrix.push_back(v);
		}
		myFile.close();
		std::cout << "\nIn addition, costs are added for each break in the tournament and soft constraint s4.";
	}
	else { std::cout << "Error while trying to open file!"; }


	variableNeighborhoodSearch(M3);
	//tabuSearch(M3, 20);

	// Print out some information on the solution
	std::cout	<< "\n-------------------------"
				<< "\nConstraints in phase 1: " << doPhaseOne 
				<< "\nConstraints in phase 2: " << doPhaseTwo
				<< "\nSolution feasible:      " << isFeasible(M3)
				<< "\nNumber of breaks:       " << numberOfBreaks(M3);
	
	// Print elapsed time
	auto end = chrono::steady_clock::now();
	cout << "\n---\nElapsed time: "
		<< "\n" << "Phase 1: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseOne - start).count() << " ms"
		<< "\n" << "Phase 2: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseTwo - postPhaseOne).count() << " ms"
		<< "\n" << "Phase 3: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - postPhaseTwo).count() << " ms\n---"
		<< "\n" << "Total  : " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"
		<< "\n" << "------------------------";

	// Play sound on completion
	std::string a1 = "av";
	PlaySoundW(LPCWSTR(a1.c_str()), NULL, 0x0000);
	return 0;
}

/*********************************************************| Refactoring |*********************************************************/

/********* General *********/
// Method for printing out an array as a matrix.
// Args: (Array, # of rows, # of columns)
void printMat(vector<long> &arr, int r, int c) {
	cout << "\n" << "\n" << "round  |";
	for (int j = 0; j < c; j++) {
		cout << setw(4) << j + 1;
	}
	cout << "\n_______|";
	for (int j = 0; j < c; j++) {
		cout << setw(4) << "____";
	}
	cout << "\n";
	for (int i = 0; i < r; i++) {
		cout << "team " << setw(2) << i + 1 << "| ";
		for (int j = 0; j < c; j++) {
			cout << setw(3) << arr[i * c + j] << ' ';
		}
		cout << "\n";
	}
}

// Method that returns number of breaks for a full tournament
int numberOfBreaks(vector<long> plan) {
	int num = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2 * m - 1; j++) {
			if (isNegative(plan[i * 2 * m + j]) == isNegative(plan[i * 2 * m + j + 1])) {
				num += 1;
			}
		}
	}
	return num;
}

/********* Phase 1 *********/


/********* Phase 2 *********/
// The modified modulo function used by Matsui & Miyashiro (2006)
int modMod(int a, int b) {
	int t;
	if (a % b > 0) {
		t = a % b;
	}
	else if (a % b < 0) {
		t = a % b + m; // You might think that this is not necessary, but for some reason, modulo can be negative here...
	}
	else {
		t = b;
	}
	return t;
}

// Method that returns 1 if the argument is 0 and returns 0 otherwise
int flipOneAndZero(int t) {
	if (t == 0) {
		return 1; 
	}
	if (t == 1){
		return 0;
	}
	return 0;
}

// Returns a 2-dimensional version of the given 1-dimensional vector
vector<vector<long>> convertVectorOneToTwoDimensions(vector<long> v, int rows, int cols) {
	vector<vector<long>> temp;
	temp.resize(rows);
	for (int i = 0; i < rows; i++) {
		temp[i].resize(cols);
		for (int j = 0; j < cols; j++) {
			temp[i][j] = v[i * cols + j];
		}
	}
	return temp;
}


/********* Phase 3 *********/
// Method for swapping two rounds in a full tournament plan (move p_1)
// Args: (Matrix, 1st round, 2nd round)
void swapRounds(vector<long>& mat, int k, int l) {
	for (int i = 0; i < 2 * n; i++) {
		swap(mat[i * m + k], mat[i * m + l]);
	}
}

// Method for swapping two rows in a full tournament plan
// Args: (Matrix, 1st row, 2nd row)
void swapRows(vector<long> &mat, int k, int l) {
	for (int i = 0; i < 2 * m; i++) {
		swap(mat[i + k * 2 * m], mat[i + l * 2 * m]);
	}
}

// Method for swapping all entries in a matrix containing k with l.
// Args: (Matrix, 1st int, 2nd int)
void swapNumbers(vector<long> &mat, int k, int l) {
	for (int i = 0; i < mat.size(); i++) {
		if (mat[i] == k) {
			mat[i] = l;
		}
		else if (mat[i] == -k) {
			mat[i] = -l;
		}
		else if (mat[i] == l) {
			mat[i] = k;
		}
		else if (mat[i] == -l) {
			mat[i] = -k;
		}
	}
}

// Method for swapping the schedules of two teams (move p_2)
// Args: (Tournament plan, team 1, team 2) 
void swapTeams(vector<long> &mat, int k, int l) {
	swapRows(mat, k, l);
	swapNumbers(mat, k+1, l+1);
}

// Returns the costs for a tournament plan as defined by a .txt file
long cost(vector<long>& plan) {
	long cost = 0;
	int constraintCode = 0;
	int round = 0;
	long team1 = 0;
	long team2 = 0;

	for (int i = 0; i < costMatrix.size(); i++) {
		constraintCode = costMatrix[i][0];
		switch (constraintCode) {
		case 1:	// Soft constraint s1
			team1 = costMatrix[i][1] - 1;
			team2 = costMatrix[i][2] - 1;
			for (int j = 0; j < 2 * m; j++) {
				if (isNegative(plan[team1 * 2 * m + j]) == isNegative(plan[team2 * n + j])
					&& !isNegative(plan[team1 * n + j])) {
					cost += costMatrix[i][3];
				}
			}
			break;
		case 2: // Soft constraint s2
			team1 = costMatrix[i][1] - 1;
			round = costMatrix[i][2] - 1;
			if (!isNegative(plan[team1 * 2 * m + round])) {
				cost += costMatrix[i][3];
			}
			break;
		case 3: // Soft constraint s3
			team1 = costMatrix[i][1] - 1;
			team2 = costMatrix[i][2] - 1;
			round = costMatrix[i][3] - 1;
			if (plan[team1 * 2 * m + round] == team2 || plan[team1 * 2 * m + round] == -team2) {
				cost += costMatrix[i][4];
			}
			break;
		default:
			break;
		}
	}

	// Costs for each break in the plan.
	cost += numberOfBreaks(plan) * 9000000;

	// Soft constraint s4 should always be in place
	long ms;
	long smallest;
	for (int col = 0; col < 2 * m; col++) {
		smallest = 10000;
		for (int row = 0; row < n; row++) {
			ms = (plan[row * 2 * m + col]) * (row + 1);
			if (isNegative(ms)) {
				ms = -ms;
			}
			if (ms < smallest) {
				smallest = ms;
			}
		}
		cost += smallest;
	}
	return cost;
}

// Returns the first improving neighbor (first accept) to a given tournament plan
// The neighborhood used is determined by a given integer, num
bool firstAcceptVariableNeighborhoodDescent(vector<long> &plan) {
	vector<long> oldPlan = plan; // Save the plan before making changes
	int oldCost = cost(plan); // Save the cost of the original plan
	int num = 0;
	while (num < 4) {
		// move p1
		cout << "\nmove p1: Swap rounds  : ";
		for (int i = 0; i < m - 1; i++)
			for (int j = i + 1; j < m; j++) {
				swapRounds(plan, i, j);
				if (isFeasible(plan) && oldCost > cost(plan)) {			// If the new plan is better, return it
					cout << "FA: found better";
					return true;
				}
				plan = oldPlan;
			}
		cout << "no improvement, trying another neighborhood... ";
		num++;
		// move p2
		cout << "\nmove p2: Swap teams   : ";
		for (int k = 0; k < n - 1; k++)
			for (int l = k + 1; l < n; l++) {
				swapTeams(plan, k, l);
				if (isFeasible(plan) && oldCost > cost(plan)) {			// If the new plan is better, return it
					cout << "FA: found better";
					return true;
				}
				plan = oldPlan;
			}
		cout << "no improvement, trying another neighborhood... ";
		num++;
		// move p1-p2
		cout << "\nmove p1-p2            : ";
		for (int i = 0; i < m - 1; i++)
			for (int j = i + 1; j < m; j++)
				for (int k = 0; k < n - 1; k++)
					for (int l = k + 1; l < n; l++) {					
						swapRounds(plan, i, j);
						swapTeams(plan, k, l);
						if (isFeasible(plan) && oldCost > cost(plan)) {			// If the new plan is better, return it
							cout << "FA: found better";
							return true;
						}
						plan = oldPlan;
					}
		cout << "no improvement, trying another neighborhood... ";
		num++;
		// move p2-p1
		cout << "\nmove p2-p1            : ";
		for (int i = 0; i < m - 1; i++)
			for (int j = i + 1; j < m; j++)
				for (int k = 0; k < n - 1; k++)
					for (int l = k + 1; l < n; l++) {					
						swapTeams(plan, k, l);
						swapRounds(plan, i, j);
						if (isFeasible(plan) && oldCost > cost(plan)) {			// If the new plan is better, return it
							cout << "FA: found better";
							return true;
						}
						plan = oldPlan;
					}
		cout << "no improvement. ";
		num++;	
	}

	return false; // If we go through all neighbors without improvement, the old plan is returned.
}

void shakeSolution(vector<long>& plan) {
	vector<long> oldPlan = plan; // Save the plan before making changes
	int num = 0;
	int k = 0;
	int l = 0;
	int i = 0;
	int j = 0;
	while (num < 100000) {
		// choose two random teams k, l and two random rounds i, j
		k = rand() % n;
		l = rand() % n;
		i = rand() % m;
		j = rand() % m;

		// Try to do move p1 and/or p2 with these values.
		swapTeams(plan, k, l);
		if (isFeasible(plan)) {			
			oldPlan = plan;
		}
		else {
			plan = oldPlan;
		}
		swapRounds(plan, i, j);
		if (isFeasible(plan)) {
			oldPlan = plan;
		}
		else {
			plan = oldPlan;
		}
		num++;
	}
}

// Uses first-accept in a variable neighborhood search
void variableNeighborhoodSearch(vector<long>& plan) {
	std::cout << "\n\n=== Variable neighborhood search\nInitial solution:";
	printMat(plan, n, 2 * m);
	std::cout << "\nWith a cost of: " << cost(plan) << "\n";
	vector<long> bestPlan = plan;
	int iteration = 0;
	int numberOfShakes = 0;
	while (numberOfShakes < 5) {
		while (firstAcceptVariableNeighborhoodDescent(plan)) {
			iteration++;
		}
		if (cost(plan) < cost(bestPlan)) {
			std::cout << "\n\nFound better plan by shaking!! Diff: " << cost(plan) - cost(bestPlan) << "\n\n ";
			bestPlan = plan;
		}
		std::cout << "\n\nShaking!\n";
		shakeSolution(plan);
		numberOfShakes++;
	}

	cout << "\n\nDone!\nNew solution:";
	printMat(bestPlan, n, 2 * m);
	cout << "\n\nWith a cost of: " << cost(bestPlan);
}

// Tabu search
void tabuSearch(vector<long>& plan, int tolerence) {
	cout << "\n\n=== Tabu search === \nInitial solution:";
	printMat(plan, n, 2 * m);
	cout << "\nWith a cost of: " << cost(plan) << "\n" << "Doing tabu search";
	int fail = 0;
	int iteration = 0;
	int a = 100000;
	int b = 1000;
	int c = 10000;
	int d = 100000;
	int t = floor(n / 4);
	vector<int> teamTabuList;
	vector<int> roundTabuList;
	vector<long> originalSolution = plan;
	vector<long> tempNabo = plan;
	vector<long> bestNabo = plan;
	vector<long> optimalSolution = plan; // (or near-optimal)

	while (fail < tolerence) {
		iteration++;
		a = 100000;
		b = 100000;
		c = 100000;
		d = 100000;
		// Reset initial solution to the best solution from the last iteration
		plan = bestNabo;
		tempNabo = bestNabo;

		// Reset the best neighbor to the original solution
		bestNabo = originalSolution;

		int switcheroo = 0;
		switcheroo = iteration % 4;
		switch (switcheroo) {
		case 0: // move p1: switch rounds
			for (int i = 0; i < m - 1; i++) {
				if (std::find(roundTabuList.begin(), roundTabuList.end(), i) == roundTabuList.end()) {
					for (int j = i + 1; j < m; j++)	{
						if (std::find(roundTabuList.begin(), roundTabuList.end(), j) == roundTabuList.end()) {
							swapRounds(tempNabo, i, j);
							if (cost(tempNabo) < cost(bestNabo) && isFeasible(tempNabo)) {
								bestNabo = tempNabo;
								a = i;
								b = j;
							}
							swapRounds(tempNabo, i, j);
						}
					}
				}
			}
			// Update tabuList
			if (std::find(roundTabuList.begin(), roundTabuList.end(), a) != roundTabuList.end()) {
				roundTabuList.push_back(a); // If a is not in the tabulist, add it
			}
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(roundTabuList.begin(), roundTabuList.end(), b) != roundTabuList.end()) {
				roundTabuList.push_back(b); // Same for b
			}			
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			break;
		case 1: // move p2: switch teams
			for (int k = 0; k < n - 1; k++)	{
				if (std::find(teamTabuList.begin(), teamTabuList.end(), k) == teamTabuList.end()) {
					for (int l = k + 1; l < n; l++) {
						if (std::find(teamTabuList.begin(), teamTabuList.end(), l) == teamTabuList.end()) {
							swapTeams(tempNabo, k, l); // Swap two teams, since neither are on the tabulist
							if (cost(tempNabo) < cost(bestNabo) && isFeasible(tempNabo)) {
								bestNabo = tempNabo; // If the cost is better than the best neighbor so far, save it as the new best neighbor
								c = k; // save the two teams that were swapped
								d = l;
							}
							swapTeams(tempNabo, k, l); // swap back, so we are in the same neighborhood
						}
					}
				}
			}
			// Update tabuList
			if (std::find(teamTabuList.begin(), teamTabuList.end(), c) != teamTabuList.end()) {
				teamTabuList.push_back(c); // If c is not in the tabulist, add it
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(teamTabuList.begin(), teamTabuList.end(), d) != teamTabuList.end()) {
				teamTabuList.push_back(d); // Same for d
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			break;
		case 2: // move p1-p2
			for (int i = 0; i < m - 1; i++) {
				if (std::find(roundTabuList.begin(), roundTabuList.end(), i) == roundTabuList.end()) {
					for (int j = i + 1; j < m; j++) {
						if (std::find(roundTabuList.begin(), roundTabuList.end(), j) == roundTabuList.end()) {
							for (int k = 0; k < n - 1; k++) {
								if (std::find(teamTabuList.begin(), teamTabuList.end(), k) == teamTabuList.end()) {
									for (int l = k + 1; l < n; l++) {
										if (std::find(teamTabuList.begin(), teamTabuList.end(), l) == teamTabuList.end()) {
											swapRounds(tempNabo, i, j);
											swapTeams(tempNabo, k, l); // Swap two teams, since neither are on the tabulist
											if (cost(tempNabo) < cost(bestNabo) && isFeasible(tempNabo)) {
												bestNabo = tempNabo; // If the cost is better than the best neighbor so far, save it as the new best neighbor
												a = i; 
												b = j;
												c = k;
												d = l;
											}
											swapTeams(tempNabo, k, l); // swap back, so we are in the same neighborhood
											swapRounds(tempNabo, i, j);
										}
									}
								}
							}
						}
					}
				}
			}
			// Update tabuList
			if (std::find(roundTabuList.begin(), roundTabuList.end(), a) != roundTabuList.end()) {
				roundTabuList.push_back(a); // If a is not in the tabulist, add it
			}
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(roundTabuList.begin(), roundTabuList.end(), b) != roundTabuList.end()) {
				roundTabuList.push_back(b); // Same for b
			}
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}			
			if (std::find(teamTabuList.begin(), teamTabuList.end(), c) != teamTabuList.end()) {
				teamTabuList.push_back(c); // If a is not in the tabulist, add it
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(teamTabuList.begin(), teamTabuList.end(), d) != teamTabuList.end()) {
				teamTabuList.push_back(d); // Same for b
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			break;
		case 3: // move p2-p1
			for (int i = 0; i < m - 1; i++) {
				if (std::find(roundTabuList.begin(), roundTabuList.end(), i) == roundTabuList.end()) {
					for (int j = i + 1; j < m; j++) {
						if (std::find(roundTabuList.begin(), roundTabuList.end(), j) == roundTabuList.end()) {
							for (int k = 0; k < n - 1; k++) {
								if (std::find(teamTabuList.begin(), teamTabuList.end(), k) == teamTabuList.end()) {
									for (int l = k + 1; l < n; l++) {
										if (std::find(teamTabuList.begin(), teamTabuList.end(), l) == teamTabuList.end()) {
											swapTeams(tempNabo, k, l);   // Swap two teams, since neither are on the tabulist
											swapRounds(tempNabo, i, j);
											if (cost(tempNabo) < cost(bestNabo) && isFeasible(tempNabo)) {
												bestNabo = tempNabo; // If the cost is better than the best neighbor so far, save it as the new best neighbor
												a = i;
												b = j;
												c = k;
												d = l;
											}
											swapRounds(tempNabo, i, j);
											swapTeams(tempNabo, k, l); // swap back, so we are in the same neighborhood
										}
									}
								}
							}
						}
					}
				}
			}
			// Update tabuList
			if (std::find(roundTabuList.begin(), roundTabuList.end(), a) != roundTabuList.end()) {
				roundTabuList.push_back(a); // If a is not in the tabulist, add it
			}
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(roundTabuList.begin(), roundTabuList.end(), b) != roundTabuList.end()) {
				roundTabuList.push_back(b); // Same for b
			}
			if (roundTabuList.size() > t) {
				roundTabuList.erase(roundTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(teamTabuList.begin(), teamTabuList.end(), c) != teamTabuList.end()) {
				teamTabuList.push_back(c); // If a is not in the tabulist, add it
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			if (std::find(teamTabuList.begin(), teamTabuList.end(), d) != teamTabuList.end()) {
				teamTabuList.push_back(d); // Same for b
			}
			if (teamTabuList.size() > t) {
				teamTabuList.erase(teamTabuList.begin()); // If the tabuList is too big, remove the first element
			}
			break;
		default:
			break;
		}

		// Update best solution and fail counter
		if (cost(optimalSolution) <= cost(bestNabo)) {
			fail++;
		}
		else {
			fail = 0;
			optimalSolution = bestNabo;
		}
		if (fail % 31 == 0 || fail == tolerence) {
			std::cout << "\rOn iteration " << fail << "/" << tolerence << " with no improvement.";
			std::cout.flush();
		}
	}
	// Finally, update the tournament plan
	plan = optimalSolution;
	cout << " Tabu search done!\nNew solution:";
	printMat(plan, n, 2 * m);
	cout << "\n\nWith a cost of: " << cost(plan)
		<< "\n---\nTabu search parameters:\nTolerence for no improvement: " << tolerence
		<< "\nSize of tabulist(s):          " << t
		<< "\nNumber of iterations:         " << iteration
		<< "\n---";
}

// Returns true if the argument is negative
bool isNegative(long t) {
	return (t < 0);
}

// Returns true, if every column in a tournament plan has values from 1, ..., n
bool columnsOk(vector<long>& plan) {
	for (int i = 0; i < 2 * m; i++)	{
		for (int j = 0; j < n - 1; j++) 
		for (int k = j + 1; k < n; k++) {
			if (plan[j * 2 * m + i] == plan[k * 2 * m + i]) {
				return false;
			}
		}
	}
	return true;
}

// Returns false, if there is a problem with the row concerning feasibility
bool rowsOk(vector<long>& plan) {
	// Check if a team plays the same team multiple times in the same half
	for (int i = 0; i < 2 * n; i++) 
		for (int j = 0; j < m - 1; j++)
			for (int k = j + 1; k < m; k++) {
				if (plan[i * m + j] == plan[i * m + k]) {
					return false;
				}
			}
	// Check if a team plays itself
	
	//cout << " Plays itself? ";
	for (int i = 0; i < n; i++)
		for (int j = 0; j < 2 * m; j++) {
			//cout << "\nIs " << plan[i * 2 * m + j] << " or " << -plan[i * 2 * m + j] << " equal to " << i+1 << "?";
			if (plan[i * 2 * m + j] == i + 1 || -plan[i * 2 * m + j] == i + 1) {
				//cout << " Team plays itself!";										//<--- Bug: If there is a cout somewhere in here, it works fine, otherwise not........
				return false;
			}
		}
	return true;

}

// Returns false, if there are 2 breaks (3 elements with the same sign) in a row
bool breaksOk(vector<long>& plan) {
	for (int i = 0; i < n; i++)
	for (int j = 0; j < 2 * m - 2; j++) {
		// Return false, if the j'th element in a row has the same sign as the next two elements
		if (isNegative(plan[i * 2 * m + j]) == isNegative(plan[i * 2 * m + j + 1]) 
			&& isNegative(plan[i * 2 * m + j + 1]) == isNegative(plan[i * 2 * m + j + 2])) {
			return false;
		}
	}
	return true;
}

// Returns false, if the top teams constraints in phase 1 are broken.
bool topTeamsOk(vector<long> plan) {
	int first;
	int second;
	// Case 6 in the switch statement (constraint h5)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2 * m - 1; j++) {
			first = plan[i * 2 * m + j];
			second = plan[i * 2 * m + j + 1];
			if (isNegative(first)) {
				first = -first;
			}
			if (isNegative(second)) {
				second = -second;
			}
			if (first <= topTeamsCaseSix && second <= topTeamsCaseSix) {
				return false;
			}
		}
	}
	// Case 7 (constraint h6)
	for (int k = 0; k < bannedRoundsForTopMatches.size(); k++) {
		int j = bannedRoundsForTopMatches[k] - 1;
		for (int i = 0; i < topTeamsCaseSeven; i++) {
			if (plan[i * 2 * m + j] <= topTeamsCaseSeven && plan[i * 2 * m + j] >= -topTeamsCaseSeven) {
				return false;
			}
		}
	}
	// Case 8 (constraint h7)
	for (int i = 0; i < topTeamsCaseEight; i++) {
		for (int j = 0; j < 2 * m - roundsBetweenTopMatch; j++) {
			first = plan[i * 2 * m + j];
			for (int k = 1; k <= roundsBetweenTopMatch; k++)	{
				second = plan[i * 2 * m + j + k];
				if (isNegative(first)) {
				first = -first;
			}
				if (isNegative(second)) {
				second = -second;
			}
				if (first <= topTeamsCaseEight && second <= topTeamsCaseEight) {
					return false;
				}
			}
		}
	}
	return true;
}

// Returns true if the tournament plan is feasible
bool isFeasible(vector<long>& plan) {
	return (columnsOk(plan) && rowsOk(plan) && breaksOk(plan) && topTeamsOk(plan));
}
