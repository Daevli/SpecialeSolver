// ---------------------------------
// Author: Johan Arendal Jørgensen
// Title:  Tournament Planning Tool
// Version: 0.3.2
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

ILOSTLBEGIN

using namespace std;

int n; // Number of teams
int m; // Number of rounds in the first half of the tournament
string phaseOneConstraintsPath;
string phaseTwoConstraintsPath;

int main();

// Declare refactoring methods
void printMat(int*, int, int);
void swapRounds(int* mat, int k, int l);
void swapRows(int* mat, int k, int l);
void swapNumbers(int* mat, int k, int l);
void swapTeams(int* mat, int k, int l);
int modMod(int a, int b);
int flipOneAndZero(int t);

int main() {
	// Greeting
	cout << "Program Start...\n" << "How many teams? (Please make it an even number) ";

	// Save number of teams, n, and number of rounds, m
	cin >> n;
	if (n % 2 != 0) { n++; }
	m = n - 1;

	// Dialog: Constraints? (Skip phase 1 or 2?)
	cout << "------\n" << "Number of teams set to: " << n << "\nNumber of rounds: " << 2 * m << "\nTotal number of matches: " << n * m << "\n------\n";
	cout << "Are there any problem-specific constraints on matchups? (1/0)\n"; bool doPhaseOne; cin >> doPhaseOne;
	if (doPhaseOne) {
		cout << "Then we are doing phase 1.\n" << "Please input path/name of constrains file:\n";
		// cin >> phaseOneConstraintsPath; // Save path for constraints file
	}
	else { cout << "Skipping phase 1 then...\n"; }
	cout << "Are there any problem-specific constraints on location? (1/0)\n"; bool doPhaseTwo; cin >> doPhaseTwo;
	if (doPhaseTwo) {
		cout << "Then we are doing phase 2.\n" << "Please input path/name of constrains file:\n";
		// cin >> phaseTwoConstraintsPath; // Save path for constraints file
	}
	else { cout << "Skipping phase 2 then...\n"; }

	auto start = chrono::steady_clock::now();
	/*************************************************************************************************/
	/**********************************| Phase 1: Edge coloring/CP |**********************************/
	/**********************************/ cout << "\n\n--- Phase 1:\n";/*******************************/
	int* M1 = new int[n * m];

	//if (doPhaseOne) {
		// Edge-coloring or latin square?
	//}
	//else {
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
	//}

	// Extends the plan to a full DRR tournament
	int* M1_2 = new int[n * 2 * m];
	int* temp = new int[m];
	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < m; i++) {
			temp[i] = M1[i + m * k];
		}
		for (int j = 0; j < m; j++) {
			M1_2[j + m * (2 * k)] = temp[j];
			M1_2[j + m * (2 * k + 1)] = temp[j];
		}
	}

	printMat(M1_2, n, 2 * m);

	auto postPhaseOne = chrono::steady_clock::now();
	/*************************************************************************************************/
	/*****************************************| Phase 2: IP |*****************************************/
	/*********************************/ cout << "\n\n--- Phase 2:\n";/********************************/
	int* M2 = new int[n * m];

	if (doPhaseTwo) {
		cout << "Problem-specific constraints detected! \nUsing Cplex to solve the IP/CP model...";
		// Build IP/CP and solve with Cplex
		// Build model
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

			// Constraint making sure that the two teams meeting, are playing opposite H/A
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					model.add(h[j * m + i] != h[(M1[j * m + i] - 1) * m + i]);
				}
			}

			// Constraints 3, 4, 5, 6
			for (int i = 0; i < n * m; i += m) {
				e1.clear();
				e1 += h[i] + h[i + 1] - h[i + m - 1];
				model.add(e1 <= 1);
				e1.clear();
				e1 += h[i] + h[i + 1] - h[i + m - 1];
				model.add(e1 >= 0);
				e1.clear();
				e1 += h[i] - h[i + m - 2] - h[i + m - 1];
				model.add(e1 <= 0);
				e1.clear();
				e1 += h[i] - h[i + m - 2] - h[i + m - 1];
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
						cout << "\nHard complementary constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " share stadium.";
						// add constraint to the model. Note, that this ensures that the two teams do not play home at the same time IN THE FIRST HALF.
						for (int i = 0; i < m; i++)	
						{
							e1.clear();
							e1 = h[(v.operator[](1) - 1) * m + i] + h[(v.operator[](2) - 1) * m + i];
							model.add(e1 <= 1);
						}
						break;
					case 2: 
						if (v.operator[](2) == 0) { H = "away"; } else { H = "home"; } // Make string for print
						e1.clear();
						e1 = h[(v.operator[](1) - 1) * m + v.operator[](3) - 1];
						if(v.operator[](3) <= m)	// If the round is in first half, add constraint
						{
							model.add(e1 == v.operator[](2));
							cout << "\nfirst half " << v.operator[](2);
						}
						else {						// Otherwise, add opposite constraint for the first half (remember: tnmt is mirrored)
							model.add(e1 == flipOneAndZero(v.operator[](2)));
							cout << "\nsecond half " << flipOneAndZero(v.operator[](2));
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
			else { cout << "Error while trying to open file!"; }

			// Extract model and set parameters for the solve
			cplex.extract(model);
			cplex.setOut(env.getNullStream()); // <-- Gets rid of cplex output
			// cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 1); // <-- Only need one solution

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

			if (cplex.solve()) {
				cout << "\n\nSolved cplex model Zuccezzfully!";
				// double OBJval = (double)cplex.getObjValue(); // <-- Get obj val
				// cout << "\nNumber of breaks: " << OBJval << endl;
			}
		}
		// end try
		catch (IloException& e) { cerr << "\nConcert exception caught: " << e << endl; }
		catch (...) { cerr << "\nUnknown exception caught" << endl; }
		env.end();
	}
	else {
		cout << "No problem-specific location constraints. \nUsing modified canonical pattern by de Werra (1981)";
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
			if (i % 2 != 0 && i <= m - 4) {
				M2[(n - 1) * m + i] = 1;
				M2[(M1[(n - 1) * m + i] - 1) * m + i] = -1;
			}
			if (i % 2 == 0 && i > m - 4) {
				M2[(n - 1) * m + i] = 1;
				M2[(M1[(n - 1) * m + i] - 1) * m + i] = -1;
			}
			if (i % 2 != 0 && i > m - 4) {
				M2[(n - 1) * m + i] = -1;
				M2[(M1[(n - 1) * m + i] - 1) * m + i] = 1;
			}
		}
	}

	// Extends the plan to a full DRR tournament
	int* M2_2 = new int[n * 2 * m];
	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < m; i++) {
			temp[i] = M2[i + m * k];
		}
		for (int j = 0; j < m; j++) {
			M2_2[j + m * (2 * k)] = temp[j];
			M2_2[j + m * (2 * k + 1)] = -temp[j];
		}
	}

	printMat(M2_2, n, 2 * m);
	auto postPhaseTwo = chrono::steady_clock::now();
	/*************************************************************************************************/
	/***********************************| Phase 3: Metaheuristics |***********************************/
	/*********************************/ cout << "\n\n--- Phase 3:\n";/********************************/
	int* M3 = new int[n * 2 * m];

	// Create M3 by multiplying entries from the solutions found in the previous phases
	for (int i = 0; i < n * 2 * m; i++) {
		M3[i] = M1_2[i] * M2_2[i];
	}
	printMat(M3, n, 2 * m);












	// Print elapsed time
	auto end = chrono::steady_clock::now();
	cout << "\n\nElapsed time: "
		<< endl << "Phase 1: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseOne - start).count() << " ms"
		<< endl << "Phase 2: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseTwo - postPhaseOne).count() << " ms"
		<< endl << "Phase 3: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - postPhaseTwo).count() << " ms\n---"
		<< endl << "Total  : " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"
		<< endl << "------------------------";

	return 0;
}

/*******************************************| Refactoring |*******************************************/

/********* General *********/
// Method for printing out an array as a matrix.
// Args: (Array, # of rows, # of columns)
void printMat(int* arr, int r, int c) {
	cout << endl << endl << "round  |";
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
		cout << endl;
	}
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

/********* Phase 3 *********/
// Method for swapping two rounds in a full tournament plan (move p_1)
// Args: (Matrix, 1st round, 2nd round) - note: uses actual round number (e.g. there is no round 0)
void swapRounds(int* mat, int k, int l) {
	for (int i = 0; i < 2 * n; i++) {
		int temp = mat[i * m + k - 1];
		mat[i * m + k - 1] = mat[i * m + l - 1];
		mat[i * m + l - 1] = temp;
	}
}

// Method for swapping two rows in a full tournament plan
// Args: (Matrix, 1st row, 2nd row)
void swapRows(int* mat, int k, int l) {
	for (int i = 0; i < 2 * m; i++) {
		int temp = mat[i + k * 2 * m];
		mat[i + k * 2 * m] = mat[i + l * 2 * m];
		mat[i + l * 2 * m] = temp;
	}
}

// Method for swapping all entries in a matrix containing k with l.
// Args: (Matrix, 1st int, 2nd int)
void swapNumbers(int* mat, int k, int l) {
	for (int i = 0; i < 2 * n * m; i++) {
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
// Args: (Tournament plan, team 1, team 2) - note: uses actual team number (e.g. there is no team 0)
void swapTeams(int* mat, int k, int l) {
	swapRows(mat, k - 1, l - 1);
	swapNumbers(mat, k, l);
}