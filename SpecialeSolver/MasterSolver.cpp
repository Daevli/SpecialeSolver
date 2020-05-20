// ---------------------------------
// Author: Johan Arendal J�rgensen
// Title:  Tournament Planning Tool
// Version: 0.3.6
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
#include <Windows.h>
ILOSTLBEGIN

using namespace std;

int n; // Number of teams
int m; // Number of rounds in the first half of the tournament
vector<int> latinSquare;
vector<int> fixedMatchups;
vector<int> bannedRounds;

int main();

// Declare refactoring methods
void printMat(vector<long>, int, int);
void swapRounds(vector<long> mat, int k, int l);
void swapRows(vector<long> mat, int k, int l);
void swapNumbers(vector<long> mat, int k, int l);
void swapTeams(vector<long> mat, int k, int l);

int modMod(int a, int b);
int flipOneAndZero(int t);

bool usedInRow(vector<int> latSq, int row, int num);
bool usedInCol(vector<int> latSq, int col, int num);
bool isSafe(vector<int> latSq, int row, int col, int num);
bool solveLatinSquare(vector<int> latSq);

pair<int, int> getUnassignedLocation(vector<int> latSq);

int main() {
	// Greeting
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
		cout << "Reading from hard constraints file...\n";
	}
	else { cout << "Using circle method then...\n"; }
	cout << "Are there any hard problem-specific constraints on location? (1/0)\n"; bool doPhaseTwo; 
	cin >> doPhaseTwo;
	if (doPhaseTwo) {
		cout << "Reading from hard constraints file...\n";
	}
	else { cout << "Using modified canonical pattern then...\n"; }

	auto start = chrono::steady_clock::now();
	/*****************************************************************************************************************************/
	/*******************************************| Phase 1: Edge-coloring/Latin Square |*******************************************/
	/**********************************************/ cout << "\n\n--- Phase 1:\n";/***********************************************/
	vector<long> M1(n * m);
	M1.resize(n * m);

	if (doPhaseOne) {
	// Edge-coloring/latin square
		latinSquare.resize(n * n);
		bannedRounds.resize(n * n);

	// fill the diagonal with 1's and makes the algorithm ignore the diagonal.
		for (int i = 0; i < n; i++) {
			latinSquare[i * n + i] = 1;
			fixedMatchups.push_back(i * n + i);
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
			cout << "\n\nFound constraints: ";
			while (getline(myFile, line)) {			// Read string 
				stringstream ss(line);				// Make stringstream from s 
				v.clear();							// Clear vector v
				while (ss >> f) {					// Read ints from ss into f 
					v.push_back(f);					// Add f to vector 
				}
				constraintCode = v.operator[](0);	// Read first int of v, which determines the type of constraint
				switch (constraintCode) {
				case 3:
					round = v.operator[](3) % m;	// If the specified round exceeds number of rounds, uses mod.
					bannedRounds[(v.operator[](1) - 1) * n + (v.operator[](2) - 1)] = round + 1;	// Adds the row/col and given round to a banned list
					bannedRounds[(v.operator[](2) - 1) * n + (v.operator[](1) - 1)] = round + 1;	
					cout << "\nHard pairing constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " must NOT play in round " << v.operator[](3);
					break;
				case 5:
					fixedMatchups.push_back((v.operator[](1) - 1) * n + (v.operator[](2) - 1));			// Adds the row/col of latinSquare to an "ignore list"
					fixedMatchups.push_back((v.operator[](2) - 1) * n + (v.operator[](1) - 1));			// Since the matrix is symmetric.
					latinSquare[(v.operator[](1) - 1) * n + (v.operator[](2) - 1)] = v.operator[](3);	// Sets the row/col of latinSquare to the specified round
					latinSquare[(v.operator[](2) - 1) * n + (v.operator[](1) - 1)] = v.operator[](3);
					cout << "\nHard pairing constraint : Team " << v.operator[](1) << " and team " << v.operator[](2) << " must play in round " << v.operator[](3);
					break;

				default:
					break;
				}
			}
			myFile.close();
		}
		else { cout << "Error while trying to open file!"; }

		

		// Print the latinSquare before solving it
		cout << "\n\nLatin square before solving:\n";
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				cout << setw(3) << latinSquare[i * n + j];
			}
			cout << "\n";
		}
		cout << "\n";

		if (solveLatinSquare(latinSquare) == true) {
			// Print latinSquare after a successful solve...
			cout << "\n";
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					cout << setw(3) << latinSquare[i * n + j];
				}
			cout << "\n";
			}

			// Convert latin square to tournament
			int row;
			int col;
			for (int i = 0; i < latinSquare.size(); i++)
			{
				latinSquare[i]--;
				if (latinSquare[i] != 0) {
					row = floor(i / n);
					col = latinSquare[i]-1;
					M1[row * m + col] = i % n + 1;
				}
			}
			cout << "\nM1:";
			printMat(M1, n, m);
		}
		else {
			cout << "No solution exists!" << "\n" << "\n";
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

	cout << "\nM1_2 post phase 1:";
	printMat(M1_2, n, 2 * m);

	// For some reason, this solves a bug concerning the construction of M3 in phase 3
	
	vector<long> M1_safety(n * 2 * m);
	M1_safety = M1_2;
		

	auto postPhaseOne = chrono::steady_clock::now();
	/*****************************************************************************************************************************/
	/*******************************************************| Phase 2: CP |*******************************************************/
	/***********************************************/ cout << "\n\n--- Phase 2:\n";/**********************************************/
	vector<long> M2 (n * m);
	M2.resize(n* m);

	
	if (doPhaseTwo) {
		cout << "Problem-specific constraints detected! \nUsing Cplex to solve the CP model...";

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
			else { cout << "Error while trying to open file!"; }
			cout << "\n---\n";
			// Extract model and set parameters for the solve
			cplex.extract(model);
			// cplex.setOut(env.getNullStream()); // <-- Gets rid of cplex output
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

			if (cplex.solve()) {
				cout << "\n\nSolved cplex model Zuccezzfully!";
				// double OBJval = (double)cplex.getObjValue(); // <-- Get obj val
				// cout << "\nNumber of breaks: " << OBJval << "\n";
			}
		}
		// end try
		catch (IloException& e) { cerr << "\nConcert exception caught: " << e << "\n"; }
		catch (...) { cerr << "\nUnknown exception caught" << "\n"; }
		env.end();
	}
	else {
		cout << "No problem-specific location constraints. \nUsing modified canonical pattern by de Werra (1981)\n";
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
	cout << "\nM2:";
	printMat(M2, n, m);

	// Extends the plan to a full DRR tournament
	vector<long> M2_2 (n * 2 * m);
	M2_2.resize(n * 2 * m);

	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < m; i++) {
			temp[i] = M2[i + m * k];
		}
		for (int j = 0; j < m; j++) {
			M2_2[j + m * (2 * k)] = temp[j];
			M2_2[j + m * (2 * k + 1)] = -temp[j];
		}
	}

	cout << "\nM2_2 post phase 2:";
	printMat(M2_2, n, 2 * m);
	auto postPhaseTwo = chrono::steady_clock::now();
	/*****************************************************************************************************************************/
	/*************************************************| Phase 3: Metaheuristics |*************************************************/
	/***********************************************/ cout << "\n\n--- Phase 3:\n";/**********************************************/
	vector<long> M3 (n * 2 * m);
	M3.resize(n * 2 * m);

	// Create M3 by multiplying entries from the solutions found in the previous phases
	for (int i = 0; i < n * 2 * m; i++) {
		M3[i] = M1_2[i] * M2_2[i];
		//cout << M1_2[i] << "  " << M2_2[i] << "  " << M3[i] << "\n";
	}
	printMat(M3, n, 2 * m);

	//printMat(M1_safety, n, 2 * m);










	// Print elapsed time
	auto end = chrono::steady_clock::now();
	cout << "\n\nElapsed time: "
		<< "\n" << "Phase 1: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseOne - start).count()-700 << " ms"
		<< "\n" << "Phase 2: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(postPhaseTwo - postPhaseOne).count() << " ms"
		<< "\n" << "Phase 3: " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - postPhaseTwo).count() << " ms\n---"
		<< "\n" << "Total  : " << setw(12) << chrono::duration_cast<chrono::milliseconds>(end - start).count()-700 << " ms"
		<< "\n" << "------------------------";

	return 0;
}

/*********************************************************| Refactoring |*********************************************************/

/********* General *********/
// Method for printing out an array as a matrix.
// Args: (Array, # of rows, # of columns)
void printMat(vector<long> arr, int r, int c) {
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

/********* Phase 1 *********/
// Returns a boolean which indicates whether any assigned entry
// in the specified row matches the given number. 
bool usedInRow(vector<int> latSq, int row, int num) {
	for (int col = 0; col < n; col++)
		if (latSq[row * n + col] == num)
		{
			return true;
		}
	return false;
}

// Returns a boolean which indicates whether any assigned entry
// in the specified column matches the given number. 
bool usedInCol(vector<int> latSq, int col, int num) {
	for (int row = 0; row < n; row++)
		if (latSq[row * n + col] == num)
		{
			return true;
		}
	return false;
}

// Returns a boolean which indicates whether the given number breaks a constraint
// either by trying to change a fixed number or assigning a number that is not allowed
bool breaksConstraint(vector<int> latSq, int col, int row, int num) {
	for (int i = 0; i < fixedMatchups.size(); i++)
	{
		if (row * n + col == fixedMatchups[i]) {
			return true;
		}
		if (num == bannedRounds[row * n + col]) {
			return true;
		}
	}
	return false;
}

// Returns a boolean which indicates whether it will be legal to assign
// num to the given row,col location.
bool isSafe(vector<int> latSq, int row, int col, int num) {
	// Check if 'num' is not already placed in current row,
	// and current column
	return !usedInRow(latSq, row, num) &&
		!usedInCol(latSq, col, num) &&
		!breaksConstraint(latSq, row, col, num);
}

// Searches the grid to find an entry that is still unassigned. If
// found, the reference parameters row, col will be set the location
// that is unassigned, and true is returned. If no unassigned entries
// remain, false is returned. 
pair<int, int> getUnassignedLocation(vector<int> latSq) {
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			if (latSq[row * n + col] == 0)
			{
				return make_pair(row, col);
			}
	return make_pair(n,n);
}

// Takes a partially filled-in grid and attempts to assign values to
// all unassigned locations in such a way to meet the requirements
// for a latin square (no duplicates in rows and columns) as well as
// the problem-specific hard constraints
bool solveLatinSquare(vector<int> latSq) {
	// Stops if the grid has been filled
	if (getUnassignedLocation(latSq) == make_pair(n,n))
	{
		cout << " Latin square found!\n";
		Sleep(700);
		latinSquare = latSq;
		return true;
	}

	// Get an unassigned grid location
	pair<int, int> rowAndCol = getUnassignedLocation(latSq);
	int row = rowAndCol.first;
	int col = rowAndCol.second;

	// Consider digits 1 to n (colors/rounds)
	for (int num = 1; num <= n; num++)
	{
		// If placing the current number in the current
		// unassigned location is valid, go ahead
		if (isSafe(latSq, row, col, num))
		{
			// Make tentative assignment
			latSq[row* n + col] = num;
			latSq[col* n + row] = num;


			std::cout << "\rProgress: " << setw(3) << row * n + col << " out of " << setw(3) << latSq.size() << ". Clicking in console will cancel!...";
			std::cout.flush();

			// Do the same thing again recursively. If we go 
			// through all of the recursions, and in the end 
			// return true, then all of our number placements 
			// are valid and we are done
			if (solveLatinSquare(latSq))
			{
				return true;
			}

			// As we were not able to validly go through all 
			// of the recursions, we must have an invalid number
			// placement somewhere. Lets go back and try a 
			// different number for this particular unassigned location
			latSq[row * n + col] = 0;
			latSq[col * n + row] = 0;
		}
	}

	// If we have gone through all possible numbers for the current unassigned
	// location, then we probably assigned a bad number early. So backtrack 
	// and try a different number for the previous unassigned locations.
	return false;
}

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
void swapRounds(vector<long> mat, int k, int l) {
	for (int i = 0; i < 2 * n; i++) {
		long temp = mat[i * m + k - 1];
		mat[i * m + k - 1] = mat[i * m + l - 1];
		mat[i * m + l - 1] = temp;
	}
}

// Method for swapping two rows in a full tournament plan
// Args: (Matrix, 1st row, 2nd row)
void swapRows(vector<long> mat, int k, int l) {
	for (int i = 0; i < 2 * m; i++) {
		long temp = mat[i + k * 2 * m];
		mat[i + k * 2 * m] = mat[i + l * 2 * m];
		mat[i + l * 2 * m] = temp;
	}
}

// Method for swapping all entries in a matrix containing k with l.
// Args: (Matrix, 1st int, 2nd int)
void swapNumbers(vector<long> mat, int k, int l) {
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
void swapTeams(vector<long> mat, int k, int l) {
	swapRows(mat, k - 1, l - 1);
	swapNumbers(mat, k, l);
}