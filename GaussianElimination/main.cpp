#include <iostream>
#include <cmath>

using namespace std;


void swap(double *a, double *b) {
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

void multDiffRowStep(double *matrix, double copyCurRowByValue[], int numOfVariables, int rowDealWith, int curColumn) {
	double factor = matrix[rowDealWith + curColumn] / copyCurRowByValue[curColumn];
	for (int i = curColumn; i < numOfVariables + 1; i++) {
		copyCurRowByValue[i] *= factor;
		matrix[rowDealWith + i] -= copyCurRowByValue[i];
	}
}

int main(void)
{
	int numOfEquations, numOfVariables;

	//	read the number of rows, columns
    cin >> numOfEquations >> numOfVariables;
    if (!(numOfEquations >= 1) || !(numOfVariables >= 1)) {
		return 1;
	}

	//matrix reading
    double *matrix = new double[numOfEquations * (numOfVariables + 1)];
	for (int i = 0; i < numOfEquations * (numOfVariables + 1); i += numOfVariables + 1) {
		for (int j = 0; j < numOfVariables + 1; j++) {
			cin >> matrix[i + j];
		}
	} 

    // direct path of method 
    double pivotEl = 0; 
	int curRow = 0;     
	int curColumn = 0;
	double *copyCurRowByValue = new double[numOfVariables + 1];
	
    while (curRow != numOfEquations - 1) {
		
		//host element selection
        int  rowPivotIndex;
        for (int i = curRow * (numOfVariables + 1); i < numOfEquations * (numOfVariables + 1); i += numOfVariables + 1) {
			if (pivotEl < abs(matrix[i + curColumn])) {
				pivotEl = abs(matrix[i + curColumn]);
				rowPivotIndex = i;
			}
		}

		//exit condition
        if (matrix[rowPivotIndex + curColumn] == 0) {
			curRow++;
			curColumn++;
			continue;
		}

		//move the line with leading elements to the place of the first line
        if (curRow * (numOfVariables + 1) != rowPivotIndex) {
			for (int i = curColumn; i < numOfVariables + 1; i++) {
				swap(&matrix[curRow * (numOfVariables + 1) + i], &matrix[rowPivotIndex + i]);
			}
		}

		//subtraction of the row with the leading element multiplied by the corresponding coefficient from the remaining ones
		for (int i = (curRow + 1) * (numOfVariables + 1); i < numOfEquations * (numOfVariables + 1); i += numOfVariables + 1) {
			if (matrix[i + curColumn] != 0) {
                for (int i = curColumn; i < numOfVariables + 1; i++) {
					copyCurRowByValue[i] = matrix[i + curRow * (numOfVariables + 1)];
				}
                multDiffRowStep(matrix, copyCurRowByValue, numOfVariables, i, curColumn);
			}
		}

		//transition to the next minor by deleting the current row and column
        curRow++;
		curColumn++;
        pivotEl = 0;
	}

	//calculation of the rank of the resulting matrix with sifting out the possibilities of the absence of roots
	int rank = 0;
	bool isNoSolution = false;

	for (int i = 0; i < (numOfVariables  + 1) * numOfEquations; i += numOfVariables + 1) {
		bool isZeroRow = true;
		for (int j = i; j < i + numOfVariables + 1; j++) {
			if (matrix[j] != 0 && j != i + numOfVariables) {
				isZeroRow = false;
				break;
			}
			if (j == i + numOfVariables && isZeroRow == true) {
				isNoSolution = true;
				break;
			}
		}
		if (isNoSolution == true) {
			break;
		}
		if (isZeroRow == false) {
			rank++;
		}
	}

	if (isNoSolution == true) {
		cout << "NO" << endl;
	}

	//checking for an infinite number of solutions using the Kronecker-Capelli theorem
	else if (rank < numOfVariables) {
		cout << "INF";
	}

	//if the solution is unique, then we apply the reverse move
	else if (rank == numOfVariables) {
		double *solutions = new double[numOfVariables];
		int numOfSolutions = numOfVariables - 1;
		while (curRow >= 0) {
			double factor = matrix[curRow * (numOfVariables + 1) + curColumn];
			double sum = 0;
			for (int i = 0; i < numOfVariables + 1; i++) {
				if (curRow * (numOfVariables + 1) + curColumn + 1 + i != (curRow + 1) * (numOfVariables + 1) - 1) {
					sum += -matrix[curRow * (numOfVariables + 1) + curColumn + 1 + i] * solutions[curColumn + 1 + i];
				}
				else {
					sum += matrix[curRow * (numOfVariables + 1) + curColumn + 1 + i];
					break;
				}
			}
			solutions[numOfSolutions] = sum / factor;
			numOfSolutions--;
			curColumn--;
			curRow--;
		}

		//decision vector output
		cout.precision(16);
		cout << "YES" << endl;
		for (int i = 0; i < numOfVariables; i++) { 
			cout << solutions[i] <<" ";
		}
		delete solutions;
	}
	delete matrix;
	delete copyCurRowByValue;
	return 0;
}


