#include "include/linearecs.h"
#include <iostream>

json j;
vector<int> marks;
int steps = 0;

void addStep(vector<vector<double>> M){
    j[to_string(steps)] = M;
    steps++;
}

void initMarks(int n){
    marks.clear();
    for (int i = 0; i < n; i++) {
        marks.push_back(i+1);
    }
}

vector<vector<double>> createMatrix(string json){
    j.parse(json);
    vector<vector<double>> matrix;
    for (int i = 0; i<j.size();i++){
        matrix.push_back(j[to_string(i)]);
    }
    j.clear();
    addStep(matrix);
    return matrix;
}

vector<vector<double>> enterMatrix(){
    int n;
    cout << "Enter the N size of the matrix:" << endl;
    cin >> n;
    for (int i = 0; i<n;i++){
        vector<double> row;
        for (int c = 0;c<n;c++){
            cout << "Enter the element i = " << i << " j = " << c << endl;
            double data;
            cin >> data;
            row.push_back(data);
        }
        j[to_string(i)] = row;
    }
    return createMatrix(j.dump());
}

void printMatrix (vector<vector<double>> A){
    int w = 20;
    for (int i = 0;i<A.size();i++){
        for(int c = 0; c<A[i].size();c++){
            if (abs(A[i][c])>9.99e-16){
                cout << '|' << setw(w) << A[i][c] << " ";
            } else {
                cout << '|' << setw(w) << 0 << " ";
            }
        }
        cout << endl;
    }
}

void printNewtonTable (vector<vector<double>> A){
    int w = 20;
    for (int i = 0;i<A.size();i++){
        for(int c = 0; c<A[i].size();c++){
            if (abs(A[i][c])>9.99e-16){
                    cout << '|' << setw(w) << A[i][c] << " ";
            } else {
                    cout << '|' << setw(w) << 0 << " ";
            }
        }
        cout << endl;
    }
}

vector<vector<double>> createAugmentedMatrix(vector<vector<double>> A, vector<double> b){
    for (int i = 0; i<A.size();i++){
        A[i].push_back(b[i]);
    }
    addStep(A);
    return A;
}

vector<vector<double>> mulElimination(vector<vector<double>> Ab, int k, int n){
    for (int i = k+1; i<n;i++){
        double multiplier = Ab[i][k]/Ab[k][k];
        for(int c = k;c<n+1;c++){
            Ab[i][c] = Ab[i][c] - multiplier*Ab[k][c]; 
        }
    }
    addStep(Ab);
    return Ab;
}

vector<vector<double>> elimination(vector<vector<double>> A, vector<double> b){
    int n = A.size();
    vector<vector<double>> Ab = createAugmentedMatrix(A, b);
    for (int k = 0; k<A.size()-1;k++){
        Ab = mulElimination(Ab, k, n);
    }
    addStep(Ab);
    return Ab;
}

vector<vector<double>> regresiveSubstitution(vector<vector<double>> Ab){
    int n = Ab.size()-1;
    vector<double> res(n);
    double xi = Ab[n][n+1]/Ab[n][n];
    res.push_back(xi);
    for (int i = n;i>-1;i--){
        double sum = 0;
        for (int p = i+1; p<n+1;p++){
            sum = sum + Ab[i][p]*res[p];
        }
        xi = (Ab[i][n+1]-sum)/Ab[i][i];
        res[i]=xi;
    }
    vector<vector<double>> mRes;
    mRes.push_back(res);
    addStep(mRes);
    return mRes;
}

vector<vector<double>> progressiveSubstitution(vector<vector<double>> Ab){
    int n = Ab.size();
    vector<double> res(n);
    double xi = Ab[0][n]/Ab[0][0];
    for (int i = 0;i<n;i++){
        double sum = 0;
        for (int p = 0; p<n;p++){
            sum = sum + Ab[i][p]*res[p];
        }
        xi = (Ab[i][n]-sum)/Ab[i][i];
        res[i]=xi;
    }
    vector<vector<double>> mRes;
    mRes.push_back(res);
    addStep(mRes);
    return mRes;
}

void printXVector(vector<vector<double>> x){
    for (int i = 0; i<x[0].size();i++){
        cout << "x" << i+1 << " = " << x[0][i] << endl;
    }
}

void simpleGaussianElimination(vector<vector<double>> A, vector<double> b){
    vector<vector<double>> U = elimination(A, b);
    vector<vector<double>> x = regresiveSubstitution(U);
}

vector<vector<double>> switchCols(vector<vector<double>> Ab, int maxCol, int k){
    for (int i = 0; i<Ab.size();i++){
        double aux = Ab[i][k];
        Ab[i][k] = Ab[i][maxCol];
        Ab[i][maxCol] = aux;
    }
    return Ab;
}

vector<vector<double>> switchRows(vector<vector<double>> Ab, int maxRow, int k){
    vector<double> rowK = Ab[k];
    vector<double> vmaxRow = Ab[maxRow];
    Ab[k] = vmaxRow;
    Ab[maxRow] = rowK;
    return Ab;
}

void switchMarks(int maxCol, int k){
    int aux = marks[k];
    marks[k] = marks[maxCol];
    marks[maxCol] = aux;
}

vector<vector<double>> partialPivoting(vector<vector<double>> Ab, double n, int k){
    double max = abs(Ab[k][k]);
    int maxRow = k;
    for (int s = k+1;s<n;s++){
        if(abs(Ab[s][k])>max) {
            max = abs(Ab[s][k]);
            maxRow = s;
        }
    }
    if (max == 0){
        j["print"] = "The system doesn't have an unique solution.";
    } else {
        if (maxRow != k){
            Ab = switchRows(Ab, maxRow, k);
        }
        addStep(Ab);
        return Ab;
    }
}

vector<vector<double>> totalPivoting(vector<vector<double>> Ab, double n, int k){
    double max = 0;
    int maxRow = k;
    int maxCol = k;
    for (int r = k; r<n;r++){
        for (int s = k; s<n;s++){
            if(abs(Ab[r][s]) > max){
                max = abs(Ab[r][s]);
                maxRow = r;
                maxCol = s;
            }
        }
    }
    if(max == 0){
        j["print"] = "The system doesn't have an unique solution.";
    } else {
        if(maxRow != k){
            Ab = switchRows(Ab, maxRow, k);
        }
        if (maxCol != k){
            Ab = switchCols(Ab, maxCol, k);
            switchMarks(maxCol, k);
        }
        addStep(Ab);
        return Ab;
    }
}

void printXVectorMarks(vector<vector<double>> x){
    for (int i = 0; i<x[0].size();i++){
        cout << "x" << marks[i] << " = " << x[0][i] << endl;
    }
}

vector<vector<double>> eliminationPP(vector<vector<double>> A, vector<double> b){
    int n = A.size();
    vector<vector<double>> Ab = createAugmentedMatrix(A, b);
    for (int k = 0; k<A.size()-1;k++){
        Ab = partialPivoting(Ab, Ab.size(), k);
        Ab = mulElimination(Ab, k, n);
    }
    return Ab;
}

vector<vector<double>> eliminationTP(vector<vector<double>> A, vector<double> b){
    vector<vector<double>> Ab = createAugmentedMatrix(A, b);
    int n = A.size();
    initMarks(Ab.size());
    for (int k = 0; k<A.size()-1;k++){
        Ab = totalPivoting(Ab, n, k);
        Ab = mulElimination(Ab, k, n);
    }
    return Ab;
}

void PPGaussianElimination(vector<vector<double>> A, vector<double> b){
    vector<vector<double>> U = eliminationPP(A, b);
    vector<vector<double>> x = regresiveSubstitution(U);
    printXVector(x);
}

void PTGaussianElimination(vector<vector<double>> A, vector<double> b){
    vector<vector<double>> U = eliminationTP(A, b);
    vector<vector<double>> x = regresiveSubstitution(U);
    printXVectorMarks(x);
}

void enablePres(){
    cout.precision(dbl::max_digits10+1);
}

vector<vector<double>> initLU(int n){
    vector<vector<double>> res;
    for (int i = 0;i<n;++i){
        vector<double> newRow;
        res.push_back(newRow);
        for(int j= 0;j<n;++j){
            if (i == j){
                res[i].push_back(1);
            } else {
                res[i].push_back(0);
            }
        }
    }
    addStep();
    return res;
}

vector<vector<vector<double>>> DoolittleFact(vector<vector<double>> A){
    int n = A.size();
    vector<vector<double>> L = initLU(n);
    vector<vector<double>> U = initLU(n);
    for (int k = 0; k<n;k++){
        double sum1 = 0;
        for (int p=0; p<k;p++){
            sum1 = sum1 + L[k][p]*U[p][k];
        }
        U[k][k]=A[k][k] - sum1;
        addStep(U);
        for (int i = k + 1; i < n; i++) {
            double suma2 = 0;
            for (int p = 0; p < k; p++) {
                suma2 += L[i][p] * U[p][k];
            }
        L[i][k] = (A[i][k] - suma2) / U[k][k];
        addStep(L);
        }

        for (int j = k + 1; j < n; j++) {
            double suma3 = 0;
            for (int p = 0; p < k; p++){
                suma3 += L[k][p] * U[p][j];
            } 
            U[k][j] = A[k][j] - suma3;
            addStep(U);
        }
    }
    vector<vector<vector<double>>> LU;
    LU.push_back(L);
    LU.push_back(U);
    return LU;
}

vector<vector<vector<double>>> CroultFact(vector<vector<double>> A){
    int n = A.size();
    vector<vector<double>> L = initLU(n);
    vector<vector<double>> U = initLU(n);
    for (int k = 0; k<n;k++){
        double sum1 = 0;
        for (int p=0; p<k;p++){
            sum1 = sum1 + L[k][p]*U[p][k];
        }
        L[k][k]=A[k][k] - sum1;
        addStep(L);
        for (int i = k + 1; i < n; i++) {
            double suma2 = 0;
            for (int p = 0; p < k; p++) {
                suma2 += L[i][p] * U[p][k];
            }
            L[i][k] = A[i][k] - suma2;
            addStep(L);
        }

        for (int j = k + 1; j < n; j++) {
            double suma3 = 0;
            for (int p = 0; p < k; p++){
                suma3 += L[k][p] * U[p][j];
            } 
            U[k][j] = (A[k][j] - suma3) / L[k][k];
            addStep(U);
        }
    }
    vector<vector<vector<double>>> LU;
    LU.push_back(L);
    LU.push_back(U);
    return LU;
}

vector<vector<vector<double>>> CholeskyFact(vector<vector<double>> A){
    int n = A.size();
    vector<vector<double>> L = initLU(n);
    vector<vector<double>> U = initLU(n);
    for (int k = 0; k<n;k++){
        double sum1 = 0;
        for (int p=0; p<k;p++){
            sum1 = sum1 + L[k][p]*U[p][k];
        }
        L[k][k]=sqrt(A[k][k] - sum1);
        U[k][k]=L[k][k];
        addStep(L);
        addStep(U);
        for (int i = k + 1; i < n; i++) {
            double sum2 = 0;
            for (int p = 0; p < k; p++) {
                sum2 = sum2 + (L[i][p] * U[p][k]);
            }
            L[i][k] = (A[i][k] - sum2)/U[k][k];
            addStep(L);
        }

        for (int j = k + 1; j < n; j++) {
            double sum3 = 0;
            for (int p = 0; p < k; p++){
                sum3 = sum3 + (L[k][p] * U[p][j]);
            } 
            U[k][j] = (A[k][j] - sum3) / L[k][k];
            addStep(U);
        }
    }
    vector<vector<vector<double>>> LU;
    LU.push_back(L);
    LU.push_back(U);
    return LU;
}

void Croult (vector<vector<double>> A, vector<double> b){
    vector<vector<vector<double>>> LU = CroultFact(A);
    vector<vector<double>> L = LU[0];
    vector<vector<double>> U = LU[1];
    vector<vector<double>> Lb = createAugmentedMatrix(L, b);
    vector<vector<double>> z = progressiveSubstitution(Lb);
    vector<vector<double>> Uz = createAugmentedMatrix(U, z[0]);
    vector<vector<double>> x = regresiveSubstitution(Uz);
    printXVector(x);
}

void Doolittle (vector<vector<double>> A, vector<double> b){
    vector<vector<vector<double>>> LU = DoolittleFact(A);
    vector<vector<double>> L = LU[0];
    vector<vector<double>> U = LU[1];
    vector<vector<double>> Lb = createAugmentedMatrix(L, b);
    vector<vector<double>> z = progressiveSubstitution(Lb);
    vector<vector<double>> Uz = createAugmentedMatrix(U, z[0]);
    vector<vector<double>> x = regresiveSubstitution(Uz);
    printXVector(x);
}

void Cholesky (vector<vector<double>> A, vector<double> b){
    vector<vector<vector<double>>> LU = CholeskyFact(A);
    vector<vector<double>> L = LU[0];
    printMatrix(L);
    vector<vector<double>> U = LU[1];
    printMatrix(U);
    vector<vector<double>> Lb = createAugmentedMatrix(L, b);
    vector<vector<double>> z = progressiveSubstitution(Lb);
    vector<vector<double>> Uz = createAugmentedMatrix(U, z[0]);
    vector<vector<double>> x = regresiveSubstitution(Uz);
    printXVector(x);
}

double normaUniforme(vector<double> x) {
  double max = 0;
  for (int i = 0; i<x.size();i++){
      if(max < abs(x[i])){
          max = abs(x[i]);
      }
  }
  return max;
};

double norma1(vector<double> x) {
  double sum = 0; 
  for (int i = 0; i<x.size();i++){
      sum = sum + abs(x[i]);
  }
  return sum;
};

double norma2(vector<double> x){
  double sum = 0;
  for (int i = 0; i<x.size();i++){
      sum = sum + pow(x[i], 2);
  }
  sum = sqrt(sum);
  return sum;
};

vector<double> calcNewJacobi(vector<double> vector0, vector<vector<double>> A,
vector<double> b, double w) {
  int n = vector0.size();
  vector<double> vector1;
  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < n; j++) {
      if (j != i){
          sum = sum + (A[i][j] * vector0[j]);
      } 
    }
    double xi = (b[i] - sum)/A[i][i];
    vector1.push_back((w*xi)+(1-w)*vector0[i]);
  }
  
  return vector1;
};

vector<double> calcNewSeidel(vector<double> vector0, vector<vector<double>> A,
vector<double> b, double w) {
  int n = vector0.size();
  vector<double> vector1;
  for(int i = 0; i<n;i++){
      vector1.push_back(vector0[i]);
  }
  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < n; j++) {
      if (j != i){
          sum = sum + (A[i][j] * vector1[j]);
      } 
    }
    double xi = (b[i] - sum)/A[i][i];
    vector1[i]=((w*xi)+(1-w)*vector0[i]);
  }
  return vector1;
};

void jacobi (double tolerance, vector<double> vector0,int niter, vector<vector<double>> A,
vector<double> b, int norma, double w){
    double dispertion = tolerance + 1;
    double counter = 0;
    
  while (dispertion > tolerance && counter < niter) {
      vector<double> vector1 = calcNewJacobi(vector0, A, b, w);
      vector<double> vecAux;
      for(int i = 0; i<vector0.size();++i){
          vecAux.push_back(vector1[i]-vector0[i]);
      }
    if (norma == 1){
        dispertion = norma1(vecAux)/norma1(vector1);
    } else if (norma == 2){
        dispertion = norma2(vecAux)/norma2(vector1);
    } else {
        dispertion = normaUniforme(vecAux)/normaUniforme(vector1);
    }
    j[to_string(counter)]= {vector0, dispertion};
    cout << j[to_string(counter)].dump() << endl;
    vector0 = vector1;
    counter++;
  }
  j[to_string(counter)]= {vector0, dispertion};
  cout << j[to_string(counter)].dump() <<endl;
  if (dispertion < tolerance) {
    vector<vector<double>> res;
    res.push_back(vector0);
    printXVector(res);
    cout << " is an aproximation with tolerance " << tolerance << "and dispertion " << dispertion <<
     " found in " << counter << " iterations."<<endl;   
  } else {
    cout << "Failure in " << niter << "iterations";
  }
}

void gaussSeidel (double tolerance, vector<double> vector0,int niter, vector<vector<double>> A,
vector<double> b, int norma, double w){
    double dispertion = tolerance + 1;
    double counter = 0;
    
  while (dispertion > tolerance && counter < niter) {
      vector<double> vector1 = calcNewSeidel(vector0, A, b, w);
      vector<double> vecAux;
      for(int i = 0; i<vector0.size();++i){
          vecAux.push_back(vector1[i]-vector0[i]);
      }
    if (norma == 1){
        dispertion = norma1(vecAux)/norma1(vector1);
    } else if (norma == 2){
        dispertion = norma2(vecAux)/norma2(vector1);
    } else {
        dispertion = normaUniforme(vecAux)/normaUniforme(vector1);
    }
    j[to_string(counter)]= {vector0, dispertion};
    cout << j[to_string(counter)].dump();
    vector0 = vector1;
    counter++;
  }
  j[to_string(counter)]= {vector0, dispertion};
  cout << j[to_string(counter)].dump();
  if (dispertion < tolerance) {
    vector<vector<double>> res;
    res.push_back(vector0);
    printXVector(res);
    cout << " is an aproximation with tolerance " << tolerance << "and dispertion " << dispertion <<
     " found in " << counter << " iterations."<<endl;   
  } else {
    cout << "Failure in " << niter << "iterations";
  }
}

void showSteps(){
    int last = 0;
    cout << "Show Steps? y/n"<<endl;
    string answer = "";
    cin >> answer;
    if(answer.compare("n")){
        for (int i=0; i<j.size();i++){
            cout << "Step " << i+1 << ":" << endl; 
            vector<vector<double>> S = j[to_string(i)];
            printMatrix(S);
            cout << endl;
            last = i;
        }
    }
}
//ESTO ESTA MALO
void evalNewton(vector<vector<double>> polTable, double x){
    int n = polTable.size();
    double res = 0;
    for(int i = 0; i<n-1;i++){
        double mul = 1;
        for (int j = i-1; j>=0;j--){
            mul = mul * (x-polTable[j][0]);
        }
        mul = mul * polTable[i+1][i];
        res = res + mul;
    }
    cout << res << endl;
}

void newtonInterpolation(vector<double> xn, vector<double> fxn){
    int n = xn.size();
    vector<double> prevstep = fxn;
    vector<vector<double>> tabla;
    tabla.push_back(xn);
    tabla.push_back(fxn);
    for(int i = 1; i<n;i++){
        vector<double> step(n);
        for(int j = i; j<n;j++){
            step[j] = (prevstep[j]-prevstep[j-1])/(xn[j]-xn[j-i]);
        }
        tabla.push_back(step);
        prevstep = step;
    }
    printNewtonTable(tabla);
    //evalNewton(tabla, 2.5);
}

double lagrangeInterpolation(vector<double> xn, vector<double> fxn, double xi){
    double result = 0; // Initialize result
    int n = fxn.size(); 
    for (int i=0; i<n; i++) { 
        // Compute individual terms of above formula 
        double term = fxn[i]; 
        for (int j=0;j<n;j++)  { 
            if (j!=i) 
                term = term*(xi - xn[j])/double(xn[i] - xn[j]); 
        } 
        // Add current term to result 
        result += term; 
    } 
    return result; 
}

int main(int argc, char *argv[]) {
    enablePres();
    //vector<vector<double>> A = enterMatrix();
    //simpleGaussianElimination(A, {12,32,-24,14});
    //PPGaussianElimination(A, {-12,13,31,-32});
    //PTGaussianElimination(A, {-12,13,31,-32});
    //Croult(A, {-20, 69, 96, -32}); //Tambien funciona para Doolittle
    //jacobi(5.0e-6,{0,0,0}, 20, A, {-23, 5, 34}, 2);
    //gaussSeidel(5.0e-6,{0,0,0}, 20, A, {-23, 5, 34}, 2, 1); //normal
    //gaussSeidel(5.0e-6,{0,0,0}, 30, A, {24, 30, -24}, 2, 1.3); //Relajado
    //showSteps();
    //newtonInterpolation({2,2.2,2.4,2.6,2.8},{-4.6109,-4.1749,-3.3768,-2.1362,-0.3553});
    cout << lagrangeInterpolation({-1,1,2,4},{7,-1,-8,2}, 2.5);
    return 0;
}