#include "include/printer.h"

void printRow(){
    for (int i = 0; i<85;i++){
        cout << "-";
    }
    cout <<endl;    
}

void secantTable(string jsonString){
    int w = 20;
    cout << "Secant Method Results:" << endl;
    j.parse(jsonString);
    cout << jsonString;
    printRow();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "Xn" 
         << '|' << setw(w) << "F(Xn)";
    string error = "\"relative\"";
    cout << j["errorType"];
    if (error.compare(j["errorType"])){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< "|" <<endl;
    printRow();
    for (int i = -1; i+1<j.size()-2;i++){
        vector<double> elements = j[to_string(i)];
        cout.precision(0);
        cout << fixed << '|' << setw(w) << elements[0];
        cout.precision(dbl::max_digits10+1);
        cout << '|' << setw(w) << elements[1];
        cout.precision(1);
        cout <<  scientific << '|' << setw(w) << elements[2];
        if (elements.size() == 3){
            cout << '|' << setw(w) << "";
        } else {
            cout << '|' << setw(w) << elements[3];
        }
        cout << '|' << endl;
    }
    printRow();
    cout << j["result"] << endl;
}
