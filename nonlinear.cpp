#include "include/nonlinear.h"
json j;
stringstream ss;
string funcF = "";
string funcG = "";
string dervF = "";
string dervsecF = "";

void enablePres(){
    ss.precision(dbl::max_digits10+1);
}

void prepareCall(){
    enablePres();
    ss.str(std::string());
    j.clear();
}

void setF(string expr){
    funcF = expr;
}

void setG(string expr){
    funcG = expr;
}

void setDervF(string expr){
    dervF = expr;
}

void setDervSecF(string expr){
    dervsecF = expr;
}

double f(double x) {
    FunctionParser fp;
    fp.Parse(funcF, "x");
    double variables[1] = {x};
    return fp.Eval(variables);
}

double fderv(double x){
    FunctionParser fp;
    fp.Parse(dervF, "x");
    double variables[1] = {x};
    return fp.Eval(variables);   
}

double g(double x){
    FunctionParser fp;
    fp.Parse(funcG, "x");
    double variables[1] = {x};
    return fp.Eval(variables);
}

double fsecderv(double x){
    FunctionParser fp;
    fp.Parse(dervsecF, "x");
    double variables[1] = {x};
    return fp.Eval(variables);   
}

string incrementalSearch(double x0, double delta, int niter)
{
    double fx0 = f(x0);
    if (fx0 == 0) {
        ss << x0 << " is a root.";
        return ss.str();
    } else {
        double x1 = x0 + delta;
        int contador = 1;
        double fx1 = f(x1);
        while ((fx0 * fx1 > 0) && (contador < niter)) {
            x0 = x1;
            fx0 = fx1;
            x1 = x0 + delta;
            fx1 = f(x1);
            contador = contador + 1;
        }
        if (fx1 == 0) {
            ss << x1 << " is a root.";
            return ss.str();
        } else if (fx0 * fx1 < 0){
            ss << "There is a root between "<< x0 << " and " << x1 <<".";
            return ss.str();
        } else { 
            ss << "Failure in "<< niter << " iterations." << x0 << " " << x1;
            return ss.str();
        }
    }
}

string bisection (double xi, double xs, double tol, double niter, bool relativeErr){
    double fxi = f(xi);
    double fxs = f(xs);
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    if (fxi == 0){
        ss << xi << " is a root.";
    } else if (fxs == 0) {
        ss << xs << " is a root.";
    } else if (fxi * fxs < 0) {
        double xm = (xi + xs)/2;
        double fxm = f(xm);
        int cont = 1;
        double error = tol + 1;
        while(error > tol && fxm != 0 && cont < niter){
            j[to_string(cont)] = {cont, xi, xs, xm, fxm, error};
            if (fxi*fxm < 0){
                xs = xm;
                fxs = fxm;
            } else {
                xi = xm;
                fxi = fxm;
            }
            double xaux = xm;
            xm = (xi + xs)/2;
            fxm = f(xm);
            error = xm - xaux;
            error = abs(error);
            if (relativeErr){
                error = abs(error/xm);
            }
            cont = cont + 1;
        }
        j[to_string(cont)] = {cont, xi, xs, xm, fxm, error};
        if (fxm == 0){
            ss << xm << " is root.";
        } else if (error < tol){
            ss << xm << " is a root approximation with a tolerance: ";
            ss.precision(2);
            ss << scientific << tol<< " and error: " << error;
            ss << fixed << " found in " << cont << " iterations.";
        } else {
            ss << "Failure in " << niter << "iterations.";
        }
    } else {
        ss << "The provided interval is unsuitable.";
    }
    j["result"] = ss.str();
    cout << j.dump();
    return j.dump();
}

double regulaFalsiXmAux(double a, double b){
    return a - ((f(a)*(b-a))/(f(b)-f(a)));
}

string regulaFalsi (double xi, double xs, double tol, double niter, bool relativeErr){
    double fxi = f(xi);
    double fxs = f(xs);
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    if (fxi == 0){
        ss << xi << " is a root.";
    } else if (fxs == 0) {
        ss << xs << " is a root.";
    } else if (fxi * fxs < 0) {
        double xm = regulaFalsiXmAux(xi, xs);
        double fxm = f(xm);
        int cont = 1;
        double error = tol + 1;
        while(error > tol && fxm != 0 && cont < niter){
            j[to_string(cont)] = {cont, xi, xs, xm, fxm, error};
            if (fxi*fxm < 0){
                xs = xm;
                fxs = fxm;
            } else {
                xi = xm;
                fxi = fxm;
            }
            double xaux = xm;
            xm = regulaFalsiXmAux(xi, xs);
            fxm = f(xm);
            error = xm - xaux;
            error = abs(error);
            if (relativeErr){
                error = abs(error/xm);
            }
            cont = cont + 1;
        }
        j[to_string(cont)] = {cont, xi, xs, xm, fxm, error};
        if (fxm == 0){
            ss << xm << " is root.";
        } else if (error < tol){
            ss << xm << " is a root approximation with a tolerance: ";
            ss.precision(2);
            ss << scientific << tol<< " and error: " << error;
            ss << fixed << " found in " << cont << " iterations.";
        } else {
            ss << "Failure in " << niter << "iterations.";
        }
    } else {
        ss << "The provided interval is unsuitable.";
    }
    j["result"] = ss.str();
    return ss.str();
}

string fixedPoint(double tol, double xa, double niter, bool relativeErr){
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    double fx = f(xa);
    int cont = 0;
    double error = tol+1;
    while (fx != 0 && error > tol && cont < niter){
        j[to_string(cont)]= {cont, xa, fx, error};
        double xn = g(xa);
        fx = f(xn);
        error = abs(xn - xa);
        if(relativeErr){
            error = abs(error/xn);
        }
        xa = xn;
        cont = cont + 1;
    }
    j[to_string(cont)] = {cont, xa, fx, error};
    if (fx == 0){
        ss << xa << " is root with error: ";
        ss.precision(2);
        ss << scientific << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else if (error < tol){
        ss << xa << " is a root approximation with a tolerance: ";
        ss.precision(2);
        ss << scientific << tol<< " and error: " << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else {
        ss << "Failure in " << niter << "iterations.";
    }
    j["result"]=ss.str();
    return j.dump();
}

string newtonMethod(double tol, double x0, double niter, bool relativeErr){
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    double x1 = 0;
    double fx = f(x0);
    double dfx = fderv(x0);
    int cont = 0;
    double error = tol+1;
    while (fx != 0 && error > tol && dfx != 0 && cont < niter){
        j[to_string(cont)]= {cont, x0, fx, dfx, error};
        double x1 = x0 -(fx/dfx);
        fx = f(x1);
        dfx = fderv(x1);
        error = abs(x1 - x0);
        if(relativeErr){
            error = abs(error/x1);
        }
        x0 = x1;
        cont = cont + 1;
    }
    j[to_string(cont)]= {cont, x0, fx, dfx, error};
    if (fx == 0){
        ss << x0 << " is root with error: ";
        ss.precision(2);
        ss << scientific << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else if (error < tol){
        ss << x0 << " is a root approximation with a tolerance: ";
        ss.precision(2);
        ss << scientific << tol<< " and error: " << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else if (dfx == 0){
        ss << x0 << " is a possible multiple root.";
    } else {
        ss << "Failure in " << niter << "iterations.";
    }
    j["result"] = ss.str();
    return j.dump();
}

string secantMethod(double tol, double x0,double x1 ,double niter, bool relativeErr){
    double fx0 = f(x0);
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    j["-1"]= {-1, x0, fx0};
    if (fx0 == 0){
        ss << x0 << " is a root. ";
        ss.str();
    }else{
        double fx1=f(x1);
        j["0"]= {0, x1, fx1};
        int cont = 0;
        double error = tol+1;
        double den = fx1-fx0;
        while(error > tol && fx1 !=0 && den!=0 && cont < niter){
            double x2 = x1 - fx1*(x1-x0)/den;
            error = abs(x2 - x1);
            if(relativeErr){
                error = abs(error/x2);
            }
            x0=x1;
            fx0=fx1;
            x1=x2;
            fx1 = f(x1);
            j[to_string(cont+1)] = {cont+1, x1, fx1, error}; 
            den=fx1-fx0;
            cont=cont+1;
        }
        if(fx1==0){
            ss << x1 << " is a root. ";
        }else if(error<tol){
            ss << x1 << " is root with error: ";
            ss << scientific << error;
            ss << fixed << " found in " << cont << " iterations.";
        }else if(den==0){
            ss << "There's a possible multiple root";
        }else{
            ss << "Failure in " << niter << " iterations.";
        }
    }
    j["result"] = ss.str();
    return j.dump();
}

string multipleRootsMethod(double tol, double x0, double niter, bool relativeErr){
    if(relativeErr){
        j["errorType"] = "relative";
    } else {
        j["errorType"] = "absolute";
    }
    double x1 = 0;
    double fx = f(x0);
    double dfx = fderv(x0);
    double dfx2 = fsecderv(x0);
    int cont = 0;
    double error = tol+1;
    while (fx != 0 && error > tol && dfx != 0 && dfx2 != 0 && cont < niter){
        j[to_string(cont)]= {cont, x0, fx, dfx,dfx2, error};
        double x1 = x0 -(fx*dfx/(pow(dfx,2)-(fx*dfx2)));
        fx = f(x1);
        dfx = fderv(x1);
        dfx2 = fsecderv(x1);
        error = abs(x1 - x0);
        if(relativeErr){
            error = abs(error/x1);
        }
        x0 = x1;
        cont = cont + 1;
    }
    j[to_string(cont)]= {cont, x0, fx, dfx,dfx2, error};
    if (fx == 0){
        ss << x0 << " is root with error: ";
        ss.precision(2);
        ss << scientific << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else if (error < tol){
        ss << x0 << " is a root approximation with a tolerance: ";
        ss.precision(2);
        ss << scientific << tol<< " and error: " << error;
        ss << fixed << " found in " << cont << " iterations.";
    } else if (dfx == 0){
        ss << x0 << " is a possible multiple root.";
    } else {
        ss << "Failure in " << niter << "iterations.";
    }
    j["result"] = ss.str();
    return j.dump();
}

void printRow4Col(){
    for (int i = 0; i<85;i++){
        cout << "-";
    }
    cout <<endl;    
}

void printRow5Col(){
    for (int i = 0; i<107;i++){
        cout << "-";
    }
    cout <<endl;    
}

void printRow6Col(){
    for (int i = 0; i<139;i++){
        cout << "-";
    }
    cout <<endl;    
}

void closedMethodsTable(string jsonString, string methodName){
    int w = 22;
    cout << methodName << " Results:" << endl;
    printRow6Col();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "xi" 
         << '|' << setw(w) << "xu"
         << '|' << setw(w) << "xm"
         << '|' << setw(w) << "F(xm)";
    string error = j["errorType"];
    if (error.compare("absolute")){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< "|" <<endl;
    printRow6Col();
    for (int i = 1; i<j.size()-1;i++){
        vector<double> elements = j[to_string(i)];
        cout.precision(0);
        cout << fixed << '|' << setw(w) << elements[0];
        cout.precision(dbl::max_digits10+1);
        cout << '|' << setw(w) << elements[1]
             << '|' << setw(w) << elements[2]
             << '|' << setw(w) << elements[3];
        cout.precision(1);
        cout <<  scientific << '|' << setw(w) << elements[4];
        if (i == 1){
            cout << '|' << setw(w) << "";
        } else {
            cout << '|' << setw(w) << elements[5];
        }
        cout << '|' << endl;
    }
    printRow6Col();
    cout << j["result"] << endl;
}

void FPTable(string jsonString){
    int w = 20;
    cout << "Fixed Point Results:" << endl;
    printRow4Col();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "Xn" 
         << '|' << setw(w) << "F(Xn)";
    string error = j["errorType"];
    if (error.compare("absolute")){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< "|" <<endl;
    printRow4Col();
    for (int i = 0; i+1<j.size()-1;i++){
        vector<double> elements = j[to_string(i)];
        cout.precision(0);
        cout << fixed << '|' << setw(w) << elements[0];
        cout.precision(dbl::max_digits10+1);
        cout << '|' << setw(w) << elements[1];
        cout.precision(1);
        cout <<  scientific << '|' << setw(w) << elements[2];
        if (i == 0){
            cout << '|' << setw(w) << "";
        } else {
            cout << '|' << setw(w) << elements[3];
        }
        cout << '|' << endl;
    }
    printRow4Col();
    cout << j["result"] << endl;
}

void NewtonTable(string jsonString){
    int w = 20;
    cout << "Newton Method Results:" << endl;
    printRow5Col();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "Xn" 
         << '|' << setw(w) << "F(Xn)"
         << '|' << setw(w) << "F'(Xn)";
    string error = j["errorType"];
    if (error.compare("absolute")){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< " |" <<endl;
    printRow5Col();
    for (int i = 0; i+1<j.size()-1;i++){
        vector<double> elements = j[to_string(i)];
        cout.precision(0);
        cout << fixed << '|' << setw(w) << elements[0];
        cout.precision(dbl::max_digits10+1);
        cout << fixed << '|' << setw(w) << elements[1];
        cout.precision(1);
        cout <<  scientific << '|' << setw(w) << elements[2];
        cout.precision(dbl::max_digits10+1);
        cout  << fixed << '|' << setw(w) << elements[3];
        cout.precision(1);
        if (i == 0){
            cout << '|' << setw(w) << "";
        } else {
            cout <<  scientific << '|' << setw(w) << elements[4];
        }
        cout << '|' << endl;
    }
    printRow5Col();
    cout << j["result"] << endl;
}

void secantTable(string jsonString){
    int w = 20;
    cout << " Secant Method Results:" << endl;
    printRow4Col();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "Xn" 
         << '|' << setw(w) << "F(Xn)";
    string error = j["errorType"];
    if (error.compare("absolute")){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< "|" <<endl;
    printRow4Col();
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
    printRow4Col();
    cout << j["result"] << endl;
}

void MRTable(string jsonString){
    int w = 22;
    cout << "Multiple Roots Method Results:" << endl;
    printRow6Col();
    cout << '|' << setw(w) << "n" 
         << '|' << setw(w) << "Xn" 
         << '|' << setw(w) << "F(Xn)"
         << '|' << setw(w) << "F'(Xn)"
         << '|' << setw(w) << "F''(Xn)";
    string error = j["errorType"];
    if (error.compare("absolute")){
        cout << '|' << setw(w) << "Relative Error";
    } else {
        cout << '|' << setw(w) << "Absolute Error";
    }
    cout<< " |" <<endl;
    printRow6Col();
    for (int i = 0; i+1<j.size()-1;i++){
        vector<double> elements = j[to_string(i)];
        cout.precision(0);
        cout << fixed << '|' << setw(w) << elements[0];
        cout.precision(dbl::max_digits10+1);
        cout << fixed << '|' << setw(w) << elements[1];
        cout.precision(1);
        cout <<  scientific << '|' << setw(w) << elements[2];
        cout.precision(dbl::max_digits10+1);
        cout  << fixed << '|' << setw(w) << elements[3];
        cout  << fixed << '|' << setw(w) << elements[4];
        cout.precision(1);
        if (i == 0){
            cout << '|' << setw(w) << "";
        } else {
            cout <<  scientific << '|' << setw(w) << elements[5];
        }
        cout << '|' << endl;
    }
    printRow6Col();
    cout << j["result"] << endl;
}

int main(int argc, char *argv[]) {
    //setF("x^3 + 4*(x^2)-10"); //FP Newton sec
    //setF("(x^4)-(18*x^2)+81"); //multiple roots
    
    //setDervF("x * (8 + (3 * x))"); //Newton
    //setDervF("3*x*sin(3*x)-cos(3*x)+3*exp(3*x-12)-2*x+4");
    //setDervSecF("(12*x^2)-36");
    prepareCall();
    //setF("cos((7*x)-8)*exp(-(x^2)+4)+log((x^4)+3)-x-15");
    //setF("exp(-(x^2)+3)-5*x^2");
    //setG("(exp(-(x^2)+3)*((2*x^2)+1)+5*x^2)/((2*x)*(exp(-(x^2)+3)+5))");
    //setF("(x^4)-(7.45*x^3)+(19.0956*x^2)-(20.0471*x)+7.14292992");
    //setF("exp(-(x^2)+3*sin((3*x+4)))*log((x^2)+4)-3");
    //setDervF("-20.0471+(38.1912*x)-(22.35*x^2)+(4*x^3)");
    //setDervSecF("(39.1912)-(44.7*x)+(12*x^2)");
    //setDervF("(5*x^4)+(16*x^3)+(3*x^2)-(20*x)-4");
    //setDervSecF("(20*x^3)+(48*x^2)+(6*x)-(20)");

    //setG("sqrt(exp((3*x)-12)-(x*cos(3*x))+(4*x))");
    //cout << incrementalSearch(-5, 1, 20);
    setF("exp(-(x^2)+1)-(4*x^3)+25"); //bisection regula-falsi
    closedMethodsTable(bisection(1, 2, 10e-7, 20, false), "Bisection Method");
    //closedMethodsTable(regulaFalsi(1, 2, 0.001, 20, false), "Regula Falsi Method");
    //NewtonTable(newtonMethod(10e-8,3.7855,20,true));
    //MRTable(multipleRootsMethod(10e-5,1.4, 40, true));
    //ss.str(std::string());
    secantTable(secantMethod(5.0*10e-8,0.8, 0.9,20,true));
    return 0;
}
