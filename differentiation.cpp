#include "include/linearecs.h"

void enablePres(){
    cout.precision(dbl::max_digits10+1);
}

void dosPuntos(double fx0, double fx1, double h){
    double res = (fx1 - fx0)/h;
    cout << "f'(x) = " << res << " - " << h/2 << "f''(E)" << endl;
}

void tresPuntosProgresiva(double fx0, double fx1, double fx2,double h){
    double res = ((-3*fx0)+(4*fx1)-(fx2))/(2*h);
    cout << "f'(x) = " << res << " - " << pow(h,2)/3 << "f'''(E)" << endl;
}

void tresPuntosRegresiva(double fx0, double fx1, double fx2,double h){
    double res = ((fx2)-(4*fx1)+(3*fx0))/(2*h);
    cout << "f'(x) = " << res << " - " << pow(h,2)/3 << "f'''(E)" << endl;
}

//asumiendo que x1 es el valor del centro
void tresPuntosCentrada(double fx0, double fx2,double h){
    double res = (-(fx0)+(fx2))/(2*h);
    cout << "f'(x) = " << res << " - " << pow(h,2)/6 << "f'''(E)" << endl;
}
//De x0, x1, x2, x3, x4, asumiendo que x2 es el valor del medio:
void cincoPuntosCentrada(double fx0, double fx1,double fx3, double fx4, double h){
    double res = (fx0-(8*fx1)+(8*fx3)-(fx4))/(12*h);
    cout << "f'(x) = " << res << " - " << pow(h,4)/30 << "f''''(E)" << endl;
}
//funciona para adelante y para atras, depende del orden de los valores
void cincoPuntos(double fx0, double fx1,double fx2, double fx3, double fx4, double h){
    double res = ((-25*fx0)+(48*fx1)-(36*fx2)+(16*fx3)-(3*fx4))/(12*h);
    cout << "f'(x) = " << res << " - " << pow(h,4)/5 << "f'''''(E)" << endl;
}

int main(int argc, char *argv[]) {
    enablePres();
    dosPuntos(-11.8646, -12.4775, 0.1);
    tresPuntosProgresiva(-13.0891, -13.6997, -14.3092,0.1);
    tresPuntosRegresiva(-13.0891, -12.4775, -11.8646,0.1);
    tresPuntosCentrada(-12.4775,-13.6997,0.1);
    cincoPuntosCentrada(-11.8646, -12.4775, -13.6997, -14.3092, 0.1);
    cincoPuntos(-11.8646, -12.4775, -13.0891, -13.6997, -14.3092, 0.1);
    return 0;
}