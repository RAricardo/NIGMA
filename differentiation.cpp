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

void trapecio(double a, double b, double fa, double fb){
    double res = (((b-a)/2)*(fa+fb));//+;
    cout << res << " + " << ((pow((b-a),3))/12) << "f''(E)" << endl;
}

void trapecioGeneralizado(vector<double> xn, vector<double> fx, double h){
    int n= xn.size();
    double sum = 0;
    for(int i = 1; i<n-1;i++){
        sum = sum + fx[i];   
    }
    sum = fx[0] + sum*2 + fx[n-1];

    double res = (h/2)*sum;
    cout << res << " + " << ((pow((h),3))/(12*pow(n,3))) << " SUM de i=1  a n de f''(Ei)" << endl;
}

void simpsonunterc(double fx0, double fx1, double fx2, double h){
    double res = (h/3)*(fx0 + (4*fx1) +fx2);
    cout << res << " - " << ((pow((h),5))/(90)) << " f''''(Ei)" << endl;
}

//par
void simpsonGeneralizado(vector<double> xn, vector<double> fx, double h){
    int n= xn.size();
    double sumImpar = 0;

    for(int i = 1; i<n-1;i++){
        if(!i%2==0){
            sumImpar = sumImpar + fx[i];
        }
    }

    sumImpar = sumImpar*4;

    double sumPar = 0;

    for(int i = 2; i<n-2;i++){
        if(i%2==0){
            sumPar = sumPar + fx[i];
        }
    }

    sumPar = sumPar*2;

    double sum = fx[0] + sumImpar + sumPar + fx[n-1];

    double res = (h/3)*sum;
    cout << res << " - (" << (n*(pow((h),5))) << " f''''(u))/180" << endl;
}

//3/8 multiplos de 3
void simpsontresOct(double fx0,double fx1,double fx2,double fx3,double h){
    double res = (3*h/8)*(fx0 + (3*fx1) +(3*fx2) + (fx3));
    cout << res << endl;
}

void simpsontresOctGeneralizado(vector<double> xn, vector<double> fx, double h){
    int n= xn.size();
    double primerSuma = 0;

    for(int i = 1; i<n-2;i+3){
        primerSuma = primerSuma + fx[i];
    }

    primerSuma = 3*primerSuma;

    double segundaSuma = 0;

    for(int i = 2; i<n-1;i+3){
        segundaSuma = segundaSuma + fx[i];
    }

    segundaSuma = segundaSuma*3;

    double terceraSuma = 0;

    for(int i = 3; i<n-3;i+3){
        terceraSuma = terceraSuma + fx[i];
    }

    double sum = fx[0] + primerSuma + segundaSuma + terceraSuma + fx[n-1];

    double res = (3*h/8)*sum;
    cout << res << endl;
}

int main(int argc, char *argv[]) {
    enablePres();
    //dosPuntos(-11.8646, -12.4775, 0.1);1
    //tresPuntosProgresiva(-13.0891, -13.6997, -14.3092,0.1);
    //tresPuntosRegresiva(-13.0891, -12.4775, -11.8646,0.1);
    //tresPuntosCentrada(-12.4775,-13.6997,0.1);
    //cincoPuntosCentrada(-11.8646, -12.4775, -13.6997, -14.3092, 0.1);
    //cincoPuntos(-11.8646, -12.4775, -13.0891, -13.6997, -14.3092, 0.1);
    trapecio(5,7.4,-3.48366e-11,7.53488e-24);
    trapecioGeneralizado({5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4},
                         {-3.48366e-11, -2.42198e-12, -3.66279e-14, 2.3182e-14, 4.96495e-15, 6.87093e-16,
                          7.59653e-17, 7.1276e18, 5.81955e-19, 4.1847e-20, 2.66514e-21, 1.50582e-22, 7.53488}, 0.1);
    simpsonunterc();
    simpsonGeneralizado();
    simpsontresOct();
    simpsontresOctGeneralizado();
    return 0;
}