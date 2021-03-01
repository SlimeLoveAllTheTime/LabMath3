#include <iostream>
#include "rkf45.h"
#include "FORSYTHE.h"

// количество вычислений
int n = 81;

// система первого порядка
void fun(double t, double *y, double *dy) {
    dy[0] = y[1];
    dy[1] = ((t + 1) * y[1] + 2 * (t - 1) * y[0]) / t;
}

// вычисление точного значения системы
void funSolver(double *y, double *dy, double h) {
    double T = 1.0;
    printf("\nfun solutions");
    for (int i = 0; i < n; ++i) {
        y[i] = exp(2 * T);
        dy[i] = 2 * exp(2 * T);
        if (i % 8 == 0) printf("\nY0[0] = %.10lf     Y0[1] = %.10lf", y[i], dy[i]);
        T += h;
    }
}

// применение процедуры rkf45 на всех точка в промежутке [1, 2]
void rkf45(double *y, double *dy, double y0, double dy0, double h) {
    double Y0[] = { y0, dy0 };
    double T = 1.0, TOUT = 1.0;
    int iflag = 1;
    int iwork[30];
    double work[15];
    double RE = 1e-5, AE = 1e-5;
    printf("RKF45");
    for (int i = 0; i < n; ++i) {
        RKF45(fun, 2, Y0, &T, &TOUT, &RE, &AE, &iflag, work, iwork);
        if (i % 8 == 0) printf("\nY[0] = %.10lf     Y0[1] = %.10lf     iflag = %d", Y0[0], Y0[1], iflag);
        y[i] = Y0[0];
        dy[i] = Y0[1];
        TOUT += h;
    }
}

// применение метода Рунге-Кутты 3ей степени точности
void RungeKutta3(double t, double *Xn, double *Xnp0, double *Xnp1, double h, int i) {

    double k[2], k2[2], k1[2], k3[2];

    fun(t, Xn, k1);

    k1[0] = h * k1[0];
    k1[1] = h * k1[1];

    k[0] = Xn[0] + k1[0] / 2;
    k[1] = Xn[1] + k1[1] / 2;

    fun(t + h / 2, k, k2);

    k2[0] = h * k2[0];
    k2[1] = h * k2[1];

    k[0] = Xn[0] - k1[0] + 2 * k2[0];
    k[1] = Xn[1] - k1[1] + 2 * k2[1];

    fun(t + h, k, k3);

    k3[0] = h * k3[0];
    k3[1] = h * k3[1];

    Xnp0[i] = Xn[0] + (k1[0] + 4 * k2[0] + k3[0]) / 6;
    Xnp1[i] = Xn[1] + (k1[1] + 4 * k2[1] + k3[1]) / 6;
}

// применение метода Рунге-Кутты 3ей степени точности на всех точка в промежутке [1, 2]
void RungeKutta(double *Xn, double *Xnp0, double *Xnp1, double h) {
    double T = 1.0;
    printf("\nRungeKutta3\n");
    for (int i = 1; i < 82; ++i) {
        RungeKutta3(T, Xn, Xnp0, Xnp1, h, i);
        if ((i - 1) % 8 == 0) printf("Y[0] = %.10lf     Y0[1] = %.10lf\n", Xn[0], Xn[1]);
        Xn[0] = Xnp0[i];
        Xn[1] = Xnp1[i];
        T += h;
    }
}

// вычисление погрешности работы процедуры rkf45
void errorRKF(double *Y0rkf, double *Y1rkf, double *F0, double *F1) {
    printf("\nRKF45 ERROR:\n");
    for (int i = 0; i < n; ++i) {
        if (i % 8 == 0) printf("Y0rkf = %.10lf     Y1rkf = %.10lf\n", abs(Y0rkf[i] - F0[i]), abs(Y1rkf[i] - F1[i]));
    }
}

// вычисление погрешности работы процедуры RungeKutta
void errorRungeKutta(double *Y0rk3, double *Y1rk3, double *F0, double *F1) {
    printf("\nRunge Kutta 3 ERROR:\n");
    for (int i = 0; i < n; ++i) {
        if (i % 8 == 0) printf("Y0rk3 = %.10lf     Y1rk3 = %.10lf\n", abs(Y0rk3[i] - F0[i]), abs(Y1rk3[i] - F1[i]));
    }
}

int main() {

    double h = 0.0125;

    double Y0rkf[81], Y1rkf[81];
    double Y0[2] = { exp(2), 2 * exp(2) };
    rkf45(Y0rkf, Y1rkf, Y0[0], Y0[1], h);

    double Xn1[]= { exp(2), 2 * exp(2) }, Xnp0[84] = { exp(2) }, Xnp1[84] { 2 * exp(2) };
    RungeKutta(Xn1, Xnp0, Xnp1, h);

    double F0[82], F1[82];
    funSolver(F0, F1, h);

    errorRKF(Y0rkf, Y1rkf, F0, F1);
    errorRungeKutta(Xnp0, Xnp1, F0, F1);

    return 0;
}