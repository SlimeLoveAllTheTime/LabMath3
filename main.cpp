#include <iostream>
#include "rkf45.h"
#include "FORSYTHE.h"

// количество вычислений
int n = 11;

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
        if (i % 1 == 0) printf("\nY0[0] = %.15lf     Y0[1] = %.15lf", y[i], dy[i]);
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
    double RE = 1e-3, AE = 1e-3;
    printf("RKF45");
    for (int i = 0; i < n; ++i) {
        RKF45(fun, 2, Y0, &T, &TOUT, &RE, &AE, &iflag, work, iwork);
        if (i % 1 == 0) printf("\nY[0] = %.15lf     Y0[1] = %.15lf     iflag = %d", Y0[0], Y0[1], iflag);
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

    k[0] = Xn[0] + 3 * k2[0] / 4;
    k[1] = Xn[1] + 3 * k2[1] / 4;

    fun(t + 3 * h / 4, k, k3);

    k3[0] = h * k3[0];
    k3[1] = h * k3[1];

    Xnp0[i] = Xn[0] + (2 * k1[0] + 3 * k2[0] + 4 * k3[0]) / 9;
    Xnp1[i] = Xn[1] + (2 * k1[1] + 3 * k2[1] + 4 * k3[1]) / 9;
}

// применение метода Рунге-Кутты 3ей степени точности на всех точка в промежутке [1, 2]
void RungeKutta(double *Xn, double *Xnp0, double *Xnp1, double h) {
    double T = 1.0;
    printf("\nRungeKutta3\n");
    for (int i = 1; i < 12; ++i) {
        RungeKutta3(T, Xn, Xnp0, Xnp1, h, i);
        if ((i - 1) % 1 == 0) printf("Y[0] = %.15lf     Y0[1] = %.15lf\n", Xn[0], Xn[1]);
        Xn[0] = Xnp0[i];
        Xn[1] = Xnp1[i];
        T += h;
    }
}

// вычисление глобальной погрешности работы процедуры rkf45
void globalErrorRKF(double *Y0rkf, double *Y1rkf, double *F0, double *F1) {
    printf("\nRKF45 GLOBAL ERROR:\n");
    for (int i = 0; i < n; ++i) {
        if (i % 1 == 0) printf("Y0rkf = %.15lf     Y1rkf = %.15lf\n", abs(Y0rkf[i] - F0[i]), abs(Y1rkf[i] - F1[i]));
    }
}

// вычисление глобальной погрешности работы процедуры RungeKutta
void globalErrorRungeKutta(double *Y0rk3, double *Y1rk3, double *F0, double *F1) {
    printf("\nRunge Kutta 3 GLOBAL ERROR:\n");
    for (int i = 0; i < n; ++i) {
        if (i % 1 == 0) printf("Y0rk3 = %.15lf     Y1rk3 = %.15lf\n", abs(Y0rk3[i] - F0[i]), abs(Y1rk3[i] - F1[i]));
    }
}

// вычисление локальной погрешности RKF45
void localErrorRKF(double *Y0rkf, double *Y1rkf,
                   double *F0, double *F1,
                   double *localError0, double *localError1,
                   double h) {
    double TOUT = 1.0;
    printf("\nRKF45 LOCAL ERROR: h = %.4lf\n", h);
    for (int i = 0; i < n; i++) {
        localError0[i] = abs(Y0rkf[i] - F0[i]);
        localError1[i] = abs(Y1rkf[i] - F1[i]);
        for (int j = i - 1; j > -1; --j) {
            if (i == 0) continue;
            localError0[i] -= localError0[j];
            localError1[i] -= localError1[j];
        }
    }
    for (int i = 0; i < n; ++i) {
        if (i % 1 == 0) {
            printf("Tout = %lf     ", TOUT);
            printf("Local Error:   RKF45[0] = %.15lf     RKF45[1] = %.15lf\n", localError0[i], localError1[i]);
            TOUT += h;
        }
    }
}

// вычисление локальной погрешности Runge Kutta 3
void localErrorRungeKutta(double *Xnp0, double *Xnp1,
                          double *F0, double *F1,
                          double *localError0, double *localError1,
                          double h) {
    double TOUT = 1.0;
    printf("\nRunge Kutta 3 LOCAL ERROR: h = %.4lf\n", h);
    for (int i = 0; i < n; i++) {
        localError0[i] = abs(Xnp0[i] - F0[i]);
        localError1[i] = abs(Xnp1[i] - F1[i]);
        for (int j = i - 1; j > -1; --j) {
            if (i == 0) continue;
            localError0[i] -= localError0[j];
            localError1[i] -= localError1[j];
        }
    }
    for (int i = 0; i < n; ++i) {
        if (i % 1 == 0) {
            printf("Tout = %lf     ", TOUT);
            printf("Local Error:   RungeKutta3[0] = %.15lf     RungeKutta3[1] = %.15lf\n", localError0[i],
                   localError1[i]);
            TOUT += h;
        }
    }
}

int main() {

    double h = 0.1;

    double Y0rkf[81], Y1rkf[81], localErrorRKF0[81], localErrorRKF1[81];
    double Y0[2] = { exp(2), 2 * exp(2) };
    rkf45(Y0rkf, Y1rkf, Y0[0], Y0[1], h);

    double Xn1[]= { exp(2), 2 * exp(2) }, Xnp0[84] = { exp(2) }, Xnp1[84] { 2 * exp(2) };
    double localErrorRK0[81], localErrorRK1[81];
    RungeKutta(Xn1, Xnp0, Xnp1, h);

    double F0[82], F1[82];
    funSolver(F0, F1, h);

    localErrorRKF(Y0rkf, Y1rkf, F0, F1, localErrorRKF0, localErrorRKF1, h);
    localErrorRungeKutta(Xnp0, Xnp1, F0, F1, localErrorRK0, localErrorRK1, h);
    globalErrorRKF(Y0rkf, Y1rkf, F0, F1);
    globalErrorRungeKutta(Xnp0, Xnp1, F0, F1);

    return 0;
}