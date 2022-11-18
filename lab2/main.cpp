#include <iostream>
#include <stdio.h>
#include <cmath>
#include <windows.h>
#define PI 3.1415926
#define RMAX 2.0 // max (1+x) as 0<=x<=1
#define ABS(x) ((x)>0?(x):-(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

// u(x,t), 0<x<1
// du/dt - po(x) * d2u/dx2 = f(x,t)
// u(x,0) = phi(x)
// du(0,t)/dx = 0, u(1,t) = 0

class Equation {

private:

  double unorm;

  //for backward and symmetric methods
  void TridiagonalMatrixAlgorithm(double* B, double* C, double* D, double* F, int n) {
    int i;
    double* lambda = new double[M + 1];
    double* delta = new double[M + 1];
    double* y = new double[M + 1];

    //прямой ход
    delta[0] = -D[0] / C[0];
    lambda[0] = F[0] / C[0];
    for (i = 1; i <= M - 1; i++) {
      delta[i] = -D[i] / (B[i] * delta[i - 1] + C[i]);
      lambda[i] = (F[i] - B[i] * lambda[i - 1]) / (B[i] * delta[i - 1] + C[i]);
    }
    delta[M] = 0;
    lambda[M] = (F[M] - B[M] * lambda[M - 1]) / (B[M] * delta[M - 1] + C[M]);

    //обратный ход
    y[M] = lambda[M];
    for (i = M - 1; i >= 0; i--) {
      y[i] = delta[i] * y[i + 1] + lambda[i];

      //write answer
      u[n + 1][i] = y[i];
    }

    delete[] lambda;
    delete[] delta;
    delete[] y;
  }

  double po(double x) {
    return 1.0 + x;
  }

  double phi(double x) {
    return x * (2.0 * x + 1.0) * sin(PI * x);
  }

  double f(double t, double x) {
    return -4.0 * exp(-2.0 * t) * x * x * sin(PI * x) - exp(-2.0 * t) * (x + 1.0) * \
      (2.0 * PI * (exp(2.0 * t) + 4.0 * x) * cos(PI * x) - \
        (PI * PI * exp(2.0 * t) * x + 2.0 * PI * PI * x * x - 4.0) * sin(PI * x));
  }

public:

  const double a = 0.0, b = 1.0;
  int N, M;
  double h, t, T;
  double* xi, * tn;

  double** u;

  Equation(int M, double T, double t) {
    this->T = T;
    this->M = M;
    h = (b - a) / M;

    this->t = h * h / (2.0 * RMAX);
    //this->t = 0.2; //backward
    //this->t = 0.8; //symmetric

    this->t += t; //to show unsteadiness

    N = (int)(T / this->t);

    xi = new double[M + 1];
    tn = new double[N + 1];

    for (int i = 0; i <= M; ++i)
      xi[i] = a + i * h;
    for (int n = 0; n <= N; ++n)
      tn[n] = n * t;

    u = new double* [N + 1];
    for (int j = 0; j <= N; ++j)
      u[j] = new double[M + 1];

    //init u(x,t) using initial and boundary conditions
    for (int i = 0; i <= M; ++i)
      u[0][i] = phi(xi[i]);

    for (int n = 0; n <= N; ++n) {
      //u[n][0] must be found
      u[n][M] = 0.0;
    }

    //find norm of exact solution
    double hnorm = 0.0, norm = 0.0;
    for (int n = 0; n <= N; ++n) {
      for (int i = 1; i <= M; ++i) {
        if (hnorm < ABS(solution(xi[i], tn[n])))
          hnorm = ABS(solution(xi[i], tn[n]));
      }
      if (norm < hnorm)
        norm = hnorm;
    }
    unorm = norm;
  }

  double solution(double x, double t) {
    return (2.0 * x * exp(-2.0 * t) + 1) * (x * sin(PI * x));
  }

  //for numeric solution
  double htNorm() {

    double hnorm = 0.0, norm = 0.0, val, diff;

    for (int n = 0; n <= N; ++n) {
      for (int i = 0; i <= M; ++i) {
        val = solution(xi[i], tn[n]);
        diff = ABS(u[n][i] - val);
        if (hnorm < diff)
          hnorm = diff;
      }
      if (norm < hnorm)
        norm = hnorm;
    }

    return norm / unorm;
  }

  double hNorm(int n) {
    double hnorm = 0.0, val, diff;
    for (int i = 0; i <= M; ++i) {
      val = solution(xi[i], tn[n]);
      diff = ABS(u[n][i] - val);

      if (hnorm < diff)
        hnorm = diff;
    }
    return hnorm;
  }

  ~Equation() {
    delete[] xi;
    delete[] tn;
    for (int n = 0; n <= N; ++n)
      delete[] u[n];
    delete[] u;
  }

  void ForwardMethod() {

    //for (int n = 0; n <= N; ++n)
      //printf("%lf ", n * t);
    //printf("\n\n\n\n");

    for (int n = 0; n <= N - 1; ++n) {
      for (int i = 1; i <= M - 1; ++i) {
        u[n + 1][i] = t * po(xi[i]) / (h * h) * \
          (u[n][i + 1] - 2.0 * u[n][i] + u[n][i - 1]) + u[n][i] + t * f(tn[n], xi[i]);
      }

      //find u[n+1][0] using boundary condition
      u[n + 1][0] = (4.0 * u[n][1] - u[n][2]) / 3.0;

      //show norm
      //printf("%lf ", hNorm(n));

    }
    //printf("%lf ", hNorm(N));
    //printf("\n");

  }

  void BackwardMethod() {

    double* F = new double[M + 1];
    double* C = new double[M + 1];
    double* D = new double[M + 1];
    double* B = new double[M + 1];

    double poi;

    for (int n = 0; n <= N - 1; ++n) {

      poi = po(xi[0]);
      B[0] = 0.0;
      C[0] = 2.0 * poi * t;
      D[0] = -2.0 * poi * t + h * h;
      F[0] = h * h * (u[n][1] + t * f(tn[n + 1], xi[1]));

      for (int i = 1; i <= M - 1; ++i) {
        poi = po(xi[i]);

        B[i] = -poi * t;
        C[i] = h * h + 2.0 * poi * t;
        D[i] = -poi * t;
        F[i] = h * h * (u[n][i] + t * f(tn[n], xi[i]));
      }
      B[M] = 0.0;
      C[M] = 1.0;
      D[M] = 0.0;
      F[M] = 0.0;

      TridiagonalMatrixAlgorithm(B, C, D, F, n);
    }

    delete[] F;
    delete[] C;
    delete[] D;
    delete[] B;
  }

  void SymmetricMethod() {

    double* F = new double[M + 1];
    double* C = new double[M + 1];
    double* D = new double[M + 1];
    double* B = new double[M + 1];

    double poi;

    for (int n = 0; n <= N - 1; ++n) {

      poi = po(xi[0]);
      B[0] = 0.0;
      C[0] = 2.0 * poi * t;
      D[0] = 2.0 * (h * h - poi * t);
      F[0] = u[n][2] * poi * t + u[n][1] * 2.0 * (-poi * t + h * h) + u[n][0] * poi * t + \
        2.0 * h * h * t * (f(tn[n] + t / 2.0, xi[1]));

      for (int i = 1; i <= M - 1; ++i) {
        poi = po(xi[i]);

        B[i] = -poi * t;
        C[i] = 2.0 * (h * h + poi * t);
        D[i] = -poi * t;
        F[i] = u[n][i + 1] * poi * t + u[n][i] * 2.0 * (-poi * t + h * h) + u[n][i - 1] * poi * t + \
          2.0 * h * h * t * (f(tn[n] + t / 2.0, xi[i]));
      }
      B[M] = 0.0;
      C[M] = 1.0;
      D[M] = 0.0;
      F[M] = 0.0;

      TridiagonalMatrixAlgorithm(B, C, D, F, n);
    }

    delete[] F;
    delete[] C;
    delete[] D;
    delete[] B;
  }

};

double PCFreq = 0.0;
__int64 CounterStart = 0;

void StartCounter()
{
  LARGE_INTEGER li;
  QueryPerformanceFrequency(&li);

  PCFreq = double(li.QuadPart) / 1000.0;

  QueryPerformanceCounter(&li);
  CounterStart = li.QuadPart;
}

double GetCounter()
{
  LARGE_INTEGER li;
  QueryPerformanceCounter(&li);
  return double(li.QuadPart - CounterStart) / PCFreq;
}

int main(void) {

  // 1) forward method
  int M = 1;
  const double T = 3.5;
  Equation* eq;
  double eps = 0.16, norm;



  /*
  // 1)
  eq = new Equation(M,T, 0.0);
  do {
    M+=5;
    delete eq;
    eq = new Equation(M, T, 0.0);
    eq->ForwardMethod();
    norm = eq->htNorm();
    printf("M = %i norm = %lf N = %i\n", M, norm, eq->N);
  } while (norm > eps);
  printf("end\n");
  */

  /*
  M = 21;
  eq = new Equation(M, T, 7.75e-5);
  eq->ForwardMethod();
  */

  //2)
  /*
  M = 1;
  eq = new Equation(M, T, 0.0);
  do {
    M+=5;
    delete eq;
    eq = new Equation(M, T, 0.0);
    eq->BackwardMethod();
    norm = eq->htNorm();
    printf("M = %i norm = %lf N = %i\n", M, norm, eq->N);
  } while (norm > eps);
  printf("end\n");
  */

  /*
  M = 21;
  double t = 0.0, step = 5e-3;
  eq = new Equation(M, T, t);
  eq->BackwardMethod();
  do {
    ++M;
    delete eq;
    t += step;
    eq = new Equation(M, T, t);
    eq->BackwardMethod();
    printf("%lf ", eq->t);
    //printf("%lf ", eq->htNorm());
  } while (eq->htNorm() < eps);
  */

  //3)
  /*
  M = 1;
  eq = new Equation(M, T, 0.0);
  do {
    M+=20;
    delete eq;
    eq = new Equation(M, T, 0.0);
    eq->SymmetricMethod();
    norm = eq->htNorm();
    printf("M = %i norm = %lf\n", M, norm);
  } while (norm > eps);
  printf("end\n");
  */

  /*
  eps = 0.16;
  M = 21;
  double t = 0.0, step = 2e-2;
  eq = new Equation(M, T, t);
  do {
    delete eq;
    t += step;
    eq = new Equation(M, T, t);
    eq->SymmetricMethod();
    printf("%lf ", eq->t);
    printf("%lf ", eq->htNorm());
  } while (eq->htNorm() < eps);

  double t0 = eq->t - step;
  delete eq;
  eq = new Equation(M, T, t0);
  eq->SymmetricMethod();
  printf("\n\n\noptimal: ");
  printf("%lf ", eq->t);
  printf("%lf ", eq->htNorm());
  delete eq;

  eq = new Equation(M, T, t0);
  eq->SymmetricMethod();
  printf("err = %lf ", eq->htNorm());
  printf("t0 = %lf", t0);
  */

  //4)
  /*
  eps = 0.21;
  eq = new Equation(M, T, 0.0);
  do {
    M = 1;
    eps -= 0.01;

    StartCounter();
    do {
      M += 20;
      delete eq;
      eq = new Equation(M, T, 0.0);
      eq->ForwardMethod();
    } while (eq->htNorm() > eps);
    printf("%lf ", GetCounter());

  } while (eps > 0.13);
  */

  /*
  eps = 0.21;
  eq = new Equation(M, T, 0.0);
  do {
    M = 1;
    eps -= 0.01;
    StartCounter();
    do {
      M += 20;
      delete eq;
      eq = new Equation(M, T, 0.0);
      eq->BackwardMethod();
    } while (eq->htNorm() > eps);
    printf("%lf ", GetCounter());

  } while (eps > 0.13);
  */

  //double t0 = 0.820567;
  double t0 = 0.8;
  eps = 0.22;
  eq = new Equation(M, T, 0.0);
  do {
    M = 1;
    eps -= 0.01;
    StartCounter();
    do {
      M += 10;
      delete eq;
      eq = new Equation(M, T, t0);
      eq->SymmetricMethod();
    } while (eq->htNorm() > eps);
    printf("%lf ", GetCounter());

  } while (eps > 0.16);



  // M = 21;
   //eq = new Equation(M, T, 0.0);
   //eq->SymmetricMethod();
   //printf("%lf", eq->htNorm());
  delete eq;
  return 0;
}