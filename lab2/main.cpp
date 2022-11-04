#include <iostream>
#include <stdio.h>
#include <cmath>
#define PI 3.1415926
#define ABS(x) ((x)>0?(x):-(x))

// u(x,t)
// du/dt - po(x) * d2u/dx2 = f(x,t)
// u(x,0) = phi(x)
// du(0,t)/dx = 0, u(1,t) = 0

class Equation {

private:

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

  const double a = 0.0, b = 1.0, T = 10.0;
  int N, M;
  double h, t;
  double* xi, * tn;

  double** u;

  Equation(int M, int N) {
    this->N = N;
    this->M = M;
    h = (b - a) / M;
    t = T / N;

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

  }

  double solution(double x, double t) {
    return (2.0 * x * exp(-2.0 * t) + 1) * (x * sin(PI * x));
  }

  ~Equation() {
    delete[] xi;
    delete[] tn;
  }

  void ForwardMethod() {

    for (int n = 0; n <= N - 1; ++n) {
      for (int i = 1; i <= M - 1; ++i) {
        u[n + 1][i] = t * po(xi[i]) / (h * h) * \
          (u[n][i + 1] - 2.0 * u[n][i] + u[n][i - 1]) + u[n][i] + t * f(tn[n], xi[i]);
      }

      //find u[n+1][0] using boundary condition
      u[n + 1][0] = (4.0 * u[n][1] - u[n][2]) / 3.0;
    }

  }

  void BackwardMethod() {

    double* F = new double[M + 1];
    double* C = new double[M + 1];
    double* D = new double[M + 1];
    double* B = new double[M + 1];

    double poi;
    int i;

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

    int i;
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

int main(void) {

  int M = 10, N = 500;
  Equation eq(M, N);

  //eq.ForwardMethod();
  //eq.BackwardMethod();
  //eq.SymmetricMethod();

  printf("x\t\texact\t\tnumeric\t\terr\n");
  for (int n = 0; n <= N; ++n) {
    printf("t = %lf\n", eq.tn[n]);
    for (int i = 0; i <= M; ++i) {
      printf("%lf\t%lf\t%lf\t%g\n", eq.xi[i], eq.solution(eq.xi[i], eq.tn[n]), eq.u[n][i], ABS(eq.solution(eq.xi[i], eq.tn[n]) - eq.u[n][i]));
    }
    printf("\n");
  }

  printf("end");
  return 0;
}