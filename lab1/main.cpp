#include <iostream>
#include <cmath>
#define ABS(x) ((x)>0?(x):-(x))
#define PI 3.1415926

class Equation {
public:
  //[0, 1]
  const double a = 0.0, b = 1.0;
  double* xk;
  int n;
  double* yk;

  Equation(int n) {
    double x = a;
    this->n = n;
    xk = new double[n];
    double h = (b - a) / n;
    for (int i = 0; i < n; ++i, x += h)
      xk[i] = x;

    //yk = new double[n];
  }
  ~Equation() {
    delete[] xk;
    delete[] yk;
  }

  double p(double x) {
    //return 1.0;
    return x;
  }

  double q(double x) {
    //return 1.0;
    return x + 1.0;
  }

  double f(double x) {
    return x * (x + 1.0) * sin(PI * x) + x * (sin(PI * x) + PI * x * cos(PI * x)) \
      - PI * (2.0 * cos(PI * x) - PI * x * sin(PI * x));
    //return x * x + 2.0 * x - 3.0;
  }

  double exact(double x) {
    return x * sin(PI * x);
    //return x*x-1.0;
  }

  double SecondNorm() {
    double ans = 0.0;
    double x = a, h = (b - a) / n;
    for (int i = 0; i <= n; ++i, x += h)
      ans += (yk[i] - exact(x))* (yk[i] - exact(x));
    ans *= h;
    return sqrt(ans);
  }
};

void FDlin(Equation* eq) {
  int n = eq->n;
  double x;
  int k;
  double* F, * C, * D, * B, * lambda, * delta, * y;
  double pk, qk, fk;
  double a = eq->a, b = eq->b;
  double h = (b - a) / n;

  y = new double[n + 1];
  F = new double[n + 1];
  C = new double[n + 1];
  D = new double[n + 1];
  B = new double[n + 1];
  lambda = new double[n + 1];
  delta = new double[n + 1];

  double p1 = eq->p(a + h), q1 = eq->q(a + h), f1 = eq->f(a + h);

  //right side derivative
  //y0' = 0
  B[0] = 0.0;
  C[0] = -1.0;
  D[0] = 1.0;
  F[0] = 0.0;

  x = a + h;
  h = (b - a) / n;
  for (k = 1; k < n; k++, x += h) {
    pk = eq->p(x);
    qk = eq->q(x);
    fk = eq->f(x);

    B[k] = -1.0;
    C[k] = 2.0 + qk * h * h - pk*h;
    D[k] = pk * h - 1.0;
    F[k] = fk * h * h;
  }

  //yn = 0
  B[n] = 0.0;
  C[n] = 1.0;
  D[n] = 0.0;
  F[n] = 0.0;

  //Метод прогонки

  //прямой ход
  delta[0] = -D[0] / C[0];
  lambda[0] = F[0] / C[0];
  for (k = 1; k < n; k++) {
    delta[k] = -D[k] / (B[k] * delta[k - 1] + C[k]);
    lambda[k] = (F[k] - B[k] * lambda[k - 1]) / (B[k] * delta[k - 1] + C[k]);
  }
  delta[n] = 0;
  lambda[n] = (F[n] - B[n] * lambda[n - 1]) / (B[n] * delta[n - 1] + C[n]);

  //обратный ход
  y[n] = lambda[n];
  for (k = n - 1; k >= 0; k--)
    y[k] = delta[k] * y[k + 1] + lambda[k];

  delete[] F;
  delete[] C;
  delete[] D;
  delete[] B;
  delete[] lambda;
  delete[] delta;

  eq->yk = y;
}

void FDquad(Equation* eq) {
  int n = eq->n;
  double x;
  int k;
  double* F, * C, * D, * B, * lambda, * delta, * y;
  double pk, qk, fk;
  double a = eq->a, b = eq->b;
  double h = (b - a) / n;

  y = new double[n + 1];
  F = new double[n + 1];
  C = new double[n + 1];
  D = new double[n + 1];
  B = new double[n + 1];
  lambda = new double[n + 1];
  delta = new double[n + 1];

  //заполняем трехдиагональную матрицу коэффициентов системы
  double p1 = eq->p(a + h), q1 = eq->q(a + h), f1 = eq->f(a + h);

  //y0' = 0
  B[0] = 0;
  C[0] = -3.0 + 2.0 / (p1*h-2.0)*(1-p1*h*h/2.0);
  D[0] = 4.0 + 2.0 / (p1*h-2.0)*(2+q1*h*h);
  F[0] = 2.0*f1*h*h/(p1*h-2.0);

  x = a+h;
  h = (b - a) / n;
  for (k = 1; k < n; k++, x+=h) {
    pk = eq->p(x);
    qk = eq->q(x);
    fk = eq->f(x);

    B[k] = - 1.0 - pk * h / 2.0; //y_(k-1)
    C[k] = 2.0 + qk * h * h; //y_k
    D[k] = pk * h / 2.0 - 1.0; //y_(k+1)
    F[k] = fk * h * h; //f_i
  }

  //double pn_1 = eq->p(b - h), qn_1 = eq->p(b - h), fn_1 = eq->p(b - h);
  //yn = 0
  B[n] = 0.0;
  C[n] = 1.0;
  D[n] = 0.0;
  F[n] = 0.0;

  //Метод прогонки

  //прямой ход
  delta[0] = -D[0] / C[0];
  lambda[0] = F[0] / C[0];
  for (k = 1; k < n; k++) {
    delta[k] = -D[k] / (B[k] * delta[k - 1] + C[k]);
    lambda[k] = (F[k] - B[k] * lambda[k - 1]) / (B[k] * delta[k - 1] + C[k]);
  }
  delta[n] = 0;
  lambda[n] = (F[n] - B[n] * lambda[n - 1]) / (B[n] * delta[n - 1] + C[n]);

  //обратный ход
  y[n] = lambda[n];
  for (k = n - 1; k >= 0; k--)
    y[k] = delta[k] * y[k + 1] + lambda[k];

  delete[] F;
  delete[] C;
  delete[] D;
  delete[] B;
  delete[] lambda;
  delete[] delta;

  eq->yk = y;
}

int main() {
  int n;
  double x;
  double h;

  //1) exact function
  /*
  h = 1e-2;
  for (x = eq.a; x <= eq.b; x += h)
    printf("%lf ", x);
  printf("\n\n");
  for (x = eq.a; x <= eq.b+h; x += h)
    printf("%lf ", eq.exact(x));
  */

  //2)solve diff equation
  void (*solver)(Equation*);
  solver = FDlin;

  double eps = 1e-8;
  for (h = 0.1; h > eps; h /= 2.0)
    printf("%lf ", -log(h));
  printf("\n\n");

  for (h = 0.1; h > eps; h /= 2.0) {
    n = (int)1.0 / h;
    Equation* eq = new Equation(n);
    solver(eq);
    double val = -log(eq->SecondNorm());
    printf("%lf ", val);
    delete eq;
  }
  

  /*
  for (n = 4; n < 10000; n *= 2) {
    h = 1.0 / n;
    printf("n = %i\n", n);
    //printf("numeric\t\t\texact\n\n");
    x = 0.0;
    Equation* eq = new Equation(n);
    FDlin(eq);
    printf("%lf\n", eq->SecondNorm());

    /*
    for (int i = 0; i <= n; ++i, x+=h) {
      //printf("f(%lf) = %lf\tu(%lf) = %lf\t diff = %lf\n", x, eq->yk[i], x, eq->exact(x), eq->SecondNorm());
    }
    
  delete eq;
  printf("\n\n");
  }
  */

  return 0;
}