#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

int
func (double t, const double y[], double f[],
      void *params1, void *params2, void *params3)
{
  (void)(t); /* avoid unused parameter warning */
  float E = *(double *)params1;
  float V = *(double *)params2;
  float z = *(double *)params3;
  f[0] = y[1];
  f[1] = -((2*pow(1, -30))/pow(1.05, -34))*(E-V)*y[0];
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

int
main (void)
{
  float E;
  float V; //potentiel*
  float z; //valeut initial de Y2 en 0
  printf("Entrer la valeur de l énergie:\n");
  scanf("%d", &E);
  printf("Entrer la valeur du potentiel:\n");
  scanf("%d", &V);
  printf("Entrer la valeur de z:\n");
  scanf("%d", &z);
  gsl_odeiv2_system sys = {func, jac, 2, &E, &V, &z};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 1.0, 0.0 };

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}
