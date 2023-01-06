#include <stdio.h>
#include "scatter_pert.h"

int main(int argc, char** argv) {

    int n = 200;
    double qmax = 7.218587863749154;
    double mass = 0.08251315259652389;
    double beta6 = 5.54125;
    double kmin = 0.01/beta6;
    double kmax = 0.3/beta6;
    int m = 30;
    double dk = (kmax - kmin) / (m - 1.0);

    FILE *fp;
    double k[m];
    double q[n];
    double wq[n];
    double v0[n*n];
    double v1[n*n];

    for (int i = 0; i < m; i++)
        k[i] = kmin + i*dk;

    fp = fopen("/Users/danielodell/vdw-nlo/datfiles/q.txt", "r");
    for (int i = 0; i < n; i++) {
        fscanf(fp, "%lf\n", &q[i]);
    }
    fclose(fp);

    fp = fopen("/Users/danielodell/vdw-nlo/datfiles/wq.txt", "r");
    for (int i = 0; i < n; i++) {
        fscanf(fp, "%lf\n", &wq[i]);
    }
    fclose(fp);

    fp = fopen("/Users/danielodell/vdw-nlo/datfiles/v0.txt", "r");
    for (int i = 0; i < n*n; i++) {
        fscanf(fp, "%lf\n", &v0[i]);
    }
    fclose(fp);

    fp = fopen("/Users/danielodell/vdw-nlo/datfiles/v1.txt", "r");
    for (int i = 0; i < n*n; i++) {
        fscanf(fp, "%lf\n", &v1[i]);
    }
    fclose(fp);

    // driving term

    for (int i = 0; i < m; i++) {
        gsl_complex t = t_on_shell_pert1(k[i], v0, v1, q, wq, n, qmax, mass);
        // printf("%.8lf %.8lf\n", GSL_REAL(t), GSL_IMAG(t));
        double kcd = kcotdelta_pert1(k[i], v0, v1, q, wq, n, qmax, 0, mass);
        printf("%.8e  %.8e  %.8e\n", GSL_REAL(t), GSL_IMAG(t), kcd);
    }
    return 0;

}