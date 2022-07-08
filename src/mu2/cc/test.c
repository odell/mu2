#include <stdio.h>
#include "scatter.h"

int main(void) {

    int i, n=200, m=30;
    double R=5.5;
    double qmax=10*2/R;
    double mass=0.12400095845605276;

    FILE *fp;
    double k[m];
    double q[n];
    double wq[n];
    double v[n*n];

    fp = fopen("../datfiles/k.txt", "r");
    for (i = 0; i < m; i++) {
        fscanf(fp, "%lf\n", &k[i]);
    }
    fclose(fp);

    fp = fopen("../datfiles/q.txt", "r");
    for (i = 0; i < n; i++) {
        fscanf(fp, "%lf\n", &q[i]);
    }
    fclose(fp);

    fp = fopen("../datfiles/wq.txt", "r");
    for (int i = 0; i < n; i++) {
        fscanf(fp, "%lf\n", &wq[i]);
    }
    fclose(fp);

    fp = fopen("../datfiles/v_tilde.txt", "r");
    for (int i = 0; i < n*n; i++) {
        fscanf(fp, "%lf\n", &v[i]);
    }
    fclose(fp);

    for (i = 0; i < m; i++) {
        printf("%.8lf\n", k3cotdelta(k[i], v, q, wq, n, qmax, mass));
    }
    return 0;
}