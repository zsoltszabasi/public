/*****************************************************************************
*                                                                            *
* Created by zsolt.szabasi                                                   *
*                                                                            *
******************************************************************************/

#include <math.h>

int gauss_jordan(double* m, int w, int h) { /* Convert matrix to Reduced Row Echelon Form. */
    int y, c, c2;
    double* py = m;
    for (y = 0; y < h; y++) {
        double* pm = py;
        double* py2 = py + w;
        double* p1;
        double* p2;
        double f;
        c = h - y - 1;
        while (c-- > 0) { /* find max pivot */
            if (fabs(py2[y]) > fabs(pm[y])) {
                pm = py2;
            }
            py2 += w;
        }
        p1 = py;
        p2 = pm;
        c = w;
        while (c-- > 0) {
            double t = *p1;
            *p1++ = *p2;
            *p2++ = t;
        }
        if (fabs(py[y]) <= 1e-15) {
            return 0; /* singular matrix */
        }
        f = 1 / py[y];
        py2 = py + w + y;
        c2 = h - y - 1;
        while (c2-- > 0) { /* eliminate column y */
            double f1 = *py2 * f;
            double* p1 = py2;
            double* p2 = py + y;
            c = w - y;
            while (c-- > 0) {
                *p1++ -= *p2++ * f1;
            }
            py2 += w;
        }
        py += w;
    }
    py = m + (h - 1) * w;
    for (y = h - 1; y >= 0; y--) { /* subst */
        double f = 1 / py[y];
        double* py2 = m;
        double* p1;
        c2 = y;
        while (c2-- > 0) {
            double* p1 = py2 + w - 1;
            double* p2 = py + w - 1;
            double f1 = py2[y] * f;
            c = w - y;
            while (c-- > 0) {
                *p1-- -= *p2-- * f1;
            }
            py2 += w;
        }
        py[y] *= f;
        p1 = py + h;
        c = w - h;
        while (c-- > 0) {
            *p1++ *= f; /* normalize row y */
        }
        py -= w;
    }
    return 1;
}
