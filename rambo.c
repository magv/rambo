/* This program integrates phase space integrals in arbitrary
 * number of dimensions by using RAMBO-like mapping from hypercube
 * coordinates into phase space points together with Vegas
 * integration.
 */

#define _POSIX_C_SOURCE 199309L

/* RAMBO IN D DIMENSIONS
 */

#include <alloca.h>
#include <assert.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#include <math.h>

#define M_PI 3.141592653589793238462643383279502884

/* This is the total integral over the phase space of N particles
 * in D dimensions defined as:
 *     \int \prod_i(d^Dp_i \delta(p_i^2)) \delta^D(q - \sum_i p_i),
 *     where q = (w, \bar 0)
 */
double
psvolume(int d, int n, double w)
{
    return
        1.0
        * pow(0.5, n-1)
        * pow(M_PI, (d-2)*(n-1)/2.0)
        * pow(w, (d-2)*n-d)
        * pow(tgamma((d-2)/2.0), n)
        / tgamma((d-2)*(n-1)/2.0)
        / tgamma((d-2)*n/2.0);
}

/* This is the dimensionality of the hypercube used by pspoint().
 */
int
pspoint_dim(int d, int n)
{
    return d*n - (d - 1);
}

/* Uniformly generate a point in d-dimensional n-particle phase
 * space from a point in pspoint_dim(d, n)-dimensional hypercube.
 * Return 0 if successfull.
 */
int
pspoint(int d, int n, double w, const double *uniform, double *point)
{
    int uidx = 0;
    // For particle k={1,...,N} we'll generate uniform momenta
    // direction in D-1 Euclidean space, and then rotate the
    // momenta so that components {k...D-1} will fold into the
    // k-th axis.
    //
    // Such rotation preserves the value of integrals over scalar
    // products of momenta because any configuration of k points
    // can be represented in k-dimensional space. The value of
    // such a rotation is the improved integration convergence
    // with Cuba's Vegas integrator. Otherwise it's optional.
    for (int pidx = 0, axis = 1; pidx < n*d; pidx += d, axis += 1) {
        if (axis == 1) {
            double v0 = gsl_cdf_gamma_Pinv(uniform[uidx++], d-2, 1.0);
            point[0 + 0] = v0;
            point[0 + 1] = v0;
            for (int i = 2; i < d; i++) {
                point[0 + i] = 0.0;
            }
        } else {
            double vv = 0.0;
            for (int i = 1; (i < d) && (i < axis); i++) {
                double vi = gsl_cdf_ugaussian_Pinv(uniform[uidx++]); 
                point[pidx + i] = vi;
                vv += vi*vi;
            }
            if (axis < d) {
                double vv2 = 0.0;
                // This loop can be replaced by a single
                // gsl_cdf_chisq_Pinv, but that is slower and
                // doesn't seem to improve convergence.
                for (int i = axis; i < d; i++) {
                    double vi = gsl_cdf_ugaussian_Pinv(uniform[uidx++]); 
                    point[pidx + i] = 0.0;
                    vv2 += vi*vi;
                }
                //vv2 += gsl_cdf_chisq_Pinv(uniform[uidx++], d-axis);
                point[pidx + axis] = sqrt(vv2);
                vv += vv2;
            }
            if (vv == 0.0) return 1;
            double v0 = gsl_cdf_gamma_Pinv(uniform[uidx++], d-2, 1.0);
            point[pidx + 0] = v0;
            // To keep p_i^2 non-negative we'll make the norm a
            // bit smaller than it should have been.
            double norm = v0*(1.0 - DBL_EPSILON/2.0)/sqrt(vv);
            for (int i = 1; i < d; i++) {
                point[pidx + i] *= norm;
            }
        }
    }
    assert(uidx == pspoint_dim(d, n));
    double *q = alloca(d*sizeof(double));
    // q = \sum_i p_i
    for (int j = 0; j < d; j++) {
        q[j] = point[0 + j];
    }
    for (int pidx = d; pidx < n*d; pidx += d) {
        for (int j = 0; j < d; j++) {
            q[j] += point[pidx + j];
        }
    }
    // qq = q^2
    double qq = q[0]*q[0];
    for (int j = 1; j < d; j++) {
        qq -= q[j]*q[j];
    }
    double im = 1.0/sqrt(qq);
    double x = w*im;
    double a = im/(1.0 + q[0]*im);
    for (int pidx = 0; pidx < n*d; pidx += d) {
        double v0 = point[pidx + 0];
        double bq = 0.0;
        for (int i = 1; i < d; i++) {
            bq += q[i]*point[pidx + i];
        }
        point[pidx + 0] = im*x*(q[0]*v0 - bq);
        for (int i = 1; i < d; i++) {
            point[pidx + i] = x*(point[pidx + i] + im*q[i]*(a*bq - v0));
        }
    }
    if (0) { // Sanity check
        // q = \sum_i p_i
        for (int j = 0; j < d; j++) {
            q[j] = point[0 + j];
        }
        for (int pidx = d; pidx < n*d; pidx += d) {
            for (int j = 0; j < d; j++) {
                q[j] += point[pidx + j];
            }
        }
        // q \approx (w, \bar 0)
        assert(fabs(q[0] - w) < 1e-14);
        for (int i = 1; i < d; i++) {
            assert(fabs(q[i]) < 1e-14);
        }
        for (int pidx = 0; pidx < n*d; pidx += d) {
            // pp = p_i^2
            double pp = point[pidx + 0]*point[pidx + 0];
            for (int i = 1; i < d; i++) {
                pp -= point[pidx + i]*point[pidx + i];
            }
            // pp \approx 0
            assert(fabs(pp) < 1e-12);
        }
    }
    return 0;
}

/* PHASE-SPACE INTEGRALS
 */

#include <cuba.h>
#include <gsl/gsl_monte_vegas.h>

double
sproduct(int d, int n, const double *p, int i, int j)
{
    assert(i < n);
    assert(j < n);
    double sp = 0.0;
    for (int k = 1; k < d; k++) {
        sp += p[i*d + k]*p[j*d + k];
    }
    return p[i*d + 0]*p[j*d + 0] - sp;
}

double
integrand(int integral, int d, int n, const double *x, double alpha)
{
    double *p = alloca(n*d*sizeof(double));
    if (pspoint(d, n, 1.0, x, p) != 0) return 0.0;
#define SP(i, j) sproduct(d, n, p, i-1, j-1)
#define Q1(i) (1.0 - 2.0*p[(i-1)*d + 0])
#define Q(i, j) (1.0 - 2.0*(p[(i-1)*d + 0] + p[(j-1)*d + 0]) + 2.0*SP(i, j))
    switch (integral) {
        case 1: // phase space volume
            return 1.0;
        case 2: // 12 14 23 34
            return (1.0/16)/(SP(1,2)*SP(1,4)*SP(2,3)*SP(3,4) + alpha);
        case 3: // q1 q5
            return (1.0)/(Q1(1)*Q1(5) + alpha);
        case 4: // 12 13 q2 q3
            return (1.0/4)/(SP(1,2)*SP(1,3)*Q1(2)*Q1(3) + alpha);
        case 5: // 12 34 q2 q3
            return (1.0/4)/(SP(1,2)*SP(3,4)*Q1(2)*Q1(3) + alpha);
        case 6: // 12 13 q2 q4 45 35
            return (1.0/16)/(SP(1,2)*SP(1,3)*Q1(2)*Q1(4)*SP(4,5)*SP(3,5) + alpha);
        case 7: // q45 q34
            return (1.0)/(Q(4,5)*Q(3,4) + alpha);
        case 8: // 13 14 q45 q35
            return (1.0/4)/(SP(1,3)*SP(1,4)*Q(4,5)*Q(3,5) + alpha);
        case 9: // 12 23 q23 q12
            return (1.0/4)/(SP(1,2)*SP(2,3)*Q(2,3)*Q(1,2) + alpha);
        case 10: // q1 q45
            return (1.0)/(Q1(1)*Q(4,5) + alpha);
        case 11: // q1 q45 q14 q5
            return (1.0)/(Q1(1)*Q(4,5)*Q(1,4)*Q1(5) + alpha);
        case 12: // 12 q23 q5
            return (1.0/2)/(SP(1,2)*Q(2,3)*Q1(5) + alpha);
        case 13: // 12 13 q35 q5 q2
            return (1.0/4)/(SP(1,2)*SP(1,3)*Q(3,5)*Q1(5)*Q1(2) + alpha);
        case 14: // q2 12 q23 q5
            return (1.0/2)/(Q1(2)*SP(1,2)*Q(2,3)*Q1(5) + alpha);
        case 15: // 12 13 q35 q25 q2 q3
            return (1.0/4)/(SP(1,2)*SP(1,3)*Q(3,5)*Q(2,5)*Q1(2)*Q1(3) + alpha);
        case 16: // 12 34 q2 q3 q35 q25
            return (1.0/4)/(SP(1,2)*SP(3,4)*Q1(2)*Q1(3)*Q(3,5)*Q(2,5) + alpha);
        case 17: // 12 34 q13
            return (1.0/4)/(SP(1,2)*SP(3,4)*Q(1,3) + alpha);
        case 18: // 12 13 25 q35 q5 q4
            return (1.0/8)/(SP(1,2)*SP(1,3)*SP(2,5)*Q(3,5)*Q1(5)*Q1(4) + alpha);
        case 19: // 12 13 25 q35 q13 q5
            return (1.0/8)/(SP(1,2)*SP(1,3)*SP(2,5)*Q(3,5)*Q(1,3)*Q1(5) + alpha);
        case 20: // 12 13 q35 q25 35 25
            return (1.0/16)/(SP(1,2)*SP(1,3)*Q(3,5)*Q(2,5)*SP(3,5)*SP(2,5) + alpha);
        case 21: // 12 13 24 q12 q5
            return (1.0/8)/(SP(1,2)*SP(1,3)*SP(2,4)*Q(1,2)*Q1(5) + alpha);
        case 22: // 12 13 24 q12 45 q5
            return (1.0/16)/(SP(1,2)*SP(1,3)*SP(2,4)*Q(1,2)*SP(4,5)*Q1(5) + alpha);
        case 23: // 12 13 25 45 q25 q12
            return (1.0/16)/(SP(1,2)*SP(1,3)*SP(2,5)*SP(4,5)*Q(2,5)*Q(1,2) + alpha);
        case 24: // 12 13 45 q25 q5 q3
            return (1.0/8)/(SP(1,2)*SP(1,3)*SP(4,5)*Q(2,5)*Q1(5)*Q1(3) + alpha);
        case 25: // 12 13 45 q25 q12 q5
            return (1.0/8)/(SP(1,2)*SP(1,3)*SP(4,5)*Q(2,5)*Q(1,2)*Q1(5) + alpha);
        case 26: // q12 q45
            return (1.0)/(Q(1,2)*Q(4,5) + alpha);
        case 27: // 12 q1 q45 q13
            return (1.0/2)/(SP(1,2)*Q1(1)*Q(4,5)*Q(1,3) + alpha);
        case 28: // q45 q12 q5 q1
            return (1.0)/(Q(4,5)*Q(1,2)*Q1(5)*Q1(1) + alpha);
        case 29: // 12 q45 25 q13 q5 q1
            return (1.0/4)/(SP(1,2)*Q(4,5)*SP(2,5)*Q(1,3)*Q1(5)*Q1(1) + alpha);
        case 30: // 12 13 q34 q25 q5 q3
            return (1.0/4)/(SP(1,2)*SP(1,3)*Q(3,4)*Q(2,5)*Q1(5)*Q1(3) + alpha);
        case 31: // 12 13 q34 q25 q5 q4
            return (1.0/4)/(SP(1,2)*SP(1,3)*Q(3,4)*Q(2,5)*Q1(5)*Q1(4) + alpha);
    }
    return 0.0;
}

typedef struct Integrand {
    int integral;
    int n;
    int d;
    double alpha;
} Integrand;

int
integrand_cuba(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
{
    Integrand *i = (Integrand*)userdata;
    int d = i->d;
    int n = i->n;
    double alpha = i->alpha;
    assert(*ndim == pspoint_dim(d, n));
    assert(*ncomp == 1);
    f[0] = integrand(i->integral, d, n, x, alpha);
    return 0;
}

double
integrand_gsl(double *x, size_t ndim, void *userdata)
{
    Integrand *i = (Integrand*)userdata;
    int d = i->d;
    int n = i->n;
    double alpha = i->alpha;
    assert((int)ndim == pspoint_dim(d, n));
    return integrand(i->integral, d, n, x, alpha);
}

/* MAIN
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define countof(arr) (sizeof(arr)/sizeof((arr)[0]))

double
timestamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

void
print_num(double x, double dx)
{
    char s1[512], s2[512], s3[512];
    snprintf(s1, countof(s1), "%.10f", x-dx);
    snprintf(s2, countof(s2), "%.10f", x);
    snprintf(s3, countof(s3), "%.10f", x+dx);
    char *d1 = strchr(s1, '.');
    char *d3 = strchr(s3, '.');
    if ((d1 - s1) != (d3 - s3)) {
        putchar('(');
        for (char *s = s2; *s; s++) putchar(*s);
        putchar(')');
    } else {
        int printedpar = 0;
        for (size_t i = 0; (i < countof(s2)) && s2[i]; i++) {
            if (!printedpar && (s2[i] != s1[i]) && (s2[i] != s3[i])) {
                printedpar = 1;
                putchar('(');
            }
            putchar(s2[i]);
        }
        if (!printedpar) putchar('(');
        putchar(')');
    }
}

void
usage(const char *progname)
{
    fprintf(stderr,
        "Usage: %s [-d num] [-e num] [-E num] integral\n"
        "Arguments:\n"
        "     integral   the number of the integral to calculate (1-31)\n"
        "Options:\n"
        "    -d num      integrate in this many dimensions (default: 6)\n"
        "    -e num      target this relative accuracy (default: 0.01)\n"
        "    -E num      start with this many evaluations (default: 1000000)\n"
        "Environment variables:\n"
        "    CUBACORES   number of cores to use for integration\n",
        progname);
    exit(1);
}

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

int
main(int argc, char *argv[])
{
    int nparticles = 5;
    int ndims = 6;
    int usecuba = 1;
    double epsilon = 0.01;
    long nevals = 1*1000*1000;
    for (int opt; ((opt = getopt(argc, argv, "d:E:e:hCG")) != -1);) {
        switch (opt) {
            case 'd': ndims = atoi(optarg); break;
            case 'e': epsilon = strtod(optarg, NULL); break;
            case 'E': nevals = atol(optarg); break;
            case 'C': usecuba = 1; break;
            case 'G': usecuba = 0; break;
            case 'h': usage(argv[0]); break;
            default: usage(argv[0]); break;
        }
    }
    if (optind == argc) usage(argv[0]);
    setvbuf(stdout, NULL, _IOLBF, 0);
    gsl_rng_env_setup();
    for (; optind < argc; optind++) {
        int integral = atoi(argv[optind]);
        if (!((1 <= integral) && (integral <= 31))) {
            printf("UNKNOWN INTEGRAL: %s\n", argv[optind]);
            continue;
        }
        int nparams = pspoint_dim(ndims, nparticles);
        printf("INTEGRAL = %d\n", integral);
        printf("       D = %d\n", ndims);
        printf("       N = %d\n", nparticles);
        printf(" NPARAMS = %d\n", nparams);
        printf("  NEVALS = %ld\n", nevals);
        printf("     EPS = %.2g\n", epsilon);
        printf("PSVOLUME = %.10g\n", psvolume(ndims, nparticles, 1.0));
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_monte_vegas_state *mc = gsl_monte_vegas_alloc(nparams);
        assert(mc != NULL);
        Integrand iparams = {
            .alpha=1.0,
            .d=ndims,
            .n=nparticles,
            .integral=integral
        };
        gsl_monte_function fn = { integrand_gsl, nparams, (void*)&iparams };
        assert(gsl_monte_vegas_init(mc) == 0);
        long maxnevals = nevals*128;
        for (int k = 30; k <= 100; k += 10) {
            double alpha = pow(2, -k);
            double t1 = timestamp();
            iparams.alpha = alpha;
            printf(" * * * * * * * * *\n");
            printf("   alpha = 2^(-%d)\n", k);
            if (usecuba) { // Integrate with Cuba
                int neval = -1;
                int fail = -1;
                double result = -1;
                double error = -1;
                double prob = -1;
                for (int iter = 0; iter < 5; iter++) {
                    Vegas(nparams, 1,
                            integrand_cuba, (void*)&iparams, 1, // integrand_t integrand, void *userdata, const int nvec,
                            epsilon, 1e-12, // const cubareal epsrel, const cubareal epsabs,
                            1, 0, // const int flags, const int seed,
                            nevals, 4*nevals, // const int mineval, const int maxeval,
                            nevals/5, nevals/11, 10000, // const int nstart, const int nincrease, const int nbatch,
                            0, NULL, NULL, // const int gridno, const char *statefile, void *spin,
                            &neval, &fail, // int *neval, int *fail,
                            &result, &error, &prob);
                    printf("integral = "); print_num(result, error); printf("\n");
                    printf("   error = %.10f\n", error);
                    printf("    prob = %.5g\n", prob);
                    // We're not testing for prob < 0.1 because
                    // Cuba's quasi-random point generation often
                    // makes the integral converge faster than
                    // the variance estimation would tell us,
                    // thus lowering chi-squared value.
                    //
                    // Strictly speaking chi-squared test in
                    // not applicable here at all, but that is
                    // what Cuba gives you.
                    if (prob > 0.9) {
                        printf("     ~ poor chi2 ~\n");
                        nevals = min(3*nevals + nevals/7, maxnevals);
                        continue;
                    }
                    if (fabs(error/result) > epsilon) {
                        printf("     ~ poor accuracy ~\n");
                        nevals = min(3*nevals + nevals/7, maxnevals);
                        continue;
                    }
                    break;
                }
                printf("totevals = %d\n", neval);
            }
            if (!usecuba) { // Integrate with GSL
                int iterations = 5;
                double result = -1;
                double abserr = -1;
                double *x0 = alloca(nparams*sizeof(double));
                double *x1 = alloca(nparams*sizeof(double));
                for (int i = 0; i < nparams; i++) {
                    x0[i] = 0.0;
                    x1[i] = 1.0;
                }
                gsl_monte_vegas_params mcp;
                gsl_monte_vegas_params_get(mc, &mcp);
                mcp.iterations = iterations;
                mcp.verbose = 1;
                mcp.ostream = stdout;
                gsl_monte_vegas_params_set(mc, &mcp);
                if (mcp.stage == 0) {
                    long evals = nevals/4/iterations;
                    printf("  nevals = %ld*%d (warmup)\n", evals, iterations);
                    gsl_monte_vegas_integrate(&fn, x0, x1, nparams, evals, rng, mc, &result, &abserr);
                }
                for (int iter = 0; iter < 10; iter++) {
                    long evals = nevals/iterations;
                    printf("  nevals = %ld*%d\n", evals, iterations);
                    gsl_monte_vegas_integrate(&fn, x0, x1, nparams, evals, rng, mc, &result, &abserr);
                    printf("integral = "); print_num(result, abserr); printf("\n");
                    printf("   error = %.10f\n", abserr);
                    double chisq = gsl_monte_vegas_chisq(mc);
                    printf("chi2/dof = %.10f\n", chisq);
                    if (fabs(chisq - 1.0) > 0.6) {
                        printf("     ~ poor chi2 ~\n");
                        nevals = min(nevals + nevals/2, maxnevals);
                        continue;
                    }
                    if (fabs(abserr/result) > epsilon) {
                        printf("     ~ poor accuracy ~\n");
                        nevals = min(nevals + nevals, maxnevals);
                        continue;
                    }
                    break;
                }
            }
            double t2 = timestamp();
            printf("    time = %.3gs\n", t2-t1);
        }
        gsl_monte_vegas_free(mc);
    }
    return 0;
}
