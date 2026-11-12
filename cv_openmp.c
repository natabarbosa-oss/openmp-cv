/*
 * Exercício Final: Análise de Dados com OpenMP
 * Coeficiente de variação (estatura e peso) para dados agrupados em
 * intervalos de classe.
 *   - Estatura: classes de 8 cm
 *   - Peso:     classes de 4 kg
 *
 * Compilação:
 *   Versão paralela:    gcc -O2 -fopenmp cv_openmp.c -o cv_par -lm
 *   Versão sequencial:  gcc -O2          cv_openmp.c -o cv_seq -lm
 *   (sem -fopenmp as diretivas #pragma são ignoradas -> código sequencial)
 *
 * Uso:
 *   ./cv_par dados.txt [fator_replicacao] [repeticoes]
 *     fator_replicacao : replica o conjunto de dados N vezes (só para
 *                        benchmark; não altera o CV). Padrão = 1.
 *     repeticoes       : quantas vezes o cálculo é repetido para medir
 *                        o tempo médio. Padrão = 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define H_ALT  8.0   /* amplitude de classe da estatura (cm) */
#define H_PESO 4.0   /* amplitude de classe do peso (kg)     */

typedef struct {
    double lo;        /* limite inferior da 1ª classe            */
    double h;         /* amplitude da classe                     */
    int    k;         /* número de classes                       */
    long  *freq;      /* frequências fi                          */
    double media;     /* média ponderada (mu)                    */
    double desvio;    /* desvio padrão populacional (sigma)      */
    double cv;        /* coeficiente de variação (%)             */
} Dist;

#ifdef _WIN32
#include <windows.h>
static double agora(void) {
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}
#else
static double agora(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
#endif

static long ler_dados(const char *arq, double **alt, double **peso) {
    FILE *f = fopen(arq, "r");
    if (!f) { perror("fopen"); exit(1); }
    long cap = 1024, n = 0;
    double *a = malloc(cap * sizeof(double));
    double *p = malloc(cap * sizeof(double));
    double x, y;
    while (fscanf(f, "%lf %lf", &x, &y) == 2) {
        if (n == cap) {
            cap *= 2;
            a = realloc(a, cap * sizeof(double));
            p = realloc(p, cap * sizeof(double));
        }
        a[n] = x; p[n] = y; n++;
    }
    fclose(f);
    *alt = a; *peso = p;
    return n;
}

/* Calcula as estatísticas agrupadas das duas variáveis. */
static void calcular(const double *alt, const double *peso, long N,
                     Dist *da, Dist *dp)
{
    double amin = alt[0],  amax = alt[0];
    double pmin = peso[0], pmax = peso[0];

    /* (1) Mínimos e máximos: redução paralela */
    #pragma omp parallel for reduction(min:amin,pmin) reduction(max:amax,pmax)
    for (long i = 0; i < N; i++) {
        if (alt[i]  < amin) amin = alt[i];
        if (alt[i]  > amax) amax = alt[i];
        if (peso[i] < pmin) pmin = peso[i];
        if (peso[i] > pmax) pmax = peso[i];
    }

    da->h = H_ALT;  da->lo = floor(amin);
    dp->h = H_PESO; dp->lo = floor(pmin);
    da->k = (int)floor((amax - da->lo) / da->h) + 1;
    dp->k = (int)floor((pmax - dp->lo) / dp->h) + 1;

    long *fa = calloc(da->k, sizeof(long));
    long *fp = calloc(dp->k, sizeof(long));
    const int ka = da->k, kp = dp->k;
    const double loa = da->lo, ha = da->h, lop = dp->lo, hp = dp->h;

    /* (2) Distribuição de frequências: cada thread acumula num histograma
     *     privado (reduction sobre array) -> sem condição de corrida.   */
    #pragma omp parallel for reduction(+:fa[:ka], fp[:kp])
    for (long i = 0; i < N; i++) {
        fa[(int)((alt[i]  - loa) / ha)]++;
        fp[(int)((peso[i] - lop) / hp)]++;
    }
    da->freq = fa; dp->freq = fp;

    /* (3) Média ponderada: soma(xi*fi) / N */
    double sa = 0.0, sp = 0.0;
    #pragma omp parallel for reduction(+:sa)
    for (int i = 0; i < ka; i++) sa += (loa + (i + 0.5) * ha) * fa[i];
    #pragma omp parallel for reduction(+:sp)
    for (int i = 0; i < kp; i++) sp += (lop + (i + 0.5) * hp) * fp[i];
    da->media = sa / N;
    dp->media = sp / N;

    /* (4) Desvio padrão: soma(fi*(xi-mu)^2) / N */
    double va = 0.0, vp = 0.0;
    const double ma = da->media, mp = dp->media;
    #pragma omp parallel for reduction(+:va)
    for (int i = 0; i < ka; i++) {
        double d = (loa + (i + 0.5) * ha) - ma;
        va += fa[i] * d * d;
    }
    #pragma omp parallel for reduction(+:vp)
    for (int i = 0; i < kp; i++) {
        double d = (lop + (i + 0.5) * hp) - mp;
        vp += fp[i] * d * d;
    }
    da->desvio = sqrt(va / N);
    dp->desvio = sqrt(vp / N);

    /* (5) Coeficiente de variação */
    da->cv = da->desvio / da->media * 100.0;
    dp->cv = dp->desvio / dp->media * 100.0;
}

static void imprimir(const char *nome, const char *un, const Dist *d, long N) {
    printf("\n=== %s ===\n", nome);
    printf("%-18s %-10s %-8s\n", "Classe", "Ponto med.", "fi");
    for (int i = 0; i < d->k; i++) {
        double li = d->lo + i * d->h, ls = li + d->h;
        printf("%6.1f |-- %6.1f %10.1f %8ld\n", li, ls, li + d->h / 2, d->freq[i]);
    }
    printf("N = %ld\n", N);
    printf("Media (mu)      = %.4f %s\n", d->media, un);
    printf("Desvio (sigma)  = %.4f %s\n", d->desvio, un);
    printf("CV              = %.4f %%\n", d->cv);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Uso: %s dados.txt [fator_replicacao] [repeticoes]\n", argv[0]);
        return 1;
    }
    long fator = (argc > 2) ? atol(argv[2]) : 1;
    int  reps  = (argc > 3) ? atoi(argv[3]) : 1;

    double *alt0, *peso0;
    long n0 = ler_dados(argv[1], &alt0, &peso0);

    /* replicação opcional (apenas para benchmark) */
    long N = n0 * fator;
    double *alt = alt0, *peso = peso0;
    if (fator > 1) {
        alt  = malloc(N * sizeof(double));
        peso = malloc(N * sizeof(double));
        for (long j = 0; j < fator; j++)
            for (long i = 0; i < n0; i++) { alt[j*n0+i] = alt0[i]; peso[j*n0+i] = peso0[i]; }
    }

#ifdef _OPENMP
    printf("Versao PARALELA (OpenMP) - threads = %d\n", omp_get_max_threads());
#else
    printf("Versao SEQUENCIAL\n");
#endif

    Dist da, dp;
    double t0 = agora();
    for (int r = 0; r < reps; r++) {
        if (r > 0) { free(da.freq); free(dp.freq); }
        calcular(alt, peso, N, &da, &dp);
    }
    double t = (agora() - t0) / reps;

    imprimir("ESTATURA (classes de 8 cm)", "cm", &da, N);
    imprimir("PESO (classes de 4 kg)", "kg", &dp, N);
    printf("\nTempo medio por execucao: %.9f s (N=%ld, reps=%d)\n", t, N, reps);
    printf("TEMPO %.9f\n", t);
    return 0;
}
