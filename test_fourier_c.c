#include <math.h>    /* mathematische Funktionen */
#include <complex.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265359
#endif // M_PI

//Minimalbeispiel in C
//hier ist das Problem, dass Werte gar nicht erst eingesetzt werden (?)

long bitrev(long n,int m)
/* Bit Umkehr von  n mit m dualen Stellen  */
{
  long ntilde; /* fuer das Ergebnis */
  int p;       /* fuer Schleife ueber die dualen Stellen */
  ldiv_t divres;

  ntilde=0;

  for(p=0;p<m;p++)
    {
      divres=ldiv (n,2);
      n=divres.quot;
      ntilde=2*ntilde+divres.rem;
    }
  return ntilde;

}

void sortbitrev(long n,int m,double complex *h)
/* Routine sortiert das Feld h der Laenge n
   in bit-umgekehrter Reihenfolge wobei die Anzahl
   der dualen Stellen m ist. Es sollte n=2^m sein. */
{
  double complex z;
  long i,j;

  for(i=0;i<n;i++)
    {
      j=bitrev(i,m);  /* i soll vertauscht werden mit j */
      if(j<i)  /* verhindert, dass zurueckgetauscht wird */
	{
	  z=h[i];
	  h[i]=h[j];
	  h[j]=z;
	}
    }

}


void fft(long n,int m,double complex *h,double sign)
/* FFT Routine ersetzt den Datensatz in h mit Laenge n
   durch die diskrete Fouriertransformation (n=pow(m,2) !!!)
   sign +1.0 ist FFT, sign=-1.0 ist inverse FFT
*/
{
  long iblock,nblock,k;
  int mtilde;
  double complex z,*w,fe,fo;  /* exp von omega Schritt, even und odd Element in Schritt */
  double fakt;


  /* Teste n,m Relation */

  if(n!=(long)(pow(2.0,(double)m)+0.1))
    {
      printf("Zusammenhang zwischen n und m nicht korrekt: %ld %d \n",n,m);
      abort();
    }

  if(!(sign==1.0 || sign==-1.0))
    {
      printf("sign ist nicht korrekt %ld  %d \n",n,m);
      abort();
    }

  /* exp mit Schrittweite fuer jedes  mtilde */

  w=malloc(sizeof(double complex)*(m+1));  /*reserviere Speicherplatz fuer das w */

  /* Schleife belegt w[i]=exp(2*M_PI*i*2^(m-1-i)/n) */

  w[0]=cexp(2*M_PI*I/((double)n));
  if(sign==-1.0)
   {
     w[0]=conj(w[0]);
   }

  for(mtilde=1;mtilde<=m;mtilde++)
    {
      w[mtilde]=w[mtilde-1]*w[mtilde-1];
    }

  /* Sortierung der Daten an FFT Schema anpassen */

  sortbitrev(n,m,h);

  /* rekursive Transformation durchfuehren */

  nblock=1;
  for(mtilde=m-1;mtilde>=0;mtilde--)
      /* diese Schleife zaehlt die mtilde so durch wie beschrieben */
    {
      for(iblock=0;iblock<n;iblock+=2*nblock)
	/* die k-Abhaengig liegt in den ersten Stellen der Dualzahl d1 ... dmtilde 0000.. + k*
	   das iblock springt die unterschiedlichen d1  ... dmtilde Bloecke an */
	{
	  z=1.+0.*I;
	  for(k=0;k<nblock;k++)
	    /* schliesslich der k-Loop, der hier durch die
               k = 2^(m-1-mtilde)/2 laeuft  (also die fuer die Elemente d1 ... dmtilde+1 */
	    {
	      fe=h[k+iblock];           // Element d1 ... dmtilde 0  + k
	      fo=h[k+iblock+nblock];    // Element d1 ... dmtilde 1  + k


	      h[k+iblock]=fe+fo*z;                // bestimme d1 ... dmtilde + k
	      h[k+iblock+nblock]=fe+fo*z*w[m-1];  // bestimme d1 ... dmtilde + k + 2^(m-1-mtilde)
	      z=z*w[mtilde];     // damit ist oben z        = W^(2^(mtilde) * k)
                                 //   und          z*w[m-1] = W^(2^(mtilde) * k + 2^(m-1))
	    }	                 //                         = W^(2^(mtilde) * (k+ 2^(m-1-mtilde)))
     	}

      nblock=nblock*2;           // nur damit nblock = 2^(m-1-mtilde) ist
    }

  if(sign==-1.0)
    {
      fakt=1.0/((double)n);
      for(k=0;k<n;k++)
	{
	  h[k]=h[k]*fakt;
	}
    }

  free(w);

}

int main() {

    int m = 13;
    long n = pow(2,m);
    double tstart=-10, tend=10;
    double complex *f_t, *f_omega;

    double tstep=(tend-tstart)/(n-1);

    f_t = malloc(n*sizeof(double complex));
    f_omega =  malloc(n*sizeof(double complex));

    for(int i=0; i<n; i++) {
        double t=tstart+i*tstep;
        f_t[i]=0;
        if (abs(t)<0.5) {       //Rechtecksfunktion
            f_t[i]=1;
        f_omega[i] = f_t[i];    //Kopie, damit f_omega transformiert werden kann während f_t so bleibt
        }
    }

    fft(n, m, f_omega, 1);  //es sollte rauskommen: N*sinc(w/2)

    int n_print = 50;       //Anzahl der Ausgaben
    int n_skip = n / n_print;   //wieviele Schritte in der Loop übersprungen werden müssen

    for(int i=0; i<n; i+=n_skip) {
        printf("t= %4.2f, f(t)= %4.2f, f(omega)=%4.2f\n", tstart+i*tstep, f_t[i], f_omega[i]);
    }

    free(f_t);
    free(f_omega);

    return 0;

}
