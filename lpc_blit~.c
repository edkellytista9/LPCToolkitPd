/*----------------------------------------------------------------------------------
Filename:		mbc.blit~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		5/15/07
Updated:		10/21/10
PD port:                Edward Kelly
Ported:                 02/01/18
Description:	band-limited impulse train generator external object for Max/MSP.  
				Uses the Sum of Windowed Sinc (SWS) method for blit generation.  See
				http://ccrma.stanford.edu/~stilti/papers/TimStilsonPhDThesis2006.pdf 
-------------------------------------------------------------------------------------*/

#include "m_pd.h"
#include <math.h>

#define DEFAULT_PINC 440
#define DEFAULT_P 64;
#define DEFAULT_ZCS 1024;
#define DEFAULT_BANDLIMIT 1;

static t_class *lpc_blit_class;

typedef struct _pulse 
{
	int	active;
	float phase;
} t_pulse;

////////////////////////// object struct
typedef struct _lpc_blit 
{
	t_object					x_obj;			// the object itself (t_pxobject in MSP)
        t_float						b_pinc;			// == frequency
	double						b_phase;		// global phase
	t_float*					b_sincTable;	        // wavetable
	t_pulse*					b_pulse;		// array of pulses
	double						b_fs;			// sampling rate
	int						b_zcs;			// samples per zero crossing in sinc table
	int						b_P;			// P in sinc calc... number of zero crossings
	long						b_length;		// table length
	int						b_bandlimit;	        // bandlimiting on/off
  t_float f;
} t_lpc_blit;

void lpc_blit_winGen(double *win, long N) {
	//this function generates a blackman harris window (low resolution/high dynamic range)
	int i;
	double a0,a1,a2,a3;
	
	a0 = 0.3635819;
	a1 = 0.4891775;
	a2 = 0.1365995;
	a3 = 0.0106411;
	
	for (i = 0; i < N; i++) {
		win[i] = a0 - a1*cos((2*M_PI*i)/(N-1)) + a2*cos((4*M_PI*i)/(N-1)) - a3*cos((6*M_PI*i)/(N-1));
	}
}

void lpc_blit_sincGen(t_lpc_blit *x) {
	int i;
	double P = (double) x->b_P;
	double zcs = (double) x->b_zcs;
	double M = 2*floor(P/2) + 1;
	long length = x->b_length;
	double win[length];
	t_float* pSincTable = x->b_sincTable;
	double*	pWin = win;
	
	//generate blackman-harris window
	lpc_blit_winGen(win,length);
	
	//since this equation makes a zerophase window, we need to split it into 2 to make it causal again
	//first half
	for (i = (floor(length/2) + 1); i <= length; i++) {
		*pSincTable = sin(M_PI*(i/zcs)*(M/P))/(P*sin(M_PI*(i/zcs)/P)) * (*pWin);
		pSincTable++;
		pWin++;
	}
	//second half
	for (i = 1; i < (floor(length/2) + 1); i++) {
		*pSincTable = sin(M_PI*(i/zcs)*(M/P))/(P*sin(M_PI*(i/zcs)/P)) * (*pWin);
		pSincTable++;
		pWin++;
	}	
}

//t_int *lpc_blit_sigperf(t_int *w)
t_int *lpc_blit_perform(t_int *w)
{
  int i;
  t_lpc_blit  *x  =   (t_lpc_blit *)(w[1]);
  t_sample  *pinc =     (t_sample *)(w[2]);
  t_sample  *out  =     (t_sample *)(w[3]);
  int           n =            (int)(w[4]);

  t_float thresh = x->b_fs;
  double phase = x->b_phase;
  t_float phasen1 = phase;
  long length = x->b_length;
  int zcs = x->b_zcs;
  t_float offset, step, idx, eta;
  int P = x->b_P;
  int PO2 = P / 2;
		
  while (n--) {
    phase += *pinc;
    if (phase >= thresh)
      {
	offset = 1.0 - ((thresh - phasen1)/(phase - phasen1));
	for (i = 0; i < PO2; i++)
	  {
	    if (!(x->b_pulse[i].active))
	      {
		x->b_pulse[i].active = 1;
		x->b_pulse[i].phase = offset;
		break;
	      }
	  }
	phase -= thresh;
      }
    phasen1 = phase;
		
		//render active pulses
    *out = 0.0;
    for (i = 0; i < PO2; i++)
      {
	if (x->b_pulse[i].active)
	  {
	    step = x->b_pulse[i].phase;
	    idx = step * zcs;
	    eta = idx - floor(idx);
	    idx = floor(idx);
	    *out += (float)(((1.0 - eta) * x->b_sincTable[(long)(idx)] + eta * x->b_sincTable[((long)(idx + 1.0)) % length]) * 0.89); //linear interpolation, the 0.89 is a scaling factor to keep peak below 1 (and therefore no aliasing... good up to 20k
	    step += 1.0;
	    if (step > (float)(P-1))
	      {
		x->b_pulse[i].active = 0;
		x->b_pulse[i].phase = 0.0;
	      }
	    else
	      {
		x->b_pulse[i].phase = step;
	      }
	  }
      }
    out++;
    pinc++;
  }	
  x->b_phase = phase;	
  return (w+5);
}

/*t_int *lpc_blit_sigperf_a(t_int *w)
{
  t_lpc_blit *x = (t_lpc_blit *)(w[1]);
  t_float *pinc = (t_float *)(w[2]);
  t_float *out = (t_float *)(w[3]);
  int n = (int)(w[4]);
  float thresh = x->b_fs;
  double phase = x->b_phase;
		
  while (n--)
    {
      phase += *pinc;
      if (phase >= thresh)
	{
	  *out = 1.0;
	  phase -= thresh;
	}
      else
	{
	  *out = 0.0;
	}

      out++;
      pinc++;
    }
  x->b_phase = phase;
  return (w+5);
}*/

/*t_int *lpc_blit_fltperf(t_int *w)
{
  int i;

  t_lpc_blit *x = (t_lpc_blit *)(w[1]);
  t_float *out = (t_float *)(w[2]);
  int n = (int)(w[3]);
  t_float thresh = x->b_fs;
  double phase = x->b_phase;
  t_float phasen1 = phase;
  t_float pinc = x->b_pinc;
  long length = x->b_length;
  int zcs = x->b_zcs;
  t_float offset, step, idx, eta;
  int P = x->b_P;
  int PO2 = P / 2;
		
  while (n--)
    {
      phase += pinc;
      if (phase >= thresh)
	{
	  //find sub-sample offset
	  offset = 1.0 - ((thresh - phasen1)/(phase - phasen1));
	  for (i = 0; i < PO2; i++)
	    {
	      if (!(x->b_pulse[i].active))
		{
		  x->b_pulse[i].active = 1;
		  x->b_pulse[i].phase = offset;
		  break;
		}
	    }
	  phase -= thresh;
	}
      phasen1 = phase;
		
      //render active pulses
      *out = 0.0;
      for (i = 0; i < PO2; i++)
	{
	  if (x->b_pulse[i].active)
	    {
	      step = x->b_pulse[i].phase;
	      idx = step * zcs;
	      eta = idx - floor(idx);
	      idx = floor(idx);
	      *out += (float)(((1.0 - eta) * x->b_sincTable[(long)(idx)] + eta * x->b_sincTable[((long)(idx + 1.0)) % length]) * 0.89); //linear interpolation, the 0.89 is a scaling factor to keep peak below 1 (and therefore no aliasing... good up to 20k
	      step += 1.0;
	      if (step > (float)(P-1))
		{
		  x->b_pulse[i].active = 0;
		  x->b_pulse[i].phase = 0.0;
		}
	      else
		{
		  x->b_pulse[i].phase = step;
		}
	    }
	}
      out++;
    }
  x->b_phase = phase;
  return (w+4);
}*/

/*t_int *lpc_blit_fltperf_a(t_int *w)
{
  t_lpc_blit *x = (t_lpc_blit *)(w[1]);
  t_float *out = (t_float *)(w[2]);
  int n = (int)(w[3]);
  t_float thresh = x->b_fs;
  double phase = x->b_phase;
  t_float pinc = x->b_pinc;
		
  while (n--)
    {
      phase += pinc;
      if (phase >= thresh)
	{
	  *out = 1.0;
	  phase -= thresh;
	}
      else
	{
	  *out = 0.0;
	}

      out++;
    }
  x->b_phase = phase;
	
  return (w+4);
  }*/

void lpc_blit_dsp(t_lpc_blit *x, t_signal **sp)
{
  x->b_fs = sys_getsr();
  x->b_phase = 0;

  /*  if (count[0])
    { // perform signal based frequency update
      if (x->b_bandlimit)
	dsp_add(lpc_blit_sigperf, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
      else
	dsp_add(lpc_blit_sigperf_a, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
    }
  else
    { //perform interrupt based frequency update
      if (x->b_bandlimit)
	dsp_add(lpc_blit_fltperf, 3, x, sp[1]->s_vec, sp[1]->s_n);
      else
	dsp_add(lpc_blit_fltperf_a, 3, x, sp[1]->s_vec, sp[1]->s_n);
	}*/

  dsp_add(lpc_blit_perform, 4, x,
  sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void lpc_blit_init(t_lpc_blit* x)
{
  int i;
  for(i = 0; i < (x->b_P / 2); i++)
    {
      x->b_pulse[i].active = 0;
      x->b_pulse[i].phase = 0.0;
    }
  x->b_phase = 0.0;
}

void lpc_blit_free(t_lpc_blit* x)
{
  //dsp_free((t_pxobject *) x);
  //sysmem_freeptr(x->b_sincTable);
  freebytes(x->b_sincTable, x->b_length * sizeof(t_float));
  freebytes(x->b_pulse, (x->b_P / 2) * sizeof(t_pulse));
}

void lpc_blit_float(t_lpc_blit *x, t_float f)
{
  x->b_pinc = f;
}

void *lpc_blit_new(t_symbol *s, long argc, t_atom *argv)
{
  //	t_blit *x = NULL;
  t_lpc_blit *x = (t_lpc_blit *)pd_new(lpc_blit_class);
  //t_bmt_tilde *x = (t_bmt_tilde *)pd_new(bmt_tilde_class);

  outlet_new(&x->x_obj, gensym("signal"));
		
  x->b_pinc = DEFAULT_PINC;
  x->b_P = DEFAULT_P;
  x->b_zcs = DEFAULT_ZCS;
  x->b_bandlimit = DEFAULT_BANDLIMIT;
  //get arguments out of gimme list
  switch(argc)
    {
    case 0:
      break;			
    case 1:
      x->b_pinc = atom_getfloatarg(0,argc,argv);
      break;		
    case 2:
      x->b_pinc = atom_getfloatarg(0,argc,argv);
      x->b_P = (int)atom_getfloatarg(1,argc,argv);
      break;		
    case 3:
      x->b_pinc = atom_getfloatarg(0,argc,argv);
      x->b_P = (int)atom_getfloatarg(1,argc,argv);
      x->b_zcs = (int)atom_getfloatarg(2,argc,argv);
      break;
    case 4:
      x->b_pinc = atom_getfloatarg(0,argc,argv);
      x->b_P = (int)atom_getfloatarg(1,argc,argv);
      x->b_zcs = (int)atom_getfloatarg(2,argc,argv);
      x->b_bandlimit = (int)atom_getfloatarg(3,argc,argv);
      break;
    default:
      x->b_pinc = atom_getfloatarg(0,argc,argv);
      x->b_P = (int)atom_getfloatarg(1,argc,argv);
      x->b_zcs = (int)atom_getfloatarg(2,argc,argv);
      x->b_bandlimit = (int)atom_getfloatarg(3,argc,argv);
      error("mbc.blit~: too many arguments");
      break;
    }		
  if ((int)(pow(2.0,round(log2((double)x->b_P)))) != x->b_P)
    {
      error("mbc.blit~: P argument must be a power of 2");
      x->b_P = pow(2.0,floor(log2(x->b_P)));
    }
	
  x->b_length = (long)(x->b_P * x->b_zcs - 1);
  x->b_sincTable = (t_float *) getzbytes(x->b_length * sizeof(t_float));
  x->b_pulse = (t_pulse *) getzbytes( (x->b_P / 2) * sizeof(t_pulse));
  lpc_blit_init(x);
  lpc_blit_sincGen(x);
  
  return (x);
}

void lpc_blit_setup(void)
{
  lpc_blit_class = class_new(gensym("lpc_blit~"),
				     (t_newmethod)lpc_blit_new,
			     (t_newmethod)lpc_blit_free, sizeof(t_lpc_blit),
				     0L, A_GIMME, 0);

  //post("|=======bmt~========|");
  //post("|=bass==mid==treble=|");
  //post("|=ed==kelly===2010==|");

  class_addmethod(lpc_blit_class, (t_method)lpc_blit_dsp,
		  gensym("dsp"), 0);
  class_addfloat(lpc_blit_class, lpc_blit_float);
  //  class_addmethod(lpc_blit_class, (t_method)lpc_blit_bandlimit, gensym("bandlimit"), A_DEFFLOAT, 0);
  CLASS_MAINSIGNALIN(lpc_blit_class, t_lpc_blit, f);
}
