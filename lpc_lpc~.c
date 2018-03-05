/*
 *
 *
 *
 *
 *
 *
 */

#include "m_pd.h"
#include <math.h>
//#include <gsl/gsl_math.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_blas.h>

static t_class *mbc_lpc_class;

#define MAX_ORDER 200
#define NLOG2(x) (ceil(log2(x)))
#define POW2(x) (1 << x);
#define DEFAULT_FS 44100
#define DEFAULT_FRAMERATE 100
#define DEFAULT_V_SIZE 64

typedef struct _DSPcomplex
{
  gsl_vector_float* real;
  gsl_vector_float* imag;
} t_DSPcomplex;

#define REAL(z,i) gsl_vector_float_set(z,2*(i))
#define IMAG(z,i) gsl_vector_float_setz(z,2*(i)+1)

typedef struct _lpc 
{
  t_object          x_obj;	           // the object itself (t_pxobject in MSP)
  //t_float*          l_frame_buff;	   //input frame buffer
  gsl_vector_float*          l_frame_buff;	   //input frame buffer
  //t_float*          l_winframe_buff;       //windowed input frame buffer
  gsl_vector_float*          l_winframe_buff;       //windowed input frame buffer
  //t_float*	    l_outCoeff_buff;       //coefficient signal
  gsl_vector_float*	    l_outCoeff_buff;       //coefficient signal
  //t_float*	    l_outParcor_buff; 	   //PARCOR coeffs
  gsl_vector_float*	    l_outParcor_buff; 	   //PARCOR coeffs
  //t_float*          l_outError_buff;	   //error signal
  gsl_vector_float*          l_outError_buff;	   //error signal
  //t_float*          l_win;	           //analysis window
  gsl_vector_float*          l_win;	           //analysis window
  //t_float*	    l_R;
  gsl_vector_float*	    l_R;
  double*	    l_W;
  double*	    l_E;
  double*	    l_K;
  double 	    l_G;
  double**	    l_A;
  double 	    l_x1;                  //last samples of pre-emphasis filter
  float 	    l_b1;	           //pre-emph coefficient
  int 	            l_order;	           //predictor order
  int 	            l_order_max;           //max order according to fs = order * frame_rate
  int 	            l_preemph;             //pre-epmhasis filter on/off
  int 	            l_frame_rate;          //analysis frame rate
  int 	            l_frame_size;          //analysis frame size, where fs = frame_rate * frame_size * 2
  int 	            l_hop_size;            //hop_size = frame_size * 2 (b/c of overlap)
  int 	            l_inframe_idx;         //current inframe buffer index
  int 	            l_outframe_idx;	   //current outframe buffer index
  long 	            l_v_size;	           //vector size
  float 	    l_fs;                  //sampling rate
  int 	            l_maxnfft;             //fft length
  int 	            l_log2nfft;	           //log2(fft length)
  //  FFTSetup          l_fftsetup;            //FFTSetup for vDSP FFT functions
  t_DSPcomplex      split;
  //  DSPSplitComplex   l_fftSplitTemp;        //temp split complex vector structure
} t_lpc;

t_int *lpc_perform(t_int *w) 
{
  t_lpc   *x          = (t_lpc *)  (w[1]);
  t_float *in         = (t_float *)(w[2]);
  t_float *out_error  = (t_float *)(w[3]);
  t_float *out_gain   = (t_float *)(w[4]);
  t_float *out_coeff  = (t_float *)(w[5]);
  t_float *out_parcor = (t_float *)(w[6]);
  t_float *out_index  = (t_float *)(w[7]);
  int n               = (int)      (w[8]);
  int p = x->l_order;
  int length = x->l_frame_size;
  int log2nfft = NLOG2(length+p+1);
  int nfft = POW2(log2nfft);
  int nfftO2 = nfft/2;
  t_float recipn = 1.0/nfft;
  int inframeidx = x->l_inframe_idx;
  int outframeidx = x->l_outframe_idx;
	
  int i, j, i1, ji;
  int in_idx = 0, out_idx = 0;
  double val;
  t_float scale, rCorr;
	
  if (x->l_preemph)
    {
      while (n--) 
	{
	  val = in[in_idx];
	  in[in_idx] = val + x->l_b1 * x->l_x1;
	  x->l_x1 = val;
	  in_idx++;
	}
      n = (int)(w[8]);
      in_idx = 0;
    }
  
  while (n--) 
    {
      if (inframeidx < length) {
	//copy input into frame buff
	//also copy into winframe buff since GSL puts the result into one of the arg vectors
	gsl_vector_float_set(&x->l_frame_buff,inframeidx,in[in_idx]);// = x->l_frame_buff[inframeidx] = in[in_idx];
	
	out_gain[out_idx] = x->l_G;
	out_error[out_idx] = gsl_vector_float_get(&x->l_outError_buff,outframeidx); //for now
	if (outframeidx < x->l_order) {
	  out_coeff[out_idx] = gsl_vector_float_get(&x->l_outCoeff_buff,outframeidx);
	  out_parcor[out_idx] = gsl_vector_float_get(&x->l_outParcor_buff,outframeidx);
	  out_index[out_idx] = outframeidx + 1;
	} else {
	  out_coeff[out_idx] = 0.0;
	  out_parcor[out_idx] = 0.0;
	  out_index[out_idx] = 0;
	}
	
	inframeidx++;
	in_idx++;
	outframeidx++;
	out_idx++;
      } else {
	gsl_vector_float_memcpy(&x->l_winframe_buff,&x->l_frame_buff);
	 
	//perform durbin-levinson - for right now, just count to the order---------------
	//clear memory, is this necessary?
	for (i=0; i < p+1; i++){
	  x->l_R[i] = 0.0;
	  x->l_W[i] = 0.0;
	  x->l_E[i] = 0.0;
	  x->l_K[i] = 0.0;			
	}
	for(i=0; i<=p; i++) {
	  for(j=0; j < p; j++) x->l_A[i][j] = 0.0;
	}
	//window frame buff
	//vDSP_vmul(x->l_frame_buff, 1, x->l_win, 1, x->l_winframe_buff, 1, length);
	gsl_vector_float_mul(&x->l_winframe_buff,&x->l_win);
#ifdef DEBUG
	for(i=0;i<length;i++)post("\nwinframe(%d) = %g;",i+1,gsl_vector_float_get(&x->l_winframe_buff,i));
#endif
	
	//create r from auto correlation
	if ((2*nfft*log2(nfft)+nfft) > length*p) { //NOTE: change this to update only when order is changed!
	  //time domain method
	  //for(i=0; i < p+1; i++) vDSP_dotpr(x->l_winframe_buff,1,x->l_winframe_buff+i,1,x->l_R+i,length-i);
	  //gsl_vector_float_memcpy(x->l_R,x->l_winframe_buff);
	  for(i=0; i < p+1; i++)
	    {
	      gsl_blas_sdot(&x->l_winframe_buff,&x->l_winframe_buff+i,&rCorr);
	      gsl_vector_float_set(&x->l_R,i,rCorr);
	    }
	} else {
	  //frequency domain method
	  // zero pad
	  //vDSP_vclr(x->l_winframe_buff+length,1,nfft-length);
	  for(i=length;i<nfft;i++)
	    {
	      gsl_vector_float_set(&x->l_winframe_buff,i,0);
	    }
	  //for(i=0;i<nfft02;i++)
	  //  {
	  //    //convert to split complex vector
	  //    gsl_vector_float_set(&x->split.real,i,gsl_vector_float_get(x->l_winframe_buff,i*2));
	  //    gsl_vector_float_set(&x->split.imag,i,gsl_vector_float_get(x->l_winframe_buff,i*2+1));
	  //  }
	  //vDSP_ctoz( ( DSPComplex * ) x->l_winframe_buff, 2, &x->l_fftSplitTemp, 1, nfftO2);
	  //perform forward in place fft
	  gsl_fft_complex_radix2_forward(&x->l_winframe_buff,1,nfft);

	  //this is as far as I (Ed Kelly) have got at 14:39pm on March 5 2018!
	  
	  vDSP_fft_zrip(x->l_fftsetup, &x->l_fftSplitTemp, 1, log2nfft, FFT_FORWARD);
	  //scaling
	  scale = 0.5;
	  vDSP_vsmul(x->l_fftSplitTemp.realp,1,&scale,x->l_fftSplitTemp.realp,1,nfftO2);
	  vDSP_vsmul(x->l_fftSplitTemp.imagp,1,&scale,x->l_fftSplitTemp.imagp,1,nfftO2);
	  //compute PSD
	  vDSP_zvmags(&x->l_fftSplitTemp,1,x->l_fftSplitTemp.realp,1,nfftO2);
	  //clear imaginary part
	  vDSP_vclr(x->l_fftSplitTemp.imagp,1,nfftO2);
	  //perform inverse in place fft
	  vDSP_fft_zrip(x->l_fftsetup, &x->l_fftSplitTemp, 1, log2nfft, FFT_INVERSE);
	  //scaling
	  vDSP_vsmul(x->l_fftSplitTemp.realp,1,&recipn,x->l_fftSplitTemp.realp,1,nfftO2);
	  vDSP_vsmul(x->l_fftSplitTemp.imagp,1,&recipn,x->l_fftSplitTemp.imagp,1,nfftO2);
	  //convert back to real number vector
	  vDSP_ztoc(&x->l_fftSplitTemp, 1, (DSPComplex *)x->l_R, 2, nfftO2);	
	}
	
	
	/*for(i=0; i < p+1; i++) {
	  x->l_R[i] = 0.0;
	  for(j=0; j < length - i; j++) {
	  x->l_R[i] += x->l_winframe_buff[j] * x->l_winframe_buff[i+j];
	  //x->l_R[i] += x->l_win[j] * x->l_win[i+j];
	  }
	  }*/
#ifdef DEBUG
	for(i=0;i< p+1; i++) post("\nR(%d) = %f;",i+1,x->l_R[i]);
#endif
	
	x->l_W[0] = x->l_R[1];
	x->l_E[0] = x->l_R[0];
	
	for (i = 1; i <= p; i++) {
	  x->l_K[i] = x->l_W[i-1] / x->l_E[i-1];
	  
	  x->l_A[i][i] = x->l_K[i];
	  i1 = i - 1;
	  if (i1 >= 1) {
	    for (j = 1; j <=i1; j++) {
	      ji = i - j;
	      x->l_A[j][i] = x->l_A[j][i1] - x->l_K[i] * x->l_A[ji][i1];
	    }				
	  }
	  
	  x->l_E[i] = x->l_E[i-1] * (1.0 - x->l_K[i] * x->l_K[i]);
	  
	  if (i != p) {
	    x->l_W[i] = x->l_R[i+1];
	    for (j = 1; j <= i; j++) {
	      x->l_W[i] -= x->l_A[j][i] * x->l_R[i-j+1];
	    }
	  }
	}
	
	x->l_G = sqrt(x->l_E[p]);
	for (i=0; i < p; i++) {
	  x->l_outCoeff_buff[i] = (float)(x->l_A[i+1][p]);
	  x->l_outParcor_buff[i] = (float)(x->l_K[i+1]);
#ifdef DEBUG
	  //post("\nParcor(%d) = %g;",i+1,x->l_K[i+1]);
	  //post("\nCoeff(%d) = %g;",i+1,x->l_A[i+1][p]);
#endif
	}
	
	//--------------------------------------------------------------------------------
	
	//copy right side to left side, move delayed input to output
	for (i=0; i < x->l_hop_size; i++) {
	  x->l_outError_buff[i] = x->l_frame_buff[i];
	  x->l_frame_buff[i] = x->l_frame_buff[i + x->l_hop_size];
	}
	
	inframeidx = x->l_hop_size;
	outframeidx = 0;
	n++; //to avoid skipping a sample (since we already decremented
	while (n--) {
	  x->l_frame_buff[inframeidx] = in[in_idx];
				
	  out_gain[out_idx] = (float)(x->l_G);
	  out_error[out_idx] = x->l_outError_buff[outframeidx]; //for now
	  if (outframeidx < x->l_order) {
	    out_coeff[out_idx] = x->l_outCoeff_buff[outframeidx];
	    out_parcor[out_idx] = x->l_outParcor_buff[outframeidx];
	    out_index[out_idx] = (float)(outframeidx + 1);
	  } else {
	    out_coeff[out_idx] = 0.0;
	    out_parcor[out_idx] = 0.0;
	    out_index[out_idx] = 0;
	  }
	  
	  inframeidx++;
	  in_idx++;
	  outframeidx++;
	  out_idx++;
	}
	break;
      }
    }
  
  x->l_inframe_idx = inframeidx;
  x->l_outframe_idx = outframeidx;
  
  return (w+9);
}





void lpc_hanning(t_lpc *x)
{
  int N = x->l_frame_size;
  int n;
	
  for (n=0; n<N; n++)
    {
      x->l_win[n] = 0.5 * (1.0 - cos((TWOPI*(n+1))/(N+1)));
    }
}

void lpc_hamming(t_lpc *x)
{
  int N = x->l_frame_size;
  int n;
	
  for (n=0; n<N; n++)
    {
      x->l_win[n] = (0.54 - 0.46*cos((TWOPI*n)/(N-1)));
    }
}

void lpc_bartlett(t_lpc *x)
{
  int N = x->l_frame_size;
  int n;

  for (n=0; n<N; n++)
    {
      x->l_win[n] = 2.0/(N - 1) * ((N-1)/2 - fabs(n - (N-1)/2.0));
    }
}

void lpc_preemph(t_lpc *x, int n)
{
  x->l_preemph = n;
}

void lpc_order(t_lpc *x, int n)
{
  if (x->l_order_max < n)
    {
      error("mbc.lpc~: framerate * order must be less than or equal to sampling rate.  At the current setting maxorder is %d.  For a higher order, lower the framerate.", x->l_order_max);
      x->l_order = x->l_order_max;
    }
  else
    {
      x->l_order = n;
    }
    //lpc_init(x) ?
}

void *lpc_new(t_symbol *s, int argc, t_atom *argv)
{
  t_lpc *x = (t_lpc *)pd_new(lpc_class);

  x->l_frame_buff = gsl_vector_float_calloc(MAX_ORDER);	   //input frame buffer
  x->l_winframe_buff = gsl_vector_float_calloc(MAX_ORDER);       //windowed input frame buffer
  x->l_outCoeff_buff = gsl_vector_float_calloc(MAX_ORDER);       //coefficient signal
  x->l_outParcor_buff = gsl_vector_float_calloc(MAX_ORDER); 	   //PARCOR coeffs
  x->l_outError_buff = gsl_vector_float_calloc(MAX_ORDER);	   //error signal
  x->l_win = gsl_vector_float_calloc(MAX_ORDER);	           //analysis window
  



}
