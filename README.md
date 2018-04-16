# LPCToolkitPd
Mark Cartwright's LPC Toolkit for Pure Data.

So far, the main piece of work is mbc_lpc~.c

I have tried to implement this library using only C routines, and I am still uncertain whether I have correctly translated them from the original vDSP Accelerate library functions for Apple Mac OS X.
mbc_lpc~ and mbc_lpc_blit~ compile, but I have not got them to work yet and in the case of mbc_lpc~ I haven't even got it to load properly in Pd.
