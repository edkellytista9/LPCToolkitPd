# LPCToolkitPd
Mark Cartwright's LPC Toolkit for Pure Data.

So far, the main piece of work is lpc_lpc~.c - actually resolves to mbc_lpc~ when compiled.
It is the 8th of April, and I have not tested this yet. The aim is to remove any platform-specific libraries (vDSP for Apple, Gnu Scientific Library for Linux etc) and just have simple elegant C code. Apart from math.h - a standard C library, it should only contain the FFT routines from m_pd.h in terms of external functions.
