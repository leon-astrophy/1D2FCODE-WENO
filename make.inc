F90=pgf90
F90FLAGS= -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium
MODINC="-I ./"
LDFLAGS= -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium
	