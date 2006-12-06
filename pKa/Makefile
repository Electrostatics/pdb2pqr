CPP = g++ -O3
FLAGS = -c
.SUFFIXES: .cpp.o
.cpp.o:
	$(CPP) $(FLAGS) $<

all:	clean pKa_MC_mult

clean:	
	rm -f $(ALLOBJS) pMC_mult_wrap* pMC_mult.o _pMC_mult.so

PKA_MCMULTOBJS = pMC_mult.o	

ALLOBJS = 	$(PKA_MCMULTOBJS)

#
# To find the correct linker option
# python
# >>> import distutils.sysconfig
#>>> distutils.sysconfig.get_config_var('LINKFORSHARED')
#'-Xlinker -export-dynamic'
#

pKa_MC_mult:	$(PKA_MCMULTOBJS)
		swig -python -c++ pMC_mult.i
		$(CPP) -c pMC_mult.cpp pMC_mult_wrap.cxx -I/usr/include/python2.3 -I/usr/include/python2.4
		$(CPP) -Xlinker -export-dynamic -shared pMC_mult.o pMC_mult_wrap.o  -o _pMC_mult.so 

depend:
	makedepend -I/usr/include/c++/3.2/ -- -E --  pMC_mult.cpp

# DO NOT DELETE THIS LINE -- make depend depends on it.

