CPP = g++ 
FLAGS = 

all:	clean pKa_MC

clean:	
	rm -f $(ALLOBJS) pMC_wrap* *.o *.so

PKA_MCOBJS =	pMC.o	

ALLOBJS = 	$(PKA_MCOBJS) pKaMC

#
# To find the correct linker option
# python
# >>> import distutils.sysconfig
#>>> distutils.sysconfig.get_config_var('LINKFORSHARED')
#'-Xlinker -export-dynamic'
#


pKa_MC:	$(PKA_MCOBJS) 
	swig -python -c++ pMC.i
	$(CPP) -c pMC.cpp pMC_wrap.cxx -I/usr/local/include/python2.4 -I/software/Python-2.3.4/include/python2.3
	$(CPP) -Xlinker -export-dynamic -shared pMC.o pMC_wrap.o  -o _pMC.so 

depend:
	makedepend -I/usr/include/c++/3.2/ -- -E --  pMC.cpp

# DO NOT DELETE THIS LINE -- make depend depends on it.

pMCo: pMC.h
