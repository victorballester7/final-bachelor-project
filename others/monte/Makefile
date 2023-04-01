include Makedefs

EXECS=test_monte

############
# Utilitats
############
# - Test Montenbruck vs nosaltres
test_monte : test_monte.cpp SAT_Force.o SAT_VecMat.o SAT_RefSys.o
	g++ -o test_monte $(OPT) $(CFLAGS) test_monte.cpp SAT_Force.o SAT_VecMat.o SAT_RefSys.o -lm

##########
# Rutines
##########
# - Camp gravitatori expandit en harmònics esfèrics

SAT_Force.o : SAT_Force.cpp
	g++ -c $(OPT) $(CFLAGS) SAT_Force.cpp

SAT_VecMat.o : SAT_VecMat.cpp
	g++ -c $(OPT) $(CFLAGS) SAT_VecMat.cpp

SAT_RefSys.o : SAT_RefSys.cpp
	g++ -c $(OPT) $(CFLAGS) SAT_RefSys.cpp

#########
# Neteja
#########

clean :
	rm -f *.o

realclean : clean
	rm -f $(EXECS)
