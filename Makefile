#Makefile for EM_PoissonMM
#by Adrianna (4/23/2017)
#_____________________________

#Define Macros:
OBJS   = EMBasins.o BasinModel.o
CC     = g++
DEBUG  = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

EM_PoissonMM : $(OBJS)
	$(CC) $(LFLAGS) -lgsl $(OBJS) -o EM_PoissonMM

EMBasins.o : EMBasins.cpp EMBasins.h BasinModel.h 
	$(CC) -O3 -c -fPIC EMBasins.cpp  

BasinModel.o : BasinModel.cpp BasinModel.h EMBasins.h
	$(CC) -O3 -c -fPIC BasinModel.cpp

clean:
	\rm $(OBJS) 

tar: 
	tar cfv EM_PoissonMM.tar EMBasins.cpp EMBasins.h BasinModel.cpp BasinModel.h
