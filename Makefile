EXEC=it_simu
DFLAG=-Wall
CFLAG=-std=c99 -lm -llapacke -lblas -march=native -Ofast

all:
	gcc main.c fonctions.c -o $(EXEC) $(CFLAG) $(DFLAG) -fopenmp

dgemm:
	gcc main.c fonctions.c $(CFLAG) $(DFLAG) -D DGEMM -o $(EXEC)

debug:
	gcc main.c fonctions.c $(CFLAG) $(DFLAG) -fopenmp -D DEBUG -o $(EXEC)

run:
	./$(EXEC)

help:
	@echo "****** Makefile help ******"
	@echo "Here is a summary of all the options"
	@echo "make\t\t This will compile all the source files."
	@echo "make debug\t This will compile all the source files with a debug flag."
	@echo "make dgemm\t This will compile all the source files with dgemm function instead of our own matrix product."
	@echo "make run\t This will execute the code."
	@echo "make clean\t This will delete the execution file."

clean:
	rm $(EXEC) 
