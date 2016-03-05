EXEC=it_simu

all:
	gcc main.c fonctions.c -std=c99 -lm -llapacke -lblas -Wall -o $(EXEC) -march=native -Ofast 

debug:
	gcc main.c fonctions.c -std=c99 -lm -llapacke -lblas -Wall -D DEBUG -o $(EXEC)

run:
	./$(EXEC)

help:
	@echo "****** Makefile help ******"
	@echo "Here is a summary of all the options"
	@echo "make\t\t This will compile all the source files."
	@echo "make debug\t This will compile all the source files with a debug flag."
	@echo "make run\t This will execute the code."
	@echo "make clean\t This will delete the execution file."

clean:
	rm $(EXEC) 
