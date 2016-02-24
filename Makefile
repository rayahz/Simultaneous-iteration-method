EXEC=it_simu

all:
	gcc main.c fonctions.c -std=c99 -lm -llapacke -lblas -Wall -o $(EXEC)

run:
	./$(EXEC)

clean:
	rm $(EXEC) 
