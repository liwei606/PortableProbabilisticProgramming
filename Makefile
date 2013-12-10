CC = gcc
CFLAGS = -lm -Wall -std=c99
src = $(shell find ./ -name "*.c")

test:
	$(CC) erp.c inference.c test.c -o test.exe -Wall -lm -std=c99

inference: ppp.h $(src)
	$(CC) $(CFLAGS) -c $(src) -std=c99
   
clean:
	rm -rf *.o *~ erp
 
