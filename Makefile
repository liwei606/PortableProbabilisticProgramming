CC = gcc
CFLAGS = -lm -Wall
src = $(shell find ./ -name "*.c")

inference: ppp.h $(src)
	$(CC) $(CFLAGS) -c $(src) -std=c99
   
clean:
	rm -rf *.o *~ erp
 
