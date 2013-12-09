CC = g++
CFLAGS = -Wall

erp: erp.c
	$(CC) $(CFLAGS) erp.c -o erp
   
clean:
	rm -rf *.o *~ erp
 
