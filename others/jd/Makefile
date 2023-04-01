# Tot

all: jd2c c2jd

# Fitxers objecte i executables
jd2c: jd2c.o jd.o
	gcc -o jd2c jd2c.o jd.o -lm
jd2c.o: jd2c.c
	gcc -c -pedantic -g -Wall jd2c.c
c2jd: c2jd.o jd.o
	gcc -o c2jd c2jd.o jd.o -lm
c2jd.o: c2jd.c
	gcc -c -pedantic -g -Wall c2jd.c
jd.o: jd.c
	gcc -c -pedantic -g -Wall jd.c

# Neteja
clean:
	rm -f c2jd.o jd2c.o jd.o
realclean: clean
	rm -f jd2c c2jd
