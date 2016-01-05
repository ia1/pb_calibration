all: get_mixture_hap get_thrRA get_mixture_dip get_thrRA8a get_thrRA8 get_mixture_auto


get_mixture_auto:get_mixture_auto.o
	gcc -lm get_mixture_auto.o -g -o get_mixture_auto -lm
	
get_mixture_hap:get_mixture_hap.o
	gcc -lm get_mixture_hap.o -g -o get_mixture_hap -lm

get_mixture_hap.o : get_mixture_hap.c
	gcc -c -g get_mixture_hap.c -o get_mixture_hap.o


get_thrRA:get_thrRA.o
	gcc -lm get_thrRA.o -g -o get_thrRA -lm

get_thrRA.o : get_thrRA.c
	gcc -c -g get_thrRA.c -o get_thrRA.o
	
get_thrRA8:get_thrRA8.o
	gcc -lm get_thrRA8.o -g -o get_thrRA8 -lm

get_thrRA8.o : get_thrRA8.c
	gcc -c -g get_thrRA8.c -o get_thrRA8.o	
	
get_thrRA8a:get_thrRA8a.o
	gcc -lm get_thrRA8a.o -g -o get_thrRA8a -lm
	
get_thrRA8a.o : get_thrRA8a.c
	gcc -c -g get_thrRA8a.c -o get_thrRA8a.o	



get_mixture_dip:get_mixture_dip.o
	gcc -lm get_mixture_dip.o -g -o get_mixture_dip -lm

get_mixture_dip.o : get_mixture_dip.c
	gcc -c -g get_mixture_dip.c -o get_mixture_dip.o
	