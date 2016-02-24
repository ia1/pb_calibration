all: get_thr_ploidy get_mixture get_thr_ploidy_vcfa get_mixture_vcfa

get_thr_ploidy_vcfa:get_thr_ploidy_vcfa.o
	gcc -lm get_thr_ploidy_vcfa.o -g -o get_thr_ploidy_vcfa -lm
	
get_thr_ploidy_vcfa.o : get_thr_ploidy_vcfa.c
	gcc -c -g get_thr_ploidy_vcfa.c -o get_thr_ploidy_vcfa.o

get_mixture_vcfa:get_mixture_vcfa.o
	gcc -lm get_mixture_vcfa.o -g -o get_mixture_vcfa -lm

get_mixture_vcfa.o : get_mixture_vcfa.c
	gcc -c -g get_mixture_vcfa.c -o get_mixture_vcfa.o

get_mixture:get_mixture.o
	gcc -lm get_mixture.o -g -o get_mixture -lm

get_mixture.o : get_mixture.c
	gcc -c -g get_mixture.c -o get_mixture.o

get_thr_ploidy:get_thr_ploidy.o
	gcc -lm get_thr_ploidy.o -g -o get_thr_ploidy -lm
	
get_thr_ploidy.o : get_thr_ploidy.c
	gcc -c -g get_thr_ploidy.c -o get_thr_ploidy.o

