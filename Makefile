CUBA_CFLAGS=-I.
CUBA_LIBS=libcuba.a

GSL_CFLAGS=
GSL_LIBS=-lgsl -lgslcblas -lm

XCFLAGS=${CFLAGS} ${CUBA_CFLAGS} ${GSL_CFLAGS} -std=c99 -O2 -g -Wall -Wextra -pedantic
XLDFLAGS=${LDFLAGS} ${CUBA_LIBS} ${GSL_LIBS} -lm

all: rambo

rambo: rambo.c
	${CC} ${XCFLAGS} -o $@ rambo.c ${XLDFLAGS}

clean:
	rm -f rambo

RESULTS=\
	results/int2-dim6.txt results/int2-dim8.txt \
	results/int3-dim4.txt results/int3-dim6.txt results/int3-dim8.txt \
	results/int4-dim6.txt results/int4-dim8.txt \
	results/int5-dim6.txt results/int5-dim8.txt \
	results/int6-dim6.txt results/int6-dim8.txt \
	results/int7-dim4.txt results/int7-dim6.txt results/int7-dim8.txt \
	results/int8-dim6.txt results/int8-dim8.txt \
	results/int9-dim6.txt results/int9-dim8.txt \
	results/int10-dim4.txt results/int10-dim6.txt results/int10-dim8.txt \
	results/int11-dim4.txt results/int11-dim6.txt results/int11-dim8.txt \
	results/int12-dim6.txt results/int12-dim8.txt \
	results/int13-dim6.txt results/int13-dim8.txt \
	results/int14-dim6.txt results/int14-dim8.txt \
	results/int15-dim6.txt results/int15-dim8.txt \
	results/int16-dim6.txt results/int16-dim8.txt \
	results/int17-dim6.txt results/int17-dim8.txt \
	results/int18-dim6.txt results/int18-dim8.txt \
	results/int19-dim6.txt results/int19-dim8.txt \
	results/int20-dim6.txt results/int20-dim8.txt \
	results/int21-dim6.txt results/int21-dim8.txt \
	results/int22-dim6.txt results/int22-dim8.txt \
	results/int23-dim6.txt results/int23-dim8.txt \
	results/int24-dim6.txt results/int24-dim8.txt \
	results/int25-dim6.txt results/int25-dim8.txt \
	results/int26-dim4.txt results/int26-dim6.txt results/int26-dim8.txt \
	results/int27-dim6.txt results/int27-dim8.txt \
	results/int28-dim4.txt results/int28-dim6.txt results/int28-dim8.txt \
	results/int29-dim6.txt results/int29-dim8.txt \
	results/int30-dim6.txt results/int30-dim8.txt \
	results/int31-dim6.txt results/int31-dim8.txt

all-results: $(RESULTS)

results/%.txt:
	@mkdir -p results
	./rambo -E 10000000 -e 0.001 $$(echo $* | sed 's/int\([0-9]*\)-dim\([0-9]*\)/-d \2 \1/') | tee $@

