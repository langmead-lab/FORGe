#!gmake

CFLAGS = -O0 -g -Wall -Wextra -Werror -pedantic

bounter: bounter.o murmur3.o
	$(CC) -o $@ $^ $(CFLAGS)

bounter.o: bounter.c cms_log8.c cms_common.c murmur3.h
	$(CC) -c $< $(CFLAGS) -DBOUNTER_TEST_MAIN

murmur3.o: murmur3.c murmur3.h
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f *.o bounter
