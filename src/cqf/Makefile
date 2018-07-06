#!/usr/bin/env make -f

CFLAGS += -std=gnu99 -DSQUEAKR_TEST_MAIN -O0 -g

cqf: cqf.o gqf.o
	$(CC) $^ -lm -o $@

cqf.o: cqf.c cqf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

gqf.o: gqf.c gqf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o cqf
