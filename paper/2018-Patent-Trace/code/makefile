DEPS:=generator.c strand.c util.c -g -fsanitize=address -Wall

all: main

main:
	gcc main.c $(DEPS)

demo:
	gcc test/test_main.c $(DEPS)

clean:
	rm a.out
