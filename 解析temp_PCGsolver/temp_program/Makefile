CC				= g++
CFLAGS			= -O3 -lm -Wall -Wextra

run:	a.cpp
	$(CC) -o a.x a.cpp $(CFLAGS)
	./a.x

run_txt:	a.cpp
	$(CC) -o a.x a.cpp $(CFLAGS)
	./a.x text.txt

run_obj:	test_obj.cpp
	$(CC) -o test_obj.x test_obj.cpp $(CFLAGS)
	./test_obj.x

donut:	donut.cpp
	$(CC) -o donut.x donut.cpp $(CFLAGS)
	./donut.x