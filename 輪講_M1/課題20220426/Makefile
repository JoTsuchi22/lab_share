compile_and_execution:	compile
	mkdir -p temp
	cd temp && ../integration_2D.x

compile:	integration_2D.c
	gcc -o ./integration_2D.x ./integration_2D.c -O3 -lm -Wall -Wextra