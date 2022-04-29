compile_and_execution:	compile
	mkdir -p temp
	cd temp && ../integration_1D.x

compile:	integration_1D.c
	gcc -o ./integration_1D.x ./integration_1D.c -O3 -lm -Wall -Wextra