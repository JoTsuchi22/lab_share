01: compaile_and_execution

02:	compaile

03:	execution

compaile_and_execution: compaile
	./s_IGA_make_input_crack.x ./input_crack_test.txt

compaile:	s_IGA_make_input_crack.c
	gcc -o ./s_IGA_make_input_crack.x ./s_IGA_make_input_crack.c -O3 -lm -Wall -Wextra

execution:	s_IGA_make_input_crack.x
	./s_IGA_make_input_crack.x ./input_crack_test.txt