CC				= gcc
CFLAGS			= -O1 -lm -Wall -Wextra

compile:
	$(CC) -o a.x s_IGA_for_debug_10x10_20211208.c $(CFLAGS)

run:
	cd ./output && ../a.x ../input.txt

view:
	../src/NURBS_viewer/NURBS_view.exe ./output/input_for_NURBS.txt ./output/view.dat ./output/Displacement.dat

testrun:	compile	run	view
	