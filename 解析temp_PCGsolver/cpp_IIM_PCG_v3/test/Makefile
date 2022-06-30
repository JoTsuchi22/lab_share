CC				= gcc
CFLAGS			= -O1 -lm

compile:
	$(CC) -o a.x s_iga_3D_analysis_GaussPt_select.c $(CFLAGS)

run:
	cd ./output && ../a.x ../input_glo.txt ../input_loc.txt

testrun:	compile	run