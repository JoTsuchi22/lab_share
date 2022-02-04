all : compile_lib compile_sample

compile_lib : 
	( cd lib; make )

compile_sample : 
	( cd sample_c; make )
	( cd more_sample_c; make )

clean :
	( cd lib; make clean )
	( cd sample_c; make clean )
	( cd more_sample_c; make clean )
