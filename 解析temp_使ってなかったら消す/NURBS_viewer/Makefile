.PHONY: all autogl lib input graph view tools clean-lib clean-input clean-graph clean-tools clean-view

all: autogl lib input graph view tools 

autogl:
	tar zxvf autogl.tar.gz && cd autogl && ../tools/modify_autogl_lib.sh lib && make all

lib:
	cd lib && make all

input:
	cd NURBS_input && make all

graph:
	cd NURBS_graph && make all

tools:
	cd tools && make all

view:
	cd NURBS_view && make all

clean: clean-lib clean-input clean-view clean-tools

clean-lib:
	cd lib && make clean

clean-input:
	cd NURBS_input && make clean

clean-graph:
	cd NURBS_graph && make clean

clean-tools:
	cd tools && make clean

clean-view:
	cd NURBS_view && make clean
