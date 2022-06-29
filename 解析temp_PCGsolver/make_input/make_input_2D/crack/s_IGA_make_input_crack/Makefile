default:
	echo Please read README.txt to use this program.

0:	step3_0
	mkdir -p inputdate_result
	cd ./temp && mv input_loc.txt ../inputdate_result -f
	cd ./temp && mv input.svg ../inputdate_result -f
	cd ./temp && mv J_int.txt ../inputdate_result -f
	rm -f -R temp

1:	step3_1
	mkdir -p inputdate_result
	cd ./temp && mv input_loc.txt ../inputdate_result -f
	cd ./temp && mv input.svg ../inputdate_result -f
	cd ./temp && mv J_int.txt ../inputdate_result -f
	rm -f -R temp

step3_0:	step2_0
	cd ./temp && ../s_IGA_make_input_crack.x ../input_crack.txt 4_edited.txt 5_edited.txt 6_edited.txt 7_edited.txt 8_edited.txt 9_edited.txt > temp.dat

step3_1:	step2_1
	cd ./temp && ../s_IGA_make_input_crack.x ../input_crack.txt 16_edited.txt 17_edited.txt 18_edited.txt 19_edited.txt 20_edited.txt 21_edited.txt 22_edited.txt 23_edited.txt 24_edited.txt 25_edited.txt 26_edited.txt 27_edited.txt 28_edited.txt 29_edited.txt 30_edited.txt 31_edited.txt 32_edited.txt 33_edited.txt 34_edited.txt 35_edited.txt 36_edited.txt 37_edited.txt 38_edited.txt 39_edited.txt > temp.dat

step2_0:	step1
	cd ./temp && ../s_IGA_make_input.x ./4.txt ./5.txt ./6.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./7.txt ./8.txt ./9.txt > ./temp.dat

step2_1:	step1
	cd ./temp && ../s_IGA_make_input.x ./16.txt ./17.txt ./18.txt ./19.txt ./20.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./21.txt ./22.txt ./23.txt ./24.txt ./25.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./26.txt ./27.txt ./28.txt ./29.txt ./30.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./31.txt ./32.txt ./33.txt ./34.txt ./35.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./36.txt ./37.txt ./38.txt ./39.txt > ./temp.dat

step1:	s_IGA_make_input_crack.x s_IGA_make_input.x
	mkdir -p temp
	cd ./temp && ../s_IGA_make_input_crack.x ../input_crack.txt > temp.dat

s_IGA_make_input_crack.x:	s_IGA_make_input_crack.c
	gcc -o ./s_IGA_make_input_crack.x ./s_IGA_make_input_crack.c -O3 -lm -Wall -Wextra

s_IGA_make_input.x:	s_IGA_make_input.c
	gcc -o ./s_IGA_make_input.x ./s_IGA_make_input.c -O3 -lm -Wall -Wextra