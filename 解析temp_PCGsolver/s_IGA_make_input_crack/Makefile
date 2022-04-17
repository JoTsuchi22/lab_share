default:
	echo Please read README.txt to use this program.

0:	step3_0

1:	step3_1
	rm -f -R temp

step3_0:	step2_0
	./s_IGA_make_input_crack.x ./input_crack_test.txt ./tmep/4_edited.txt ./tmep/5_edited.txt ./tmep/6_edited.txt ./tmep/7_edited.txt ./tmep/8_edited.txt ./tmep/9_edited.txt

step3_1:	step2_1
	./s_IGA_make_input_crack.x ./input_crack_test.txt ./tmep/16_edited.txt ./tmep/17_edited.txt ./tmep/18_edited.txt ./tmep/19_edited.txt ./tmep/20_edited.txt ./tmep/21_edited.txt ./tmep/22_edited.txt ./tmep/23_edited.txt ./tmep/24_edited.txt ./tmep/25_edited.txt ./tmep/26_edited.txt ./tmep/27_edited.txt ./tmep/28_edited.txt ./tmep/29_edited.txt ./tmep/30_edited.txt ./tmep/31_edited.txt ./tmep/32_edited.txt ./tmep/33_edited.txt ./tmep/34_edited.txt ./tmep/35_edited.txt ./tmep/36_edited.txt ./tmep/37_edited.txt ./tmep/38_edited.txt ./tmep/39_edited.txt > ./temp/temp.dat

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
	cd ./temp && ../s_IGA_make_input_crack.x ../input_crack_test.txt > ./temp.dat

s_IGA_make_input_crack.x:	s_IGA_make_input_crack.c
	gcc -o ./s_IGA_make_input_crack.x ./s_IGA_make_input_crack.c -O3 -lm -Wall -Wextra

s_IGA_make_input.x:	s_IGA_make_input.c
	gcc -o ./s_IGA_make_input.x ./s_IGA_make_input.c -O3 -lm -Wall -Wextra