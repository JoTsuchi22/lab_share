CC			= gcc
CFLAGS		= -O3 -lm -Wall -Wextra

run:	step3_1
	mkdir -p inputdate_result
	cd ./temp && mv input_loc.txt ../inputdate_result -f
	cd ./temp && mv input.svg ../inputdate_result -f
	cd ./temp && mv J_int.txt ../inputdate_result -f
	rm -f -R temp

step3_1:	step2_1
	cd ./temp && ../IGA_make_input_crack.x ../input_crack.txt 16_edited.txt 17_edited.txt 18_edited.txt 19_edited.txt 20_edited.txt 21_edited.txt 22_edited.txt 23_edited.txt 24_edited.txt 25_edited.txt 26_edited.txt 27_edited.txt 28_edited.txt 29_edited.txt 30_edited.txt 31_edited.txt 32_edited.txt 33_edited.txt 34_edited.txt 35_edited.txt 36_edited.txt 37_edited.txt 38_edited.txt 39_edited.txt 40_edited.txt 41_edited.txt 42_edited.txt 43_edited.txt 44_edited.txt 45_edited.txt 46_edited.txt 47_edited.txt 48_edited.txt 49_edited.txt 50_edited.txt 51_edited.txt 52_edited.txt 53_edited.txt 54_edited.txt 55_edited.txt 56_edited.txt 57_edited.txt 58_edited.txt 59_edited.txt 60_edited.txt 61_edited.txt 62_edited.txt 63_edited.txt > temp.dat

step2_1:	step1
	cd ./temp && ../s_IGA_make_input.x ./16.txt ./17.txt ./18.txt ./19.txt ./20.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./21.txt ./22.txt ./23.txt ./24.txt ./25.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./26.txt ./27.txt ./28.txt ./29.txt ./30.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./31.txt ./32.txt ./33.txt ./34.txt ./35.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./36.txt ./37.txt ./38.txt ./39.txt ./40.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./41.txt ./42.txt ./43.txt ./44.txt ./45.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./46.txt ./47.txt ./48.txt ./49.txt ./50.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./51.txt ./52.txt ./53.txt ./54.txt ./55.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./56.txt ./57.txt ./58.txt ./59.txt ./60.txt > ./temp.dat
	cd ./temp && ../s_IGA_make_input.x ./61.txt ./62.txt ./63.txt> ./temp.dat

step1:	IGA_make_input_crack.x s_IGA_make_input.x
	mkdir -p temp
	cd ./temp && ../IGA_make_input_crack.x ../input_crack.txt > temp.dat

IGA_make_input_crack.x:	IGA_make_input_crack.c
	$(CC) -o ./IGA_make_input_crack.x ./IGA_make_input_crack_v1_backup.c $(CFLAGS)
	# $(CC) -o ./IGA_make_input_crack.x ./IGA_make_input_crack.c $(CFLAGS)
	# ↑パッチの形状をやや変更したversion

s_IGA_make_input.x:	s_IGA_make_input.c
	$(CC) -o ./s_IGA_make_input.x ./s_IGA_make_input.c $(CFLAGS)