
ノットインサーションやオーダーエレベーションを行い，
inputdataを編集するためのプログラムです．

2次元2パラメータ空間のNURBSのみ実装しています．

編集可能なパッチは1つだけです．

入力データのテキストファイルはtemplate_input.txtを参照してください．

コンパイル
gcc -g -o s_IGA_make_input.x s_IGA_make_input.c -Og -lm -Wall
gcc -g -o s_IGA_make_input.x s_IGA_make_input.c -O3 -lm -Wall

実行
./s_IGA_make_input.x test_input_01.txt > result.dat
./s_IGA_make_input.x test_input_02.txt > result.dat
./s_IGA_make_input.x test_input_03.txt > result.dat