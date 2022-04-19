
基底関数2次の特異パッチを用いたs_IGA解析用インプットデータを作成します．

インプットデータ(input_crack.txt)についてはinput_crack_template.txtを参照してください．

s_IGA_make_input.c と s_IGA_make_input_crack.c，input_crack.txt，Makefile が
同じディレクトリにあることを確認してmakeを利用して実行してください．

1/4 モデルの場合
$   make 0

full モデルの場合
$   make 1

make Error 1 が出た場合，インプットデータに問題がある可能性があるので，
./temp/temp.dat の最終行を参照してください