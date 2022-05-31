
基底関数2次の特異パッチを用いたIGA解析用インプットデータを作成します．

インプットデータ(input_crack.txt)についてはinput_crack_template.txtを参照してください．

s_IGA_make_input.c と s_IGA_make_input_crack.c，input_crack.txt，Makefile が
同じディレクトリにあることを確認してmakeを利用して実行してください．
$   make

パッチの作り方を少し変更したモードの実行
Makefie の line 32のコメントアウトを外す