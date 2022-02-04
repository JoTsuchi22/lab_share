# NURBSviewer
Zarusoba Viewerを参考に作られています。
AutoGL ライブラリ（OpenGL と GTK+）に依存しています。

## ビルド方法
	% cd /path/to/NURBSviewer
	% chmod +x ./tools/modify_autogl_lib.sh
	% make

## プログラムの説明と使用方法
### NURBS_input:可視化用点群データ作成  
使用法:  

		./NURBS_input input.txt Displacement.dat view.dat div_xi div_eta

入出力:  

		input.txt          IGAの入力ファイル
		Displacement.dat   IGAの変位ファイル
		view.dat           NURBS_viewに読み込ませる為の出力ファイル
		div_xi             ξ方向の一要素あたりの分割数(int)
		div_eta            η方向の一要素あたりの分割数(int)

### NURBS_view:可視化ソフト  
使用法(NURBS_inputでDisplacement.dat未使用):  

		./NURBS_view input.txt view.dat

入力:  

		input.txt          IGAの入力ファイル
		view.dat           NURBS_inputで作成した出力ファイル

使用法(NURBS_inputでDisplacement.dat使用):  

		./NURBS_view input.txt view.dat Displacement.dat

入力:  

		input.txt          IGAの入力ファイル
		view.dat           NURBS_inputで作成した出力ファイル
		Displacement.dat   IGAの変位ファイル

