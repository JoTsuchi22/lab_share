インプットデータのファイルフォーマットが以前のバージョンと違うので注意してください．
MI3D -> MC3D の順で操作すると問題なく作成されるはずです．
詳しくはmake_inpout/make_input_3D/MI3D, MC3D の textfile/template...を参照してください．

コンパイル方法(-O3)

$ cd (Makefileがあるディレクトリ)
$ make

デバッグ用コンパイル方法

$ cd (Makefileがあるディレクトリ)
$ make compile_debug

実行方法

Makefile 内の run のコマンドライン引数のファイル名や数を適宜変更
$ cd (Makefileがあるディレクトリ)
$ make run

データ確認方法


可視化方法
global_patch.xmf
global_patch_boundary_line.xmf
global_patch_control_point.xmf
local_patch.xmf
local_patch_boundary_line.xmf
local_patch_control_point.xmf
(localはS-IGAの場合のみ)

上記のファイルをparaviewのアプリケーションにドラッグしてXDMF Readerを選択すると可視化できます．
paraview:   https://www.paraview.org/download/