c     
c     メインルーチン
c     
      implicit none

c     ビュー上に表示されるモデルを表現するための変数群 

c     原点のラベル : 文字列 
      character*80 origin_label

c     座標軸を表示するかどうか : 0または1の整数（論理値） 
      integer axis_is_shown

c     モデルの形のタイプ : 0 - 2 の整数 (SHAPE_) 
      integer shape

c     整数変数 
      integer number

c     実数変数 
      real*8 value

      common value, origin_label, axis_is_shown, shape, number

c     ビューおよびパネル１つの例 



c     モデル用変数を初期化する。 
      origin_label = 'hello, world !'
      axis_is_shown = 1
      shape = 2
      number = 20
      value = 0.5d0



c     ビューを準備する。 

c     ビューイングパラメータ（サイズ,視線方向など）を設定する。
c     ビューの大きさ 
      call autogl_set_view_size(70.0d0)
c     注視点 
      call autogl_set_view_center(10.0d0, 10.0d0, 10.0d0)
c     視線方向 
      call autogl_set_view_direction(3.0d0, 1.0d0, 2.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode3d()
      call autogl_enable_drag_in_mode3d()

c     パネルを準備する。 
c     ボタン用サブルーチン番号やコントロール用変数を登録する。
      call autogl_set_panel_in_mode3d()

c     コメントの登録 : ラベル 
      call autogl_add_comment()
      call autogl_set_label('Model Parameters')

c     文字列変数の登録 : テキスト入力 
c     フィールド origin_label 
c     ???
c      call autogl_add_string(origin_label, 'origin_label', 80)
c     ???
c     文字列、変数名、文字列の長さ 

c     整数変数の登録 : テキスト入力 
c     フィールド number 
      call autogl_add_integer(number, 'number')
c     integer型変数と変数名 

c     実数変数の登録 : テキスト入力 
c     フィールド value 
      call autogl_add_real(value, 'value')
c     real*8型変数と変数名 

c     論理変数の登録 : トグル入力 
c     トグル show axis 
      call autogl_add_boolean(axis_is_shown, 'axis_is_shown')
c     integer型変数と変数名 
c     axis_is_shownは0または1の値を取る。 
      call autogl_set_label('show axis')
c     ラベルを指定してもよい。 
c     デフォルトでは変数名がラベルとして使用される。 

c     整数変数の登録 : ラジオボタン入力 
c     チョイス shape 
      call autogl_add_integer(shape, 'shape')
c     integer型変数と変数名
c     shapeは0から2までの値を取る。 
c     値0のアイテム 
      call autogl_add_integer_item('line')
c     値1のアイテム 
      call autogl_add_integer_item('triangle')
c     値2のアイテム 
      call autogl_add_integer_item('box')     

c     コールバックサブルーチンの登録 : ボタン 
c     描画アップデートボタン Draw 
      call autogl_set_button(0, 'Draw')
c     コールバックサブルーチンの番号とボタン名
c     サブルーチンdraw_buttonを指定
      
c     コールバックサブルーチンの登録 : ボタン 
c     印刷ボタン Print 
      call autogl_set_button(1, 'Print')
c     コールバックサブルーチンの番号とボタン名
c     サブルーチンprint_buttonを指定

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      real*8 value
      character*80 origin_label
      integer axis_is_shown
      integer shape
      integer number
      common value, origin_label, axis_is_shown, shape, number

      character*80 buf

      if (axis_is_shown .eq. 1) then

c     X軸を描画する。
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 0.0d0,
     c        50.0d0, 0.0d0, 0.0d0)
         
c     Y軸を描画する。
         call autogl_set_color(0.0d0, 1.0d0, 0.0d0) 
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 0.0d0,
     c        0.0d0, 50.0d0, 0.0d0)
         
c     Z軸を描画する。
         call autogl_set_color(0.0d0, 0.0d0, 1.0d0) 
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 0.0d0,
     c        0.0d0, 0.0d0, 50.0d0)
      endif

c     原点ラベルを表示する。 
      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, origin_label)

c     形のタイプに応じて, 
      if (shape .eq. 0) then
c     SHAPE_LINE
c     直方体の１２本の稜線を線分として描画する。 
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 0.0d0,
     c        10.0d0, 0.0d0, 0.0d0)
         call autogl_draw_line(
     c        10.0d0, 0.0d0, 0.0d0,
     c        10.0d0, 20.0d0, 0.0d0)
         call autogl_draw_line(
     c        10.0d0, 20.0d0, 0.0d0,
     c        0.0d0, 20.0d0, 0.0d0)
         call autogl_draw_line(
     c        0.0d0, 20.0d0, 0.0d0,
     c        0.0d0, 0.0d0, 0.0d0)
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 30.0d0,
     c        10.0d0, 0.0d0, 30.0d0)
         call autogl_draw_line(
     c        10.0d0, 0.0d0, 30.0d0,
     c        10.0d0, 20.0d0, 30.0d0)
         call autogl_draw_line(
     c        10.0d0, 20.0d0, 30.0d0,
     c        0.0d0, 20.0d0, 30.0d0)
         call autogl_draw_line(
     c        0.0d0, 20.0d0, 30.0d0,
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_line(
     c        0.0d0, 0.0d0, 0.0d0,
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_line(
     c        0.0d0, 20.0d0, 0.0d0,
     c        0.0d0, 20.0d0, 30.0d0)
         call autogl_draw_line(
     c        10.0d0, 20.0d0, 0.0d0,
     c        10.0d0, 20.0d0, 30.0d0)
         call autogl_draw_line(
     c        10.0d0, 0.0d0, 0.0d0,
     c        10.0d0, 0.0d0, 30.0d0)

      elseif (shape .eq. 1) then
c     SHAPE_TRIANGLE
c     原点を通る３枚の三角形を描画する。
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
         call autogl_draw_triangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_triangle(
     c        10.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_triangle(
     c        10.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 0.0d0, 0.0d0)

      elseif (shape .eq. 2) then
c     SHAPE_BOX
c     直方体(0, 0, 0) -（10, 20, 30)の6枚の面（四辺形）を描画する。
         call autogl_set_color(0.0d0, 1.0d0, 1.0d0) 
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 30.0d0, 
     c        10.0d0, 0.0d0, 30.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 20.0d0, 30.0d0)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_quadrangle(
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        10.0d0, 0.0d0, 30.0d0)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 30.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         call autogl_draw_quadrangle(
     c        0.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 20.0d0, 30.0d0)

      endif

c     位置(10, 20, 30)に整数変数numberの値を表示する。 
      write(buf,*) 'number:', number
      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_string(10.0d0, 20.0d0, 30.0d0, buf)

c     画面の左下隅から比率で(0.1, 0.1)の場所に実数変数valueの値を表示する。 
      write(buf,*) 'value = ', value
      call autogl_set_color(0.0d0, 0.0d0, 1.0d0) 
      call autogl_draw_string_at_scaled_position(0.1d0, 0.1d0, buf)

      end



c     描画アップデートコマンドのためのコールバックサブルーチン
c     "Draw"ボタンを選ぶと呼ばれる。 
      subroutine draw_button
      implicit none

c     ビューを再描画する。 
      call autogl_draw_view()

      end



c     印刷コマンドのためのコールバックサブルーチン
c     "Print"ボタンを選ぶと呼ばれる。 
      subroutine print_button
      implicit none

      real*8 value
      character*80 origin_label
      integer axis_is_shown
      integer shape
      integer number
      common value, origin_label, axis_is_shown, shape, number

      write(*,*) 'PRINT'

c     現在の原点ラベルを表示する。 
      write(*,*) ' origin_label is ', origin_label

c     現在座標軸が表示されているかどうかを表示する。 
      write(*,*) ' axis_is_shown is ', axis_is_shown

c     現在のモデルの形のタイプを表示する。 
      if (shape .eq. 0) then
         write(*,*) ' shape is SHAPE_LINE (0)'
      elseif (shape .eq. 1) then
         write(*,*) ' shape is SHAPE_TRIANGLE (1)'
      elseif (shape .eq. 2) then
         write(*,*) ' shape is SHAPE_BOX (2)'
      endif

c     整数変数numberの値を表示する。 
      write(*,*) ' number is ', number

c     実数変数valueの値を表示する。 
      write(*,*) ' value is ', value

      end



c     コールバックサブルーチン
      subroutine autogl_callback(id)
      implicit none

      integer id

      if (id .eq. 0) then
         call draw_button()
      elseif (id .eq. 1) then
         call print_button()
      endif

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

