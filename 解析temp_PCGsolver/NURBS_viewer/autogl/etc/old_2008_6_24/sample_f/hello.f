c     
c     メインルーチン
c     
      implicit none

c     サイズや視線方向などを設定する。

c     ビューの大きさ 
      call autogl_set_view_size(70.0d0)

c     注視点 
      call autogl_set_view_center(10.0d0, 10.0d0, 10.0d0)

c     視線方向 
      call autogl_set_view_direction(1.0d0, 2.0d0, 3.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode3d()
      call autogl_set_panel_in_mode3d()
      call autogl_enable_drag_in_mode3d()

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン */
c     ビューが再表示されるごとに呼ばれる。 */
      subroutine autogl_redraw
      implicit none

c     座標軸を白で描画する。

c     白 (r=1, g=1, b=1)
      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 

c     原点に原点ラベルの文字列 
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, 'O')

c     原点から(50, 0, 0)への線分
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     50.0d0, 0.0d0, 0.0d0)

c     X軸 
      call autogl_draw_string(50.0d0, 0.0d0, 0.0d0, 'X')

c     原点から(0, 50, 0)への線分 
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 50.0d0, 0.0d0)

c     Y軸 
      call autogl_draw_string(0.0d0, 50.0d0, 0.0d0, 'Y')

c     原点から(0, 0, 50)への線分 
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 0.0d0, 50.0d0)

c     Z軸 
      call autogl_draw_string(0.0d0, 0.0d0, 50.0d0, 'Z')

c     三角形を赤で描画する。 

c     赤 (r=1, g=0, b=0)
      call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 

c     三点を指定して三角形 
      call autogl_draw_triangle(
     c     -20.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, -20.0d0, 0.0d0, 
     c     0.0d0, 0.0d0, -20.0d0)

c     直方体(0, 0, 0) -（10, 20, 30)をシアンで描画する。 

c     シアン (r=0, g=1, b=1)
      call autogl_set_color(0.0d0, 1.0d0, 1.0d0) 

c     四点を指定して四辺形 
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0)
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)
      call autogl_draw_quadrangle(
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0)
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)
      call autogl_draw_quadrangle(
     c     0.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)

c     いわゆるお決まりのセリフ 

c     黄色 (r=1, g=1, b=0)
      call autogl_set_color(1.0d0, 1.0d0, 0.0d0) 

c     文字列Hello, world ! の表示
      call autogl_draw_string(
     c     10.0d0, 20.0d0, 30.0d0, 
     c     'Hello, world !')

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
