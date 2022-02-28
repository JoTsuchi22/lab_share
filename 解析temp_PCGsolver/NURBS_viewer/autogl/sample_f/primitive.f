c     
c     メインルーチン
c     
      implicit none

c     ビューイングパラメータ（サイズ,視線方向など）を設定する。
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

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 

c     原点ラベルを表示する。
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, 'O')

c     X軸を描画する。
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     50.0d0, 0.0d0, 0.0d0)
      call autogl_draw_string(50.0d0, 0.0d0, 0.0d0, 'X')

c     Y軸を描画する。
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 50.0d0, 0.0d0)
      call autogl_draw_string(0.0d0, 50.0d0, 0.0d0, 'Y')

c     Z軸を描画する。
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 0.0d0, 50.0d0)
      call autogl_draw_string(0.0d0, 0.0d0, 50.0d0, 'Z')

c     三角形を描画する。 

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
c     三点を指定して三角形 
      call autogl_draw_triangle(
     c     -20.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, -20.0d0, 0.0d0, 
     c     0.0d0, 0.0d0, -20.0d0)

c     直方体(0, 0, 0) -（10, 20, 30)の6枚の面（四辺形）を描画する。

      call autogl_set_color(0.0d0, 1.0d0, 1.0d0) 
c     四点を指定して四辺形 
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0)
      call autogl_set_color(1.0d0, 1.0d0, 0.0d0) 
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)
      call autogl_set_color(1.0d0, 0.0d0, 1.0d0) 
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)
      call autogl_set_color(0.0d0, 0.0d0, 1.0d0) 
      call autogl_draw_quadrangle(
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0)
      call autogl_set_color(0.0d0, 1.0d0, 0.0d0) 
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)
      call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      call autogl_draw_quadrangle(
     c     0.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
