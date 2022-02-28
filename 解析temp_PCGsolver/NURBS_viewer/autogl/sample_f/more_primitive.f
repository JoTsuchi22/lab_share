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

c     透視投影モードにする。
      call autogl_set_perspective_view_flag(1)

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

      real*8 radius
      integer n_divisions
      
      radius = 10.0d0
      n_divisions = 4

c     鏡面反射をオンにする。
      call autogl_turn_on_specular

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, 'O')
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     50.0d0, 0.0d0, 0.0d0)
      call autogl_draw_string(50.0d0, 0.0d0, 0.0d0, 'X')
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 50.0d0, 0.0d0)
      call autogl_draw_string(0.0d0, 50.0d0, 0.0d0, 'Y')
      call autogl_draw_line(
     c     0.0d0, 0.0d0, 0.0d0,
     c     0.0d0, 0.0d0, 50.0d0)
      call autogl_draw_string(0.0d0, 0.0d0, 50.0d0, 'Z')

c     球を描画する。
      call autogl_set_color(1.0d0, 1.0d0, 0.0d0) 
      call autogl_draw_sphere3d(
     c     10.0d0, 10.0d0, 0.0d0, 
     c     radius, n_divisions)

c     円柱を描画する。
      call autogl_set_color(0.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_cylinder3d(
     c     -10.0d0, -20.0d0, 30.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     radius, n_divisions)

c     円錐を描画する。
      call autogl_set_color(1.0d0, 0.0d0, 1.0d0) 
      call autogl_draw_cone3d(
     c     -10.0d0, -20.0d0, -30.0d0, 
     c     10.0d0, 20.0d0, -30.0d0, 
     c     radius, n_divisions)

c     円を描画する。
      call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      call autogl_draw_circle3d(
     c     -10.0d0, -10.0d0, 0.0d0, 
     c     1.0d0, 1.0d0, 0.0d0, 
     c     radius, n_divisions)

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
