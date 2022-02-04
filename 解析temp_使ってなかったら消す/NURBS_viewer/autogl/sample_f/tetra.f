c     
c     メインルーチン
c     
      implicit none

      integer i_draw_axis, i_color
      real*8 x_size, y_size, z_size
      common i_draw_axis, i_color, x_size, y_size, z_size

      write(*,*) 'set up'
      i_draw_axis = 1
      i_color = 0
      x_size = 20.0
      y_size = 30.0
      z_size = 40.0
      call autogl_set_view_size(70.0d0)

      call autogl_set_default_callback_in_mode3d()
      call autogl_set_panel_in_mode3d()
      call autogl_enable_drag_in_mode3d()

      call autogl_add_boolean(i_draw_axis, 'i_draw_axis')
      call autogl_set_label('draw axis')

      call autogl_add_integer(i_color, 'i_color')
      call autogl_set_label('color')
      call autogl_add_integer_item('red')
      call autogl_add_integer_item('yellow')
      call autogl_add_integer_item('blue')

      call autogl_add_real(x_size, 'x_size')
      call autogl_add_real(y_size, 'y_size')
      call autogl_add_real(z_size, 'z_size')

      call autogl_set_button(0, 'draw')
      call autogl_set_button(1, 'print')

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      integer i_draw_axis, i_color
      real*8 x_size, y_size, z_size
      common i_draw_axis, i_color, x_size, y_size, z_size

      if (i_draw_axis .eq. 1) then
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
         call autogl_draw_line(0.0d0, 0.0d0, 0.0d0, 
     c        50.0d0, 0.0d0, 0.0d0)
         call autogl_draw_string(50.0d0, 0.0d0, 0.0d0, 'X')
         call autogl_draw_line(0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 50.0d0, 0.0d0)
         call autogl_draw_string(0.0d0, 50.0d0, 0.0d0, 'Y')
         call autogl_draw_line(0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 0.0d0, 50.0d0)
         call autogl_draw_string(0.0d0, 0.0d0, 50.0d0, 'Z')
      endif
      if (i_color .eq. 0) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      elseif (i_color .eq. 1) then
         call autogl_set_color(1.0d0, 1.0d0, 0.0d0) 
      elseif (i_color .eq. 2) then
         call autogl_set_color(0.0d0, 0.0d0, 1.0d0) 
      endif
      call autogl_draw_triangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     x_size, 0.0d0, 0.0d0, 
     c     0.0d0, y_size, 0.0d0)
      call autogl_draw_triangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, y_size, 0.0d0, 
     c     0.0d0, 0.0d0, z_size)
      call autogl_draw_triangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, 0.0d0, z_size, 
     c     x_size, 0.0d0, 0.0d0)
      call autogl_draw_triangle(
     c     x_size, 0.0d0, 0.0d0, 
     c     0.0d0, y_size, 0.0d0, 
     c     0.0d0, 0.0d0, z_size)
      end



c     描画アップデートコマンドのためのコールバックサブルーチン
c     "Draw"ボタンを選ぶと呼ばれる。 
      subroutine draw_button

      call autogl_draw_view()

      end



c     印刷コマンドのためのコールバックサブルーチン
c     "Print"ボタンを選ぶと呼ばれる。 
      subroutine print_button
      implicit none

      integer i_draw_axis, i_color
      real*8 x_size, y_size, z_size
      common i_draw_axis, i_color, x_size, y_size, z_size

      write(*,*) 'print'
      write(*,*) 'i_draw_axis', i_draw_axis
      write(*,*) 'i_color', i_color
      write(*,*) 'x_size, y_size, z_size', x_size, y_size, z_size

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

