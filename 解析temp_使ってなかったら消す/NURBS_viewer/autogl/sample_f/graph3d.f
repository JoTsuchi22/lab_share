c     
c     メインルーチン
c     
      implicit none

      integer CELLS
      parameter (CELLS = 100)
      real*8 grid_values(CELLS + 1, CELLS + 1)
      real*8 z_scale
      integer model_is_rendered
      common grid_values, z_scale, model_is_rendered

      integer i, j
      real*8 scale, x, y

      do i = 1, CELLS + 1
         do j = 1, CELLS + 1
            scale = 4.0d0 / CELLS
            x = (i - 1 - CELLS / 2) * scale
            y = (j - 1 - CELLS / 2) * scale
            grid_values(i, j) = sin(x * 3.1415) + y * y
         enddo
      enddo

      z_scale = 5.0d0
      model_is_rendered = 0

c     ビューの大きさ 
      call autogl_set_view_size(60.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode3d()
      call autogl_set_panel_in_mode3d()
      call autogl_enable_drag_in_mode3d()

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



      subroutine plot_graph3d
      implicit none

      integer CELLS
      parameter (CELLS = 100)
      real*8 grid_values(CELLS + 1, CELLS + 1)
      real*8 z_scale
      integer model_is_rendered
      common grid_values, z_scale, model_is_rendered

      real*8 cell_size
      integer i, j
      real*8 x00, y00, value00
      real*8 x10, y10, value10
      real*8 x11, y11, value11
      real*8 x01, y01, value01

      cell_size = 60.0d0 / CELLS

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0)

      do i = 1, CELLS
         do j = 1, CELLS
            x00 = (i - 1 - CELLS / 2) * cell_size
            y00 = (j - 1 - CELLS / 2) * cell_size
            value00 = grid_values(i + 0, j + 0)
            
            x10 = x00 + cell_size
            y10 = y00
            value10 = grid_values(i + 1, j + 0)
            
            x11 = x00 + cell_size
            y11 = y00 + cell_size
            value11 = grid_values(i + 1, j + 1)
            
            x01 = x00
            y01 = y00 + cell_size
            value01 = grid_values(i + 0, j + 1)
            
            call autogl_draw_triangle(
     c           x00, y00, value00 * z_scale,
     c           x10, y10, value10 * z_scale,
     c           x11, y11, value11 * z_scale)
            call autogl_draw_triangle(
     c           x11, y11, value11 * z_scale,
     c           x01, y01, value01 * z_scale,
     c           x00, y00, value00 * z_scale)
         enddo
      enddo

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン */
c     ビューが再表示されるごとに呼ばれる。 */
      subroutine autogl_redraw
      implicit none

      integer CELLS
      parameter (CELLS = 100)
      real*8 grid_values(CELLS + 1, CELLS + 1)
      real*8 z_scale
      integer model_is_rendered
      common grid_values, z_scale, model_is_rendered

c     まず、一度だけモデルをディスプレイリストへ描画しておく。 
      if (model_is_rendered .eq. 0) then

c     ディスプレイリストを開き、これを空にする。 
         call autogl_open_display_list()

         call plot_graph3d()

c     ディスプレイリストを閉じる。 
         call autogl_close_display_list()

         model_is_rendered = 1
      endif

c     鏡面反射をオンにする。
      call autogl_turn_on_specular()

c     画面左下に三次元座標系を描画する。 
      call autogl_draw_coordinate_system3d(20)

c     ディスプレイリストに入れた描画命令を実際に描画する。 
      call autogl_draw_display_list()

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
