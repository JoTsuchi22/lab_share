c     
c     メインルーチン
c     
      implicit none

      integer CELLS
      parameter (CELLS = 100)
      real*8 grid_values(CELLS + 1, CELLS + 1)
      real*8 min_range, max_range
      integer n_grades
      integer model_is_rendered
      common grid_values, min_range, max_range, n_grades, 
     c     model_is_rendered

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

      min_range = -1.0d0
      max_range = 5.0d0
      n_grades = 20
      model_is_rendered = 0
   
c     ビューの大きさ 
      call autogl_set_view_size(60.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode2d()
      call autogl_set_mode2d(2)
      call autogl_set_panel_in_mode2d()
      call autogl_enable_drag_in_mode2d()

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



      subroutine plot_line_contour
      implicit none

      integer CELLS
      parameter (CELLS = 100)
      real*8 grid_values(CELLS + 1, CELLS + 1)
      real*8 min_range, max_range
      integer n_grades
      integer model_is_rendered
      common grid_values, min_range, max_range, n_grades, 
     c     model_is_rendered

      real*8 cell_size
      integer i, j
      real*8 x00, y00, value00
      real*8 x10, y10, value10
      real*8 x11, y11, value11
      real*8 x01, y01, value01

      cell_size = 60.0d0 / CELLS

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

c     線コンタ-つきの四辺形を表示する。 
            call autogl_draw_line_contour_quadrangle(
     c           min_range, max_range, n_grades,
     c           x00, y00, 0.0d0, value00,
     c           x10, y10, 0.0d0, value10,
     c           x11, y11, 0.0d0, value11,
     c           x01, y01, 0.0d0, value01)
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
      real*8 min_range, max_range
      integer n_grades
      integer model_is_rendered
      common grid_values, min_range, max_range, n_grades, 
     c     model_is_rendered

c     まず、一度だけモデルをディスプレイリストへ描画しておく。 
      if (model_is_rendered .eq. 0) then

c     ディスプレイリストを開き、これを空にする。 
         call autogl_open_display_list()

         call plot_line_contour()

c     ディスプレイリストを閉じる。 
         call autogl_close_display_list()

         model_is_rendered = 1
      endif

c     画面右側にコンタ-のバーを描画する。 
      call autogl_draw_line_contour_map(
     c     80, min_range, max_range, n_grades)

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
