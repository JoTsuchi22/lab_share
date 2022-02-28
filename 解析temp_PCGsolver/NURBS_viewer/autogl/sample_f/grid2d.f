c     
c     メインルーチン
c     
      implicit none

c     格子の分割数 
      integer CELLS
      parameter (CELLS = 100)

c     格子上のスカラー値 
      real*8 grid_values(CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 z_scale
      integer n_grades, contour_type
      common/parameter/ min_range, max_range, z_scale,
     c     contour_type, n_grades

      integer grid_is_shown, contour_is_shown, graph3d_is_shown
      common/flag/ grid_is_shown, contour_is_shown, graph3d_is_shown

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

c     コンタ-の最大・最小値、グレード数、表示形式の種類
      min_range = -1.0d0
      max_range = 5.0d0
      n_grades = 20
      contour_type = 0

c     三次元グラフのZ方向の拡大率 
      z_scale = 5.0d0

c     格子を表示するか 
      grid_is_shown = 0

c     コンタ-を表示するか 
      contour_is_shown = 0

c     三次元グラフを表示するか 
      graph3d_is_shown = 0
   


c     ビューを準備する。 

c     ビューの背景色 
      call autogl_set_background_color(0.7d0, 0.7d0, 0.7d0)
c     ビューの大きさ 
      call autogl_set_view_size(80.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode3d()
      call autogl_enable_drag_in_mode3d()

c     パネルを準備する。 
c     ボタン用サブルーチン番号やコントロール用変数を登録する。
      call autogl_add_group('model')
      call autogl_set_panel_in_mode3d()

c     コメントの登録 : ラベル 
      call autogl_add_comment()
      call autogl_set_label('Model Parameters')

c     格子 
      call autogl_add_boolean(grid_is_shown, 'grid_is_shown')
      call autogl_set_label('grid')

c     コンタ- 
      call autogl_add_boolean(contour_is_shown, 'contour_is_shown')
      call autogl_set_label('contour')
      call autogl_add_integer(contour_type, 'contour_type')
      call autogl_set_label('type')
      call autogl_add_integer_item('smooth')
      call autogl_add_integer_item('line')
      call autogl_add_integer_item('band')
      call autogl_add_real(min_range, 'min_range')
      call autogl_set_label('min')
      call autogl_add_real(max_range, 'max_range')
      call autogl_set_label('max')
      call autogl_add_integer(n_grades, 'n_grades')
      call autogl_set_label('num. grades')

c     三次元グラフ 
      call autogl_add_boolean(graph3d_is_shown, 'graph3d_is_shown')
      call autogl_set_label('3-D graph')
      call autogl_add_real(z_scale, 'z_scale')

      call autogl_set_button(0, 'UpdateModel')

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     格子を描画する。 
      subroutine draw_grid
      implicit none

      integer CELLS
      parameter (CELLS = 100)

      real*8 grid_values(CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 z_scale
      integer n_grades, contour_type
      common/parameter/ min_range, max_range, z_scale,
     c     contour_type, n_grades

      real*8 cell_size
      integer i, j
      real*8 x00, y00
      real*8 x10, y10
      real*8 x11, y11
      real*8 x01, y01

      cell_size = 60.0d0 / CELLS

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0)

      do i = 1, CELLS
         do j = 1, CELLS
            x00 = (i - 1 - CELLS / 2) * cell_size
            y00 = (j - 1 - CELLS / 2) * cell_size

            x10 = x00 + cell_size
            y10 = y00

            x11 = x00 + cell_size
            y11 = y00 + cell_size

            x01 = x00
            y01 = y00 + cell_size

            call autogl_draw_line(
     c           x00, y00, 0.0d0,
     c           x10, y10, 0.0d0)
            call autogl_draw_line(
     c           x10, y10, 0.0d0,
     c           x11, y11, 0.0d0)
            call autogl_draw_line(
     c           x11, y11, 0.0d0,
     c           x01, y01, 0.0d0)
            call autogl_draw_line(
     c           x01, y01, 0.0d0,
     c           x00, y00, 0.0d0)
         enddo
      enddo

      end



c     コンタ-を描画する。 
      subroutine draw_contour
      implicit none

      integer CELLS
      parameter (CELLS = 100)

      real*8 grid_values(CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 z_scale
      integer n_grades, contour_type
      common/parameter/ min_range, max_range, z_scale,
     c     contour_type, n_grades

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

            if (contour_type .eq. 0) then
c     スムーズコンタ-つきの四辺形を表示する。 
               call autogl_draw_smooth_contour_quadrangle(
     c              min_range, max_range,
     c              x00, y00, 0.0d0, value00,
     c              x10, y10, 0.0d0, value10,
     c              x11, y11, 0.0d0, value11,
     c              x01, y01, 0.0d0, value01)
            else if (contour_type .eq. 1) then
c     線コンタ-つきの四辺形を表示する。 
               call autogl_draw_line_contour_quadrangle(
     c              min_range, max_range, n_grades,
     c              x00, y00, 0.0d0, value00,
     c              x10, y10, 0.0d0, value10,
     c              x11, y11, 0.0d0, value11,
     c              x01, y01, 0.0d0, value01)
            else if (contour_type .eq. 2) then
c     バンドコンタ-つきの四辺形を表示する。 
               call autogl_draw_band_contour_quadrangle(
     c              min_range, max_range, n_grades,
     c              x00, y00, 0.0d0, value00,
     c              x10, y10, 0.0d0, value10,
     c              x11, y11, 0.0d0, value11,
     c              x01, y01, 0.0d0, value01)
            endif
         enddo
      enddo

      end



c     三次元グラフを描画する。 
      subroutine draw_graph3d
      implicit none

      integer CELLS
      parameter (CELLS = 100)

      real*8 grid_values(CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 z_scale
      integer n_grades, contour_type
      common/parameter/ min_range, max_range, z_scale,
     c     contour_type, n_grades

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
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      real*8 min_range, max_range
      real*8 z_scale
      integer n_grades, contour_type
      common/parameter/ min_range, max_range, z_scale,
     c     contour_type, n_grades

      integer grid_is_shown, contour_is_shown, graph3d_is_shown
      common/flag/ grid_is_shown, contour_is_shown, graph3d_is_shown

c     画面左下に三次元座標系を描画する。 
      call autogl_draw_coordinate_system3d(20)

      if (contour_is_shown .eq. 1) then
c     画面右側にコンタ-のバーを描画する。 
         if (contour_type .eq. 0) then
            call autogl_draw_contour_map(80, min_range, max_range)
         else if (contour_type .eq. 1) then
            call autogl_draw_line_contour_map
     c           (80, min_range, max_range, n_grades)
         else if (contour_type .eq. 2) then
            call autogl_draw_band_contour_map
     c           (80, min_range, max_range, n_grades)
         endif
      endif

c     ディスプレイリストに入れた描画命令を実際に描画する。 
      call autogl_draw_display_list()

      end



c     描画アップデートコマンドのためのコールバックサブルーチン
c     "UpdateModel"ボタンを選ぶと呼ばれる。 
      subroutine update_model_button
      implicit none

      integer grid_is_shown, contour_is_shown, graph3d_is_shown
      common/flag/ grid_is_shown, contour_is_shown, graph3d_is_shown

c     ディスプレイリストを開き、これを空にする。 
      call autogl_open_display_list()

c     以下の描画命令をディスプレイリストに順次入れていく。 
c     実際の描画はここではまだ行われない。 

      if (grid_is_shown .eq. 1) then
         call draw_grid()
      endif

      if (contour_is_shown .eq. 1) then
         call draw_contour()
      endif

      if (graph3d_is_shown .eq. 1) then
         call draw_graph3d()
      endif

c     ディスプレイリストを閉じる。 
      call autogl_close_display_list()

c     ビューを再描画する。
      call autogl_draw_view()
c     登録されたビューの再描画サブルーチン（この場合はautogl_redraw）が呼ばれる。 

      end



c     コールバックサブルーチン
      subroutine autogl_callback(id)
      implicit none

      integer id

      if (id .eq. 0) then
         call update_model_button()
      endif

      end



      subroutine autogl_handle_event
      end

      subroutine autogl_idle
      end

