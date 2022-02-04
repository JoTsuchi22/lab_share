c     
c     メインルーチン
c     
      implicit none

c     格子の分割数 
      integer CELLS
      parameter (CELLS = 40)

c     格子上のスカラー値 
      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 level_value, level_width
      real*8 a, b, c, d
      common/parameter/ min_range, max_range, level_value, level_width, 
     c     a, b, c, d

      integer section_is_shown, isosurface_is_shown
      common/flag/ section_is_shown, isosurface_is_shown

      integer i, j, k
      real*8 scale, x, y, z

      do i = 1, CELLS + 1
         do j = 1, CELLS + 1
            do k = 1, CELLS + 1
               scale = 4.0d0 / CELLS
               x = (i - 1 - CELLS / 2) * scale
               y = (j - 1 - CELLS / 2) * scale
               z = (k - 1 - CELLS / 2) * scale
               grid_values(i, j, k) = sin(x * 3.1415) + y * y + z * z
            enddo
         enddo
      enddo

c     コンタ-の最大・最小値 
      min_range = -1.0d0
      max_range = 5.0d0

c     真ん中の等値面の値 
      level_value = 2.0d0

c     幅（他の二枚の等値面は、level_value +- level_width となる。） 
      level_width = 2.0d0

c     断面の方程式 Ax+By+Cz+D=0 の係数
      a = 1.0d0
      b = 1.0d0
      c = 1.0d0
      d = 0.0d0

c     断面を表示するか 
      section_is_shown = 0

c     等値面を表示するか 
      isosurface_is_shown = 0
   


c     ビューを準備する。 

c     ビューの大きさ 
      call autogl_set_view_size(60.0d0)

c     おまじない 
      call autogl_set_default_callback_in_mode3d()
      call autogl_enable_drag_in_mode3d()

c     パネルを準備する。 
c     ボタン用サブルーチン番号やコントロール用変数を登録する。
      call autogl_add_group('model')
      call autogl_set_panel_in_mode3d()
      call autogl_add_comment()
      call autogl_set_label('Model Parameters')

c     断面コンタ- 
      call autogl_add_boolean(section_is_shown, 'section_is_shown')
      call autogl_set_label('section')
      call autogl_add_comment()
      call autogl_set_label('equation Ax+By+Cz+D=0')
      call autogl_add_real(a, 'a')
      call autogl_add_real(b, 'b')
      call autogl_add_real(c, 'c')
      call autogl_add_real(d, 'd')
      call autogl_add_real(min_range, 'min_range')
      call autogl_add_real(max_range, 'max_range')

c     等値面 
      call autogl_add_boolean(isosurface_is_shown, 
     c     'isosurface_is_shown')
      call autogl_set_label('iso-surface')
      call autogl_add_real(level_value, 'level_value')
      call autogl_add_real(level_width, 'level_width')

      call autogl_set_button(0, 'UpdateModel')

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     格子の断面コンタ-を描画する。
      subroutine plot_section
      implicit none

      integer CELLS
      parameter (CELLS = 40)

      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 level_value, level_width
      real*8 a, b, c, d
      common/parameter/ min_range, max_range, level_value, level_width, 
     c     a, b, c, d

      real*8 cell_size
      integer gx, gy, gz
      real*8 x000, y000, z000, value000
      real*8 x100, y100, z100, value100
      real*8 x110, y110, z110, value110
      real*8 x010, y010, z010, value010
      real*8 x001, y001, z001, value001
      real*8 x101, y101, z101, value101
      real*8 x111, y111, z111, value111
      real*8 x011, y011, z011, value011

      cell_size = 60.0d0 / CELLS

      do gx = 1, CELLS
         do gy = 1, CELLS
            do gz = 1, CELLS
               x000 = (gx - 1 - CELLS / 2) * cell_size
               y000 = (gy - 1 - CELLS / 2) * cell_size
               z000 = (gz - 1 - CELLS / 2) * cell_size
               value000 = grid_values(gx, gy, gz)

               x100 = x000 + cell_size
               y100 = y000
               z100 = z000
               value100 = grid_values(gx + 1, gy, gz)

               x110 = x000 + cell_size
               y110 = y000 + cell_size
               z110 = z000
               value110 = grid_values(gx + 1, gy + 1, gz)

               x010 = x000
               y010 = y000 + cell_size
               z010 = z000
               value010 = grid_values(gx, gy + 1, gz)

               x001 = x000
               y001 = y000
               z001 = z000 + cell_size
               value001 = grid_values(gx, gy, gz + 1)

               x101 = x000 + cell_size
               y101 = y000
               z101 = z000 + cell_size
               value101 = grid_values(gx + 1, gy, gz + 1)

               x111 = x000 + cell_size
               y111 = y000 + cell_size
               z111 = z000 + cell_size
               value111 = grid_values(gx + 1, gy + 1, gz + 1)

               x011 = x000
               y011 = y000 + cell_size
               z011 = z000 + cell_size
               value011 = grid_values(gx, gy + 1, gz + 1)

c     六面体の断面コンターを描画する。 
               call autogl_draw_section_smooth_contour_hexahedron(
     c              a, b, c, d, min_range, max_range,
     c              x000, y000, z000, value000,
     c              x100, y100, z100, value100,
     c              x110, y110, z110, value110,
     c              x010, y010, z010, value010,
     c              x001, y001, z001, value001,
     c              x101, y101, z101, value101,
     c              x111, y111, z111, value111,
     c              x011, y011, z011, value011)
            enddo
         enddo
      enddo

      end



c     格子の等値面を描画する。 
      subroutine plot_isosurface_of_level(value)
      implicit none

      real*8 value

      integer CELLS
      parameter (CELLS = 40)

      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      common/value/ grid_values

      real*8 min_range, max_range
      real*8 level_value, level_width
      real*8 a, b, c, d
      common/parameter/ min_range, max_range, level_value, level_width, 
     c     a, b, c, d

      real*8 cell_size
      real*8 grade
      real*8 red, green, blue
      integer gx, gy, gz
      real*8 x000, y000, z000, value000
      real*8 x100, y100, z100, value100
      real*8 x110, y110, z110, value110
      real*8 x010, y010, z010, value010
      real*8 x001, y001, z001, value001
      real*8 x101, y101, z101, value101
      real*8 x111, y111, z111, value111
      real*8 x011, y011, z011, value011

      cell_size = 60.0d0 / CELLS
      grade = (value - min_range) / (max_range - min_range)
      call autogl_get_contour_color(red, green, blue,
     c     grade)
      call autogl_set_color(red, green, blue)

      do gx = 1, CELLS
         do gy = 1, CELLS
            do gz = 1, CELLS
               x000 = (gx - 1 - CELLS / 2) * cell_size
               y000 = (gy - 1 - CELLS / 2) * cell_size
               z000 = (gz - 1 - CELLS / 2) * cell_size
               value000 = grid_values(gx, gy, gz)

               x100 = x000 + cell_size
               y100 = y000
               z100 = z000
               value100 = grid_values(gx + 1, gy, gz)

               x110 = x000 + cell_size
               y110 = y000 + cell_size
               z110 = z000
               value110 = grid_values(gx + 1, gy + 1, gz)

               x010 = x000
               y010 = y000 + cell_size
               z010 = z000
               value010 = grid_values(gx, gy + 1, gz)

               x001 = x000
               y001 = y000
               z001 = z000 + cell_size
               value001 = grid_values(gx, gy, gz + 1)

               x101 = x000 + cell_size
               y101 = y000
               z101 = z000 + cell_size
               value101 = grid_values(gx + 1, gy, gz + 1)

               x111 = x000 + cell_size
               y111 = y000 + cell_size
               z111 = z000 + cell_size
               value111 = grid_values(gx + 1, gy + 1, gz + 1)

               x011 = x000
               y011 = y000 + cell_size
               z011 = z000 + cell_size
               value011 = grid_values(gx, gy + 1, gz + 1)

c     六面体の等値面を描画する。 
               call autogl_draw_isosurface_hexahedron(
     c              value, value, 0,
     c              x000, y000, z000, value000,
     c              x100, y100, z100, value100,
     c              x110, y110, z110, value110,
     c              x010, y010, z010, value010,
     c              x001, y001, z001, value001,
     c              x101, y101, z101, value101,
     c              x111, y111, z111, value111,
     c              x011, y011, z011, value011)
            enddo
         enddo
      enddo

      end



c     等値面を描画する。 
      subroutine plot_isosurface
      implicit none

      real*8 min_range, max_range
      real*8 level_value, level_width
      real*8 a, b, c, d
      common/parameter/ min_range, max_range, level_value, level_width, 
     c     a, b, c, d

      call plot_isosurface_of_level(level_value - level_width)
      call plot_isosurface_of_level(level_value)
      call plot_isosurface_of_level(level_value + level_width)

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      real*8 min_range, max_range
      real*8 level_value, level_width
      real*8 a, b, c, d
      common/parameter/ min_range, max_range, level_value, level_width, 
     c     a, b, c, d

      integer section_is_shown, isosurface_is_shown
      common/flag/ section_is_shown, isosurface_is_shown

      real*8 size

c     画面左下に三次元座標系を描画する。 
      call autogl_draw_coordinate_system3d(20)

      if (section_is_shown .eq. 1) then
c     画面右側にコンタ-のバーを描画する。 
         call autogl_draw_contour_map(80, min_range, max_range)
      endif
            
c     格子の枠を描画する。 
      size = 30.0d0
      call autogl_set_color(1.0d0, 1.0d0, 1.0d0)
      call autogl_draw_box3d(-size, -size, -size, size, size, size)

c     半透明表示を有効にする。 
      call autogl_turn_on_transparency()

c     ディスプレイリストに入れた描画命令を実際に描画する。 
      call autogl_draw_display_list()

      end



c     描画アップデートコマンドのためのコールバックサブルーチン
c     "UpdateModel"ボタンを選ぶと呼ばれる。 
      subroutine update_model_button
      implicit none

      integer section_is_shown, isosurface_is_shown
      common/flag/ section_is_shown, isosurface_is_shown

c     ディスプレイリストを開き、これを空にする。 
      call autogl_open_display_list()

c     以下の描画命令をディスプレイリストに順次入れていく。 
c     実際の描画はここではまだ行われない。 

      if (section_is_shown .eq. 1) then
         call plot_section()
      endif

      if (isosurface_is_shown .eq. 1) then
         call plot_isosurface()
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

