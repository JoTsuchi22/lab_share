c     
c     メインルーチン
c     
      implicit none

      integer CELLS
      parameter (CELLS = 40)
      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      real*8 min_range, max_range, level_value, level_width
      integer model_is_rendered
      common grid_values, min_range, max_range, 
     c     level_value, level_width, model_is_rendered

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

      min_range = -1.0d0
      max_range = 5.0d0
      level_value = 2.0d0
      level_width = 2.0d0
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



      subroutine plot_isosurface_of_level(value)
      implicit none

      real*8 value

      integer CELLS
      parameter (CELLS = 40)
      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      real*8 min_range, max_range, level_value, level_width
      integer model_is_rendered
      common grid_values, min_range, max_range, 
     c     level_value, level_width, model_is_rendered

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



      subroutine plot_isosurface
      implicit none

      integer CELLS
      parameter (CELLS = 40)
      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      real*8 min_range, max_range, level_value, level_width
      integer model_is_rendered
      common grid_values, min_range, max_range, 
     c     level_value, level_width, model_is_rendered

      call plot_isosurface_of_level(level_value - level_width)
      call plot_isosurface_of_level(level_value)
      call plot_isosurface_of_level(level_value + level_width)

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン */
c     ビューが再表示されるごとに呼ばれる。 */
      subroutine autogl_redraw
      implicit none

      integer CELLS
      parameter (CELLS = 40)
      real*8 grid_values(CELLS + 1, CELLS + 1, CELLS + 1)
      real*8 min_range, max_range, level_value, level_width
      integer model_is_rendered
      common grid_values, min_range, max_range, 
     c     level_value, level_width, model_is_rendered

      real*8 size

c     画面左下に三次元座標系を描画する。 
      call autogl_draw_coordinate_system3d(20)

c     格子の枠を描画する。 
      size = 30.0d0
      call autogl_set_color(1.0d0, 1.0d0, 1.0d0)
      call autogl_draw_box3d(-size, -size, -size, size, size, size)

c     半透明表示を有効にする。 
      call autogl_turn_on_transparency()

c     まず、一度だけモデルをディスプレイリストへ描画しておく。 
      if (model_is_rendered .eq. 0) then

c     ディスプレイリストを開き、これを空にする。 
         call autogl_open_display_list()

         call plot_isosurface()

c     ディスプレイリストを閉じる。 
         call autogl_close_display_list()

         model_is_rendered = 1
      endif

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
