c     
c     メインルーチン
c     
      implicit none

      real*8 vertexs(6, 2)
      integer selected_edge_id
      common vertexs, selected_edge_id

c     頂点の座標 
      vertexs(1, 1) = 10.0d0
      vertexs(1, 2) = 0.0d0

      vertexs(2, 1) = 10.0d0
      vertexs(2, 2) = 20.0d0

      vertexs(3, 1) = 0.0d0
      vertexs(3, 2) = 20.0d0

      vertexs(4, 1) = -10.0d0
      vertexs(4, 2) = 0.0d0

      vertexs(5, 1) = 0.0d0
      vertexs(5, 2) = -20.0d0

      vertexs(6, 1) = 30.0d0
      vertexs(6, 2) = -20.0d0

c     現在選択されている稜線の番号 
      selected_edge_id = -1



      call autogl_set_view_size(70.0d0)

c     イベント処理コールバック関数としてautogl_handle_eventを指定。 
      call autogl_set_default_callback_in_mode2d()

c     初期状態では、ユーザー定義イベント処理を行う。 
      call autogl_set_mode2d(0)

c     イベント処理中, ドラッグを有効にする。
      call autogl_enable_drag_in_mode2d()

c     パネル上に二次元アプリ用イベント処理機能を準備する。 
      call autogl_set_panel_in_mode2d()

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     サブルーチンautoglredraw
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      real*8 vertexs(6, 2)
      integer selected_edge_id
      common vertexs, selected_edge_id

      integer i

      call autogl_set_color(0.0d0, 1.0d0, 0.0d0) 

      call autogl_draw_line(
     c     -50.0d0, 0.0d0, 0.0d0,
     c     50.0d0, 0.0d0, 0.0d0)
      call autogl_draw_string(50.0d0, 0.0d0, 0.0d0, "X = 50") 
      call autogl_draw_string(-50.0d0, 0.0d0, 0.0d0, "X = -50") 
         
      call autogl_draw_line(
     c     0.0d0, -50.0d0, 0.0d0,
     c     0.0d0, 50.0d0, 0.0d0)
      call autogl_draw_string(0.0d0, 50.0d0, 0.0d0, "Y = 50") 
      call autogl_draw_string(0.0d0, -50.0d0, 0.0d0, "Y = -50") 

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, "O")

      do i = 1, 6 - 1
         if (i .eq. selected_edge_id) then
            call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
         else
            call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
         endif
         call autogl_draw_line(
     c        vertexs(i, 1), vertexs(i, 2), 0.0d0,
     c        vertexs(i + 1, 1), vertexs(i + 1, 2), 0.0d0)
      enddo

      end



      subroutine autogl_handle_event
      implicit none

      real*8 vertexs(6, 2)
      integer selected_edge_id
      common vertexs, selected_edge_id

      integer event
      integer dc_x, dc_y, dc_z
      integer tolerance_dc
      integer i
      integer key_char
      
      call autogl_get_view_event(event)

c     もし、マウスクリックがあれば 
      if (event .eq. 2) then 
         
         tolerance_dc = 10

         call autogl_get_pointing_device_position_dc(
     c        dc_x, dc_y, dc_z)
         
         call autogl_start_selection(dc_x, dc_y, tolerance_dc)

         do i = 1, 6 - 1
            call autogl_set_selection_id(i)
            call autogl_draw_line(
     c           vertexs(i, 1), vertexs(i, 2), 0.0d0,
     c           vertexs(i + 1, 1), vertexs(i + 1, 2), 0.0d0)
         enddo
         
         call autogl_end_selection()
         call autogl_get_selected_id(
     c        selected_edge_id)
         
         write(*,*) ' selected edge ID is ', selected_edge_id

c     ビューを再描画する。 
         call autogl_draw_view()

c     もし、キーボード入力があれば 
      elseif (event .eq. 3) then
         call autogl_get_key_char(key_char)
         if (key_char .eq. ichar(' ')) then
            write(*,*) ' hit space '
         else
            call autogl_handle_default_key_event_in_mode2d()
         endif
      endif
      end



      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
