c     
c     メインルーチン
c     
      implicit none

      real*8 vertexs(6, 2)
      integer selected_vertex_id
      common vertexs, selected_vertex_id

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

c     現在選択されている頂点の番号 
      selected_vertex_id = -1



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
      integer selected_vertex_id
      common vertexs, selected_vertex_id

      integer i
      integer mark_size_dc
      integer is_filled

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

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      do i = 1, 6 - 1
         call autogl_draw_line(
     c        vertexs(i, 1), vertexs(i, 2), 0.0d0,
     c        vertexs(i + 1, 1), vertexs(i + 1, 2), 0.0d0)
      enddo
      
c     各頂点について、それが選択されているかどうかを描画する。 
      call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      do i = 1, 6
         is_filled = 0
         if (i .eq. selected_vertex_id) then
            is_filled = 1
         endif
         
         mark_size_dc = 10
c     ダイアモンドのマークを表示する。 
         call autogl_draw_diamond_mark2d (
     c        vertexs(i, 1), vertexs(i, 2), 0.0d0,
     c        mark_size_dc, is_filled)
c     mark_size_dcにより、10 pixelぐらいの大きさとなる。 
c     is_filled .eq. 1なら、この頂点は選択されている。 
      enddo

      end



      subroutine autogl_handle_event
      implicit none

      real*8 vertexs(6, 2)
      integer selected_vertex_id
      common vertexs, selected_vertex_id

      integer event
      integer i
      integer tolerance_dc
      integer flag
      
      call autogl_get_view_event(event)

c     もし、マウスクリックがあれば 
      if (event .eq. 2) then 
         
c     各頂点について、マウスポインタがヒットしたかどうか調べる。 
         selected_vertex_id = -1
         do i = 1, 6
            tolerance_dc = 10
            
c     この頂点がマウスポインタのそばにあれば 
            call autogl_pointing_device_is_hit2d(
     c           flag, 
     c           vertexs(i, 1), vertexs(i, 2),
     c           tolerance_dc)
c     点の座標は世界座標（ただし二次元）で与える。 
            if (flag .eq. 1) then
               selected_vertex_id = i
            endif
         enddo
         
         if (selected_vertex_id .ne. -1) then
            write(*,*) ' vertex ', selected_vertex_id, ' is selected '
         endif

c     ビューを再描画する。 
         call autogl_draw_view()

c     もし、キーボード入力があれば 
      elseif (event .eq. 3) then
         call autogl_handle_default_key_event_in_mode2d()
      endif

      end



      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
