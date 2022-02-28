c     
c     メインルーチン
c     
      implicit none

      integer selected_face_id
      common selected_face_id

c     現在選択されている表面の番号 
      selected_face_id = -1



      call autogl_set_view_size(70.0d0)

c     イベント処理コールバック関数としてautogl_handle_eventを指定。 
      call autogl_set_default_callback_in_mode3d()

c     初期状態では、ユーザー定義イベント処理を行う。 
      call autogl_set_mode3d(0)

c     イベント処理中, ドラッグを有効にする。
      call autogl_enable_drag_in_mode3d()

c     パネル上に三次元アプリ用イベント処理機能を準備する。 
      call autogl_set_panel_in_mode3d()

c     AutoGLのメインループ
c     二度と戻らない
      call autogl_main()

      end



c     サブルーチンautogl_redraw
c     ビューの再描画のためのコールバックサブルーチン 
c     ビューが再表示されるごとに呼ばれる。 
      subroutine autogl_redraw
      implicit none

      integer selected_face_id
      common selected_face_id

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
      
      call autogl_draw_line(
     c     0.0d0, 0.0d0, -50.0d0,
     c     0.0d0, 0.0d0, 50.0d0)
      call autogl_draw_string(0.0d0, 0.0d0, 50.0d0, "Z = 50") 
      call autogl_draw_string(0.0d0, 0.0d0, -50.0d0, "Z = -50") 

      call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      call autogl_draw_string(0.0d0, 0.0d0, 0.0d0, "O")

      if (selected_face_id .eq. 1) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0)

      if (selected_face_id .eq. 2) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)

      if (selected_face_id .eq. 3) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 0.0d0, 
     c     0.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)

      if (selected_face_id .eq. 4) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     10.0d0, 0.0d0, 30.0d0)

      if (selected_face_id .eq. 5) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     0.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 0.0d0, 
     c     10.0d0, 0.0d0, 30.0d0, 
     c     0.0d0, 0.0d0, 30.0d0)

      if (selected_face_id .eq. 6) then
         call autogl_set_color(1.0d0, 0.0d0, 0.0d0) 
      else
         call autogl_set_color(1.0d0, 1.0d0, 1.0d0) 
      endif
      call autogl_draw_quadrangle(
     c     0.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 0.0d0, 
     c     10.0d0, 20.0d0, 30.0d0, 
     c     0.0d0, 20.0d0, 30.0d0)
      
      call autogl_set_color(0.0d0, 0.0d0, 1.0d0) 
      call autogl_draw_string(5.0d0, 10.0d0, -2.0d0, "1") 
      call autogl_draw_string(5.0d0, 10.0d0, 32.0d0, "2") 
      call autogl_draw_string(-2.0d0, 10.0d0, 15.0d0, "3") 
      call autogl_draw_string(12.0d0, 10.0d0, 15.0d0, "4") 
      call autogl_draw_string(5.0d0, -2.0d0, 15.0d0, "5") 
      call autogl_draw_string(5.0d0, 22.0d0, 15.0d0, "6") 

      end


      subroutine autogl_handle_event
      implicit none

      integer selected_face_id
      common selected_face_id

      integer event
      integer dc_x, dc_y, dc_z
      integer tolerance_dc
      integer key_char
      
      call autogl_get_view_event(event)

c     もし、マウスクリックがあれば 
      if (event .eq. 2) then 
         tolerance_dc = 10
         
         call autogl_get_pointing_device_position_dc(
     c        dc_x, dc_y, dc_z)

         call autogl_start_selection(dc_x, dc_y, tolerance_dc)

         call autogl_set_selection_id(1)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0)
         
         call autogl_set_selection_id(2)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 30.0d0, 
     c        10.0d0, 0.0d0, 30.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 20.0d0, 30.0d0)
         
         call autogl_set_selection_id(3)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 0.0d0, 
     c        0.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         
         call autogl_set_selection_id(4)
         call autogl_draw_quadrangle(
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        10.0d0, 0.0d0, 30.0d0)
         
         call autogl_set_selection_id(5)
         call autogl_draw_quadrangle(
     c        0.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 0.0d0, 
     c        10.0d0, 0.0d0, 30.0d0, 
     c        0.0d0, 0.0d0, 30.0d0)
         
         call autogl_set_selection_id(6)
         call autogl_draw_quadrangle(
     c        0.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 0.0d0, 
     c        10.0d0, 20.0d0, 30.0d0, 
     c        0.0d0, 20.0d0, 30.0d0)
         
         call autogl_end_selection()
         call autogl_get_selected_id(
     c        selected_face_id)

         write(*,*) ' selected face ID is ', selected_face_id

c     ビューを再描画する。 
         call autogl_draw_view()

c     もし、キーボード入力があれば 
      elseif (event .eq. 3) then
         call autogl_get_key_char(key_char)
         if (key_char .eq. ichar(' ')) then
            write(*,*) ' hit space '
         else
            call autogl_handle_default_key_event_in_mode3d()
         endif
      endif

      end



      subroutine autogl_idle
      end

      subroutine autogl_callback(id)
      implicit none
      integer id
      end
