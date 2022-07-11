program test_portion

implicit none
  integer, parameter      :: frames = 500
  integer, parameter      :: n=100
  real*8 , dimension(0:N) :: x, y
  real*8 , parameter      :: pi = 4.d0*atan(1.d0)
  real*8 , parameter      :: xmax = 4.0*pi, xmin = 0.d0
  real*8                  :: dx
  character(len=12)       :: data_name,frame_name
  integer                 :: i
  
  call creations()
  
  dx     = (xmax-xmin)/n
  x(0:N) = [(i*dx,i=0,N)]    !independent variable
  
  DO i = 1, frames
     write(data_name,770) 'data',i,'.dat'
     write(frame_name,770) 'plot',i,'.png'
     y = exp(-x*i*0.01/4.0)*sin(x*0.01*i)
     call generate_data_files(N,x,y,data_name)
     call gnuplot_each_frame(i,frame_name,data_name)
  ENDDO
  
  call system('ffmpeg -i Frames/plot%04d.png animation.avi')
  call system('mplayer animation.avi')
  
  call annihilation()
  
770 FORMAT(a,i4.4,a)  
end program test_portion


subroutine gnuplot_each_frame(i,frame_name,data_name)
  implicit none
  integer, intent(in)         :: i
  character(12), intent(in)   :: frame_name, data_name
  character(50)               :: output, data_file
  character(len=*),parameter  :: file_name = 'anima_plot.gp'
  integer                     :: fu
  
!   output    = 'Frames/'//frame_name
!   data_file = 'data/'//trim(data_name)
    
  open(newunit=fu,action='write',file=file_name,status='replace')
  write(fu,*)   'set term pngcairo font "Arial,12"'
  write(fu,770) "set output 'Frames/plot",i,".png'" 
  write(fu,*)   'set nokey'
  write(fu,*)   'set xtics nomirror'
  write(fu,*)   'set ytics nomirror'
  write(fu,*)   'set border lw 2.0'
  write(fu,*)   'set grid lw 2.0'
  write(fu,*)   'set xlabel "Independent Variable" '
  write(fu,*)   'set ylabel "exp(-x*iframe*0.01/4.0)*sin(x*0.01*iframe)" '
  write(fu,*)   'set xrange [0.0:2.0*pi]'
  write(fu,*)   'set yrange [-1.0:1.0]'
  write(fu,*)   'set title "Simple Animation" '
!   write(fu,770) "plot 'data/data",i,".dat' u 1:2 w lp lw 2.0 pt 7 ps 2 "   
  write(fu,770) "plot 'data/data",i,".dat' u 1:2 w lp lw 2.0 lc rgb 'red' pt 7 ps 2 "   
  close(fu)
  
  call execute_command_line('gnuplot -p anima_plot.gp')
  
770 FORMAT(1x,a,i4.4,a)  
return
end subroutine gnuplot_each_frame

subroutine creations()
!Here data and frames directory going to be created
!for a temporary storage of data and frames.  
 implicit none
  logical             :: there
  character(len=125)  :: dirName, Frames, data
  
  dirName = trim(data)
  
  
  inquire(exist=there,file='animation.avi')
  IF(there) call system('rm animation.avi')  !At the beginning old movie file is deleted
  
  inquire(exist=there,file=dirName//'/') !   // is for cancationation  
  IF(.not.there) call system('mkdir data')
  
  inquire(exist=there,file=trim(Frames)//'/.')  
  IF(.not.there) call system('mkdir Frames')
  
return
end subroutine creations

subroutine annihilation()
implicit none
  call system('rm -r data')
  call system('rm -r Frames')
return
end subroutine annihilation

subroutine generate_data_files(N,x,y,data_name)
!This subroutine is being called for each frames
implicit none
 integer,       intent(in)          :: N
 character(12), intent(in)          :: data_name
 real*8, dimension(0:N), intent(in) :: x(0:N), y(0:N)
 integer                            :: i, fu1
 
 open(newunit=fu1,action='write',file='data/'//data_name,status='replace')
 DO i = 0, N
    WRITE(fu1,*) x(i), y(i)
 ENDDO
 close(fu1)
 
return
end subroutine generate_data_files
