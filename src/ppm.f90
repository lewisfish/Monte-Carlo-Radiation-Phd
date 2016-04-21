MODULE ppm

implicit none

CONTAINS
   
   subroutine write_ppm(filename, array)
   
   implicit none
   
   DOUBLE PRECISION, intent(IN) :: array(:, :, :)
   character(*),     intent(IN) :: filename
   integer                      :: i, j, io, w, h
   
   w = size(array, 1)
   h = size(array, 2)
   open(20, file = trim(filename), iostat = io)
   if(io/=0)then
      print*,'Cant open file'//trim(filename)
      print*,'Exiting write_ppm'
   else
      write(20, '(a2)') 'P3'
      write(20, '(i0, " ",i0)') w, h
      write(20, '(i0)') 255
      
      do i = h, 1, -1
         do j = 1, w
!            print*,int(array(i,j,:))
            write(20, *) int(array(i, j, :))
!            write(20, '(3a1)', advance='no') achar(int(array(i, j, :)))
            
         end do
      end do
      close(20)
   end if
   end subroutine write_ppm
   
   
   subroutine wave_to_RGB(wave, xcell, ycell , array, gamm)
   
      use constants, only : nxg, nyg
      
      implicit none
      
      DOUBLE PRECISION, intent(INOUT) :: wave, array(nxg, nyg, 3)
      integer,          intent(IN) :: xcell, ycell
      DOUBLE PRECISION             :: SSS, gamm, R, G, B
            
            if((wave .lt. 380.d0) .and. (wave .gt. 780.0d0))then
               R = 0.0d0
               G = 0.0d0
               B = 0.0d0
            end if
            if ((wave .ge. 380.d0).and.(wave .le. 440.d0))then 
              R = -1.0d0 * (wave - 440.0d0) / (440.0d0 - 380.0d0)
              G = 0.0d0
              B = 1.0d0
            end if
            if ((wave .ge. 440.0d0) .and. (wave .le. 490.0d0))then
              R = 0.0d0
              G = (wave - 440.0d0) / (490.0d0 - 440.0d0)
              B = 1.0d0
            end if
            if ((wave .ge. 490.0d0) .and. (wave .le. 510.0d0))then 
              R = 0.0d0
              G = 1.0d0
              B = -1.0d0 * (wave - 510.0d0) / (510.0d0 - 490.0d0)
            end if
            if ((wave .ge. 510.0d0) .and. (wave .le. 580.0d0))then 
              R = (wave - 510.0d0) / (580.0d0 - 510.0d0)
              G = 1.0d0
              B = 0.0d0
            end if
            if ((wave .ge. 580.0d0) .and. (wave .le. 645.0d0))then
              R = 1.0d0
              G = -1.0d0 * (wave - 645.0d0) / (645.0d0 - 580.0d0)
              B = 0.0d0
            end if
            if ((wave .ge. 645.0d0) .and. (wave .le. 780.0d0)) then
              R = 1.0d0
              G = 0.0d0
              B = 0.0d0
            end if
            
!      intensity falls off near to vision edge
         if (wave .gt. 700.0d0)then
            SSS = 0.3d0 + 0.7d0 * (780.0d0 - wave) / (780.0d0 - 700.0d0)
         else if (wave .lt. 420.0d0)then
            SSS = 0.3d0 + 0.7d0 * (wave - 380.0d0) / (420.0d0 - 380.0d0)
         else
            SSS = 1.0d0
         end if

!      gamma adjust and write image to an array
         array(xcell, ycell , 1) = array(xcell, ycell , 1) + (SSS * R)**gamm
         array(xcell, ycell , 2) = array(xcell, ycell , 2) + (SSS * G)**gamm
         array(xcell, ycell , 3) = array(xcell, ycell , 3) + (SSS * B)**gamm
   
   end subroutine wave_to_RGB
end MODULE ppm
