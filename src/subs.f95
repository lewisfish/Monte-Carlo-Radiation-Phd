MODULE subs

implicit none
save

CONTAINS
   SUBROUTINE directory
   use constants,only : cwd,homedir,fileplace,resdir
   implicit none

   integer :: io

   !get current working directory
   CALL getcwd(cwd)

   !get 'home' dir from cwd
   homedir=trim(cwd(1:index(trim(cwd),'/bin')))
   !get data dir
   fileplace=trim(homedir)//'data/'
   
   !checks to see if data folder exists, if not creates it.
   call chdir(fileplace,io)
   if(io.ne.0)then
      print*,'directory does not exist'
      print*, 'creating directory'
      call system("mkdir"//fileplace)
   end if
   
   !get res dir
   resdir=trim(homedir)//'res/'
   
   end SUBROUTINE directory

end MODULE subs
