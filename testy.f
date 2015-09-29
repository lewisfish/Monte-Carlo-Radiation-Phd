      program testy

      real v(3),x,y,z

      v=0.
      x=1.
      y=7.
      z=99.

      v=x+y+z

      print*, v

      v(1)=x
      v(2)=y
      v(3)=z
      print*,v+x*v
      end program
