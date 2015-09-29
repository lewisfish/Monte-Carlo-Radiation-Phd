      subroutine gridedge(x1,y1,z1,zmax,dist,v,dist2)
      
      implicit none
      
      real v(3),x1,x2,y1,y2,z1,z2,x,y,z,t,zmax,dist
      real dist2
      
      
      !set postion to peel off to.
      x2=0.
      y2=0.
      z2=zmax+5.
      
      
      !eqn of line r_vec=r1+t*v_vec
      !create v vector
      v(1)=x2-x1      
      v(2)=y2-y1      
      v(3)=z2-z1 
      
      t=zmax-z1/v(3)
      y=(v(2)*t)+y1
      x=(v(1)*t)+x1
      
      dist=sqrt((x-x1)**2+(y-y1)**2+(zmax-z1)**2)
      dist2=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
  
      end subroutine gridedge
