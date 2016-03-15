MODULE photon_vars
!
! Module containing photon vars:
!           %p:current position of photon
!           n%p:current direction of photon
!           sin/cos(%)/phi:various angles related to photons flight 
!           angles measured from ''north'' in thetas case.

implicit none
save

real :: xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi

end MODULE photon_vars
