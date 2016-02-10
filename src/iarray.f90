MODULE iarray
!
!  Contains all array var names.
!
implicit none
save

real, allocatable :: xface(:),yface(:),zface(:)
real, allocatable :: rhokap(:,:,:,:)
real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)
real, allocatable :: fluro_pos(:,:,:),fluro_posGLOBAL(:,:,:)

real, allocatable :: image(:,:,:),deposit(:,:,:),dep(:)
real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
real, allocatable :: depositGLOBAL(:,:,:),depGLOBAL(:)

real, allocatable :: f_cdf(:),mua_array(:,:),excite_array(:,:)
real, allocatable :: mus_array(:,:),fluro_array(:,:),e_cdf(:)
integer, allocatable :: fluroexit(:),fluroexitGLOBAL(:)
end MODULE iarray
