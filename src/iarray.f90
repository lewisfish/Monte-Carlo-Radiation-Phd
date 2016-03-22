MODULE iarray
!
!  Contains all array var names.
!
implicit none
save

DOUBLE PRECISION, allocatable :: xface(:),yface(:),zface(:)
DOUBLE PRECISION, allocatable :: rhokap(:,:,:,:)
DOUBLE PRECISION, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

!real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)
!real, allocatable :: fluro_pos(:,:,:),fluro_posGLOBAL(:,:,:)

!real, allocatable :: image(:,:,:),deposit(:,:,:),dep(:)
!real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
!real, allocatable :: depositGLOBAL(:,:,:),depGLOBAL(:)

DOUBLE PRECISION, allocatable :: f_cdf(:),mua_array(:,:),excite_array(:,:)
DOUBLE PRECISION, allocatable :: mus_array(:,:),fluro_array(:,:),e_cdf(:)
integer, allocatable :: fluroexit(:),fluroexitGLOBAL(:)
end MODULE iarray
