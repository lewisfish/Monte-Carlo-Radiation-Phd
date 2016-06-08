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

DOUBLE PRECISION, allocatable :: mus_array(:,:),mua_array(:,:)
DOUBLE PRECISION, allocatable :: nadh_array(:,:), tyrosine_array(:,:)
DOUBLE PRECISION, allocatable :: riboflavin_array(:,:)!,conc(:)

!fluorphores emission spectrum
DOUBLE PRECISION, allocatable :: nadh_fluro(:,:), nadh_cdf(:), ribo_fluro(:,:), ribo_cdf(:)
DOUBLE PRECISION, allocatable :: tyro_fluro(:,:), tyro_cdf(:)
integer, allocatable :: fluroexit(:),fluroexitGLOBAL(:)
end MODULE iarray
