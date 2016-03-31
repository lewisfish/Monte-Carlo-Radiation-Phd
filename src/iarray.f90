MODULE iarray
!
!  Contains all array var names.
!
implicit none
save

DOUBLE PRECISION, allocatable :: xface(:),yface(:),zface(:)
DOUBLE PRECISION, allocatable :: rhokap(:,:,:,:)
DOUBLE PRECISION, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

DOUBLE PRECISION, allocatable :: Carotene_array(:,:), Bilirubin_array(:,:)
DOUBLE PRECISION, allocatable :: Deoxy_Hb_array(:,:), Oxy_Hb_array(:,:)

DOUBLE PRECISION, allocatable :: Carotene_cdf(:), Bilirubin_cdf(:)
DOUBLE PRECISION, allocatable :: Oxy_Hb_cdf(:), Deoxy_Hb_cdf(:) 

!real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)
!real, allocatable :: fluro_pos(:,:,:),fluro_posGLOBAL(:,:,:)

!real, allocatable :: image(:,:,:),deposit(:,:,:),dep(:)
!real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
!real, allocatable :: depositGLOBAL(:,:,:),depGLOBAL(:)

DOUBLE PRECISION, allocatable :: mus_array(:,:),mua_array(:,:)
DOUBLE PRECISION, allocatable :: f_cdf_c(:),e_cdf_c(:),excite_array_c(:,:),fluro_array_c(:,:)
DOUBLE PRECISION, allocatable :: f_cdf_n(:),e_cdf_n(:),excite_array_n(:,:),fluro_array_n(:,:)
integer, allocatable :: fluroexit(:),fluroexitGLOBAL(:)
end MODULE iarray
