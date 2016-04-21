MODULE iarray
!
!  Contains all array var names.
!
implicit none
save

DOUBLE PRECISION, allocatable :: xface(:),yface(:),zface(:)
DOUBLE PRECISION, allocatable :: rhokap(:,:,:,:),albedo(:,:,:,:)
DOUBLE PRECISION, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)
DOUBLE PRECISION, allocatable :: rgb(:,:,:),rgbGLOBAL(:,:,:)

!skin absorbers
DOUBLE PRECISION, allocatable :: Carotene_array(:,:), Bilirubin_array(:,:)
DOUBLE PRECISION, allocatable :: Deoxy_Hb_array(:,:), Oxy_Hb_array(:,:)
DOUBLE PRECISION, allocatable :: water_array(:,:), fat_array(:,:)

!fluorphores optical properties
DOUBLE PRECISION, allocatable :: fad_array(:,:), nad_array(:,:), nadh_array(:,:)
DOUBLE PRECISION, allocatable :: tryptophan_array(:,:), tyrosine_array(:,:)
DOUBLE PRECISION, allocatable :: riboflavin_array(:,:), conc(:,:)

!fluorphores emission spectrum
DOUBLE PRECISION, allocatable :: fad_fluro(:,:), fad_cdf(:)!, nad_fluro(:,:), nad_cdf(:)
DOUBLE PRECISION, allocatable :: nadh_fluro(:,:), nadh_cdf(:), ribo_fluro(:,:), ribo_cdf(:)
DOUBLE PRECISION, allocatable :: tyro_fluro(:,:), tyro_cdf(:), try_fluro(:,:), try_cdf(:)

!real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)
!real, allocatable :: fluro_pos(:,:,:),fluro_posGLOBAL(:,:,:)

!real, allocatable :: image(:,:,:),deposit(:,:,:),dep(:)
!real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
!real, allocatable :: depositGLOBAL(:,:,:),depGLOBAL(:)

integer, allocatable :: fluroexit(:,:),fluroexitGLOBAL(:,:)
end MODULE iarray
