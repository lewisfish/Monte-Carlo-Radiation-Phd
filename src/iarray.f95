MODULE iarray

implicit none
save

real, allocatable :: xface(:),yface(:),zface(:)
real, allocatable :: rhokap(:,:,:,:)
real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)

real, allocatable :: image(:,:,:),deposit(:,:),dep(:)
real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
real, allocatable :: depositGLOBAL(:,:),depGLOBAL(:)

real, allocatable :: cdf(:),mua_array(:,:)
real, allocatable :: mus_array(:,:),fluro_array(:,:)
end MODULE iarray
