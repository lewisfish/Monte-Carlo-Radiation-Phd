MODULE iarray

implicit none
save

real, allocatable :: xface(:),yface(:),zface(:)
real, allocatable :: rhokap(:,:,:,:)
real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

real, allocatable :: noise(:,:),reflc(:,:),trans(:,:)
real, allocatable :: image(:,:,:),deposit(:,:)
real, allocatable :: kappa(:),albedo(:),dep(:)
real, allocatable :: imageGLOBAL(:,:,:),transGLOBAL(:,:)
real, allocatable :: depositGLOBAL(:,:),depGLOBAL(:)

end MODULE iarray
