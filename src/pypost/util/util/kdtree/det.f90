module det
implicit none
contains
!! -----------------------
subroutine dsydet(A, IPIV, det)
real(kind=8), dimension(:,:), intent(in) :: A
integer, dimension(:), intent(in) :: IPIV
real(kind=8), intent(out) :: det
integer :: n, i

n = size(IPIV, 1)

det = 1.d0
do i = 1, n
    det = det * A(i,i)
    if(IPIV(i) /= i) det = -det
enddo


end subroutine

end module
