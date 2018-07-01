module MNumRecipe
implicit none

! add interface block for f77 function
interface
    real(8) function rtbis(func,x1,x2,xacc)
        real(8),intent(in) :: x1,x2
        real(8),intent(out) :: xacc
        interface
            real(8) function func(x)
                real(8), intent(in) :: x
            end function func
        end interface
    end function rtbis
end interface

end module MNumRecipe
