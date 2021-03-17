!> Main Program for Video 4 of VSCode for Fortran
!!
!! @author Lukas Lamm
program main
    
    implicit none
    
    real(kind=8), dimension(6) :: v

    v = (/1, 2, 3, 5, 6, 7/)
    
    ! Plot results to screen
    write(*,'(a, f10.7)') 'Self-written routine: ', norm(v)
    write(*,'(a, 3x, f10.7)') 'Intrinsic routine: ', norm2(v)
    write(*,'(a, 10x, f10.7)') 'Difference: ', norm2(v) - norm(v)

contains

    ! Calculate the euclidean norm of an n-dimensional vector
    pure function norm(vector) result(res)

        real(kind=8), dimension(:), intent(in) :: vector
        real(kind=8)                           :: res
        integer                                :: i

        res = 0
        do i = 1,size(vector),2
            res = res + vector(i)**2
        end do

        res = sqrt(res)

    end function norm

end program main