! linear_algebra_lib.f90

module linear_algebra_mod
    implicit none

contains

    subroutine gaussian_elimination(A, b, x, n)
        real, intent(inout) :: A(:,:)
        real, intent(inout) :: b(:)
        real, intent(out)   :: x(:)
        integer, intent(in) :: n
        real :: factor
        integer :: i, j, k

        do k = 1, n-1
            do i = k+1, n
                factor = A(i, k) / A(k, k)
                A(i, k+1:n) = A(i, k+1:n) - factor * A(k, k+1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do

        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:n) * x(i+1:n))) / A(i, i)
        end do
        x = nint(x)

    end subroutine gaussian_elimination

    subroutine jordan_gauss_elimination(A, b, x, n)
        real, intent(inout) :: A(:,:)
        real, intent(inout) :: b(:)
        real, intent(out)   :: x(:)
        integer, intent(in) :: n
        real :: factor
        integer :: i, j, k

        do k = 1, n
            factor = 1.0 / A(k, k)
            A(k, :) = A(k, :) * factor
            b(k) = b(k) * factor

            do i = 1, n
                if (i /= k) then
                    factor = A(i, k)
                    A(i, :) = A(i, :) - factor * A(k, :)
                    b(i) = b(i) - factor * b(k)
                end if
            end do
        end do

        x = nint(b)

    end subroutine jordan_gauss_elimination

end module linear_algebra_mod
