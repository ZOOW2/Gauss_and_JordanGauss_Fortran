! linear_algebra_lib.f90

module linear_algebra_mod
    implicit none

contains

    ! Метод Гаусса для решения системы линейных уравнений
    subroutine gaussian_elimination(A, b, x, n)
        real, intent(inout) :: A(:,:)
        real, intent(inout) :: b(:)
        real, intent(out)   :: x(:)
        integer, intent(in) :: n
        real :: factor
        integer :: i, j, k

        ! Проход по каждому столбцу
        do k = 1, n-1
            ! Приведение к верхнетреугольному виду
            do i = k+1, n
                factor = A(i, k) / A(k, k)
                A(i, k+1:n) = A(i, k+1:n) - factor * A(k, k+1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do

        ! Обратная подстановка для нахождения решения
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:n) * x(i+1:n))) / A(i, i)
        end do

    end subroutine gaussian_elimination

    ! Метод Жордана-Гаусса для решения системы линейных уравнений
    subroutine jordan_gauss_elimination(A, b, x, n)
        real, intent(inout) :: A(:,:)
        real, intent(inout) :: b(:)
        real, intent(out)   :: x(:)
        integer, intent(in) :: n
        real :: factor
        integer :: i, j, k

        ! Прямой ход
        do k = 1, n
            ! Нормализация текущей строки
            factor = 1.0 / A(k, k)
            A(k, :) = A(k, :) * factor
            b(k) = b(k) * factor

            ! Вычитание текущей строки из всех остальных строк
            do i = 1, n
                if (i /= k) then
                    factor = A(i, k)
                    A(i, :) = A(i, :) - factor * A(k, :)
                    b(i) = b(i) - factor * b(k)
                end if
            end do
        end do

        ! Присвоение решения
        x = b

    end subroutine jordan_gauss_elimination

end module linear_algebra_mod
