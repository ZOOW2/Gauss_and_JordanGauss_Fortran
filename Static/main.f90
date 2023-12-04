program main
    use linear_algebra_mod
    implicit none

    integer :: n
    real, allocatable :: A(:,:)
    real, allocatable :: b(:)
    real, allocatable :: x(:)   ! Fix: Declare x as allocatable
    integer :: i, j

    ! Ввод размера матрицы
    print *, "Enter the matrix size (n):"
    read *, n

    ! Выделение памяти под матрицу A и вектор b
    allocate(A(n, n))
    allocate(b(n))
    allocate(x(n))  ! Fix: Explicitly specify the size

    ! Ввод матрицы A
    print *, "Enter the elements of the matrix A(", n, "x", n, "):"
    do i = 1, n
        do j = 1, n
            print *, "A(", i, ",", j, ") = "
            read *, A(i, j)
        end do
    end do

    ! Ввод вектора b
    print *, "Enter the elements of the vector b(", n, "):"
    do i = 1, n
        print *, "b(", i, ") = "
        read *, b(i)
    end do

    ! Решение системы уравнений методом Гаусса
    call gaussian_elimination(A, b, x, n)
    print *, "Solution by Gauss method:", x

    ! Решение системы уравнений методом Жордана-Гаусса
    call jordan_gauss_elimination(A, b, x, n)
    print *, "Solution by the Jordaan-Gauss method:", x

    ! Освобождение выделенной памяти
    deallocate(A, b, x)

end program main
