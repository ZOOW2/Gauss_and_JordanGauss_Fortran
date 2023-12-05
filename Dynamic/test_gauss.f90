program main
  use gaussian_module
  implicit none
 
  integer, parameter :: n = 3
  real(8) :: A(n, n), b(n), x_gauss(n), x_jordan(n)
  real(8) :: C(n, n), D(n, n)
  integer :: i, j
  
 
  ! Инициализация матрицы A и вектора b согласно вашим условиям
  do i = 1, n
    do j = 1, n
      if (i == j) then
        A(i, j) = 2 * n
      else
        A(i, j) = 1
      end if
    end do
    b(i) = (n * (n + 1) / 2) + i * (2 * n - 1)
  end do
  
  C = A
  D = A
 
  ! Вызов метода Гаусса
  x_gauss = b  ! Сохраняем оригинальный вектор b
  call gaussian_elimination(C, x_gauss, n)
 
  ! Вывод результатов метода Гаусса
  print *, "Gaussian Elimination:"
  do i = 1, n
    print *, "x(", i, ") = ", x_gauss(i)
  end do
 
  ! Восстанавливаем оригинальный вектор b для метода Жордана-Гаусса
  x_jordan = b
 
  ! Вызов метода Жордана-Гаусса
  call jordan_gauss_elimination(D, x_jordan, n)
 
  ! Вывод результатов метода Жордана-Гаусса
  print *, "Jordan-Gauss Elimination:"
  do i = 1, n
    print *, "x(", i, ") = ", x_jordan(i)
  end do
 
end program main
