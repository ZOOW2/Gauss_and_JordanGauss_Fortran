program main
  use gaussian_module
  implicit none
 
  integer, parameter :: n = 3
  real(8) :: A(n, n), b(n), x_gauss(n), x_jordan(n)
  real(8) :: C(n, n), D(n, n)
  integer :: i, j
  
 
  ! ������������� ������� A � ������� b �������� ����� ��������
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
 
  ! ����� ������ ������
  x_gauss = b  ! ��������� ������������ ������ b
  call gaussian_elimination(C, x_gauss, n)
 
  ! ����� ����������� ������ ������
  print *, "Gaussian Elimination:"
  do i = 1, n
    print *, "x(", i, ") = ", x_gauss(i)
  end do
 
  ! ��������������� ������������ ������ b ��� ������ �������-������
  x_jordan = b
 
  ! ����� ������ �������-������
  call jordan_gauss_elimination(D, x_jordan, n)
 
  ! ����� ����������� ������ �������-������
  print *, "Jordan-Gauss Elimination:"
  do i = 1, n
    print *, "x(", i, ") = ", x_jordan(i)
  end do
 
end program main
