module gaussian_module
  implicit none
 
contains
 
  subroutine gaussian_elimination(A, b, n)
    real(8), dimension(:,:), intent(inout) :: A
    real(8), dimension(:), intent(inout) :: b
    integer, intent(in) :: n
 
    integer :: i, k
    real(8) :: factor
 
    ! Прямой ход
    do k = 1, n-1
      do i = k+1, n
        factor = A(i, k) / A(k, k)
        A(i, k:n+1) = A(i, k:n+1) - factor * A(k, k:n+1)
        b(i) = b(i) - factor * b(k)
      end do
    end do
 
    ! Обратный ход
    do k = n, 1, -1
      b(k) = b(k) / A(k, k)
      do i = k-1, 1, -1
        b(i) = b(i) - A(i, k) * b(k)
      end do
    end do
    
    ! Округление
    do i = 1, n
      b(i) = nint(b(i))
    end do
 
  end subroutine gaussian_elimination
 
  subroutine jordan_gauss_elimination(A, b, n)
    real(8), dimension(:,:), intent(inout) :: A
    real(8), dimension(:), intent(inout) :: b
    integer, intent(in) :: n
 
    integer :: i, k
    real(8) :: factor
 
    ! Прямой ход
    do k = 1, n
        ! Нормализация текущей строки
        factor = 1.0 / A(k, k)
        A(k, :) = A(k, :) * factor
        b(k) = b(k) * factor

        ! Вычитание текущей строки из остальных строк
        do i = 1, n
          if (i /= k) then
            factor = A(i, k)
            A(i, :) = A(i, :) - factor * A(k, :)
            b(i) = b(i) - factor * b(k)
          end if
        end do
    end do
    
    ! Округление
    do i = 1, n
      b(i) = nint(b(i))
    end do
 
  end subroutine jordan_gauss_elimination
 
end module gaussian_module
