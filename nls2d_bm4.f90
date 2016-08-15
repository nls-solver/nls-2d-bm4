module bm4_m
  implicit none

contains

  subroutine set_forward_difference_matrix(n, d, FD)
    ! MATLAB: FD=1/d*spdiags([ones(n,1),-ones(n,1),ones(n,1)],[-n+1,0,1],n,n);
    integer, intent(in) :: n
    real(8), intent(in) :: d
    real(8), intent(out) :: FD(n, n)

    integer :: i

    FD(1, 1) = -1d0
    FD(n, 1) = 1d0
    do i = 2, n
      FD(i - 1, i) = 1d0
      FD(i, i) = -1d0
    end do
    FD(:, :) = FD(:, :) / d
  end subroutine set_forward_difference_matrix

  subroutine set_second_difference_matrix(n, d, D2)
    ! MATLAB: D2=1/(d^2)*spdiags([ones(n,1),ones(n,1),-2*ones(n,1),ones(n,1),ones(n,1)],[-n+1,-1,0,1,n-1],n,n);
    integer, intent(in) :: n
    real(8), intent(in) :: d
    real(8), intent(out) :: D2(n, n)

    integer :: i

    D2(1, 1) = -2d0
    D2(2, 1) = 1d0
    D2(n, 1) = 1d0
    do i = 2, n - 1
      D2(i - 1, i) = 1d0
      D2(i, i) = -2d0
      D2(i + 1, i) = 1d0
    end do
    D2(1, n) = 1d0
    D2(n - 1, n) = 1d0
    D2(n, n) = -2d0
    D2(:, :) = D2(:, :) / (d ** 2d0)
  end subroutine set_second_difference_matrix

  subroutine set_difference_xy_matrix(nx, ny, D2x, D2y, A)
    ! MATLAB: A = -(kron(eye(ny),D2x)+kron(D2y,eye(nx)));
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: D2x(nx, nx), D2y(ny, ny)
    real(8), intent(out) :: A(nx * ny, nx * ny)

    integer :: i1, i2, i, j, iy, ix, jy

    A(:, :) = 0d0
    ! kron(eye(ny),D2x)
    do iy = 1, ny
      i1 = nx * (iy - 1) + 1
      i2 = nx * (iy - 1) + nx
      A(i1 : i2, i1 : i2) = D2x(1 : nx, 1 : nx)
    end do
    ! kron(D2y,eye(nx))
    do iy = 1, ny
      do jy = 1, ny
        i = nx * (iy - 1)
        j = nx * (jy - 1)
        do ix = 1, nx
          A(i + ix, j + ix) = A(i + ix, j + ix) + D2y(iy, jy)
        end do
      end do
    end do
    A(:, :) = -1d0 * A(:, :)
  end subroutine set_difference_xy_matrix

  subroutine set_initial(nx, ny, dx, dy, U)
    ! MATLAB: U = ones(M1,M2) + 2*repmat(cos(x'),[1,M2])+2*repmat(cos(y'),[1,M1])';
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: dx, dy
    complex(kind(0d0)), intent(out) :: U(nx, ny)

    real(8) :: x, y
    integer :: i, j

    do i = 1, nx
      do j = 1, ny
        x = dx * (i - 1)
        y = dy * (j - 1)
        U(i, j) = dcmplx(1d0 + 2d0 * cos(x) + 2d0 * cos(y), 0d0)
      end do
    end do
  end subroutine set_initial

  subroutine integrate_gen_matE_integrand(c1, c2, c3, Mmat, Emat)
    real(8), intent(in) :: c1, c2, c3, Mmat(3, 3)
    real(8), intent(out) :: Emat(3, 3)

    integer, parameter :: n = 1000
    real(8), parameter :: h = 1d0 / real(n)
    integer :: i, j, k
    real(8) :: acc

    do i = 1, 3
      do j = 1, 3
        acc = (gen_matE_integrand(c1, c2, c3, Mmat, i, j, 0d0) + &
             gen_matE_integrand(c1, c2, c3, Mmat, i, j, 1d0)) / 2d0
        do k = 1, n - 1
          acc = acc + gen_matE_integrand(c1, c2, c3, Mmat, i, j, h * k)
        end do
        Emat(i, j) = h * acc
      end do
    end do
  end subroutine integrate_gen_matE_integrand

  real(8) function gen_matE_integrand(c1, c2, c3, Mmat, i, j, s)
    ! MATLAB (Mmat by eta, a, b, c, d, e, f is pre-calculated in gen_matE):
    !l1 = @(s) s./c1.*(s-c2)./(c1-c2).*(s-c3)./(c1-c3);
    !l2 = @(s) s./c2.*(s-c1)./(c2-c1).*(s-c3)./(c2-c3);
    !l3 = @(s) s./c3.*(s-c1)./(c3-c1).*(s-c2)./(c3-c2);
    !
    !eta = 300*theta-1;
    !a = -12-12*eta;
    !b = 6 + 6*eta;
    !c = -2 - 2*eta;
    !d = -12 - 18*eta;
    !e = 3*eta;
    !f = 3-eta;
    !
    !A = @(t,s) a.*t.^3.*s.^2 + 2.*b.*t.^3.*s + c.*t.^3 + 3.*b.*t.^2.*s.^2 + d.*t.^2.*s + e.*t.^2 ...
    !    + 3.*c.*t.*s.^2 + 2.*e.*t.*s + f.*t;
    !
    !c11 = @(s) A(c1,s).*l1(s);
    !c12 = @(s) A(c1,s).*l2(s);
    !c13 = @(s) A(c1,s).*l3(s);
    !
    !c21 = @(s) A(c2,s).*l1(s);
    !c22 = @(s) A(c2,s).*l2(s);
    !c23 = @(s) A(c2,s).*l3(s);
    !
    !c31 = @(s) A(c3,s).*l1(s);
    !c32 = @(s) A(c3,s).*l2(s);
    !c33 = @(s) A(c3,s).*l3(s);
    real(8), intent(in) :: c1, c2, c3, Mmat(3, 3), s
    integer, intent(in) :: i, j

    real(8) :: lj, taus(3), ss(3), c

    if (j == 1) then
      lj = (s / c1) * (s - c2) / (c1 - c2) * (s - c3) / (c1 - c3)
    else if (j == 2) then
      lj = (s / c2) * (s - c1) / (c2 - c1) * (s - c3) / (c2 - c3)
    else if (j == 3) then
      lj = (s / c3) * (s - c1) / (c3 - c1) * (s - c2) / (c3 - c2)
    end if

    if (i == 1) then
      c = c1
    else if (i == 2) then
      c = c2
    else if (i == 3) then
      c = c3
    end if

    taus(1) = c
    taus(2) = c * c / 2d0
    taus(3) = c * c * c / 3d0
    ss(1) = 1d0
    ss(2) = s
    ss(3) = s * s
    gen_matE_integrand = dot_product(taus, matmul(Mmat, ss)) * lj
  end function gen_matE_integrand

  subroutine pow(x, k, y)
    ! MATLAB: y = x .^ k
    real(8), intent(in) :: x(:)
    integer, intent(in) :: k
    real(8), intent(out) :: y(:)

    integer :: i

    y(:) = 1d0
    do i = 1, k
      y(:) = y(:) * x(:)
    end do
  end subroutine pow

  subroutine bm4func(n, A, eps, dt, theta, p0, q0, P1, Q1, P2, Q2, P3, Q3, Phi)
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n), eps, dt, theta
    real(8), intent(in) :: p0(n), q0(n), P1(n), Q1(n), P2(n), Q2(n), P3(n), Q3(n)
    real(8), intent(out) :: Phi(n * 6)

    real(8) :: q0pow2(n), q0pow3(n), p0pow2(n), p0pow3(n)
    real(8) :: P1pow2(n), P1pow3(n), P2pow2(n), P2pow3(n), P3pow2(n), P3pow3(n)
    real(8) :: Q1pow2(n), Q1pow3(n), Q2pow2(n), Q2pow3(n), Q3pow2(n), Q3pow3(n)
    real(8) :: X1(n), X2(n), X3(n), X4(n), X5(n), X6(n)
    real(8) :: Y11(n), Y21(n), Y31(n), Y41(n), Y51(n), Y61(n)
    real(8) :: Y12(n), Y22(n), Y32(n), Y42(n), Y52(n), Y62(n)

    call pow(q0, 2, q0pow2)
    call pow(q0, 3, q0pow3)
    call pow(p0, 2, p0pow2)
    call pow(p0, 3, p0pow3)
    call pow(P1, 2, P1pow2)
    call pow(P1, 3, P1pow3)
    call pow(P2, 2, P2pow2)
    call pow(P2, 3, P2pow3)
    call pow(P3, 2, P3pow2)
    call pow(P3, 3, P3pow3)
    call pow(Q1, 2, Q1pow2)
    call pow(Q1, 3, Q1pow3)
    call pow(Q2, 2, Q2pow2)
    call pow(Q2, 3, Q2pow3)
    call pow(Q3, 2, Q3pow2)
    call pow(Q3, 3, Q3pow3)

    Y11 = (1108800d0 * q0 + 1108800d0 * Q3 - 1108800d0 * Q2 - 1108800d0 * Q1) * theta - &
         68376d0 * q0 + 12936d0 * Q3 + 16632d0 * Q2 - 182952d0 * Q1
    Y12 = (582600d0 * q0pow3 + &
         (70800d0 * Q3 - 307800d0 * Q2 + 680400d0 * Q1) * q0pow2 + &
         (70800d0 * Q3pow2 + (43200d0 * Q2 + 43200d0 * Q1) * Q3 + &
         243000d0 * Q2pow2 - 388800d0 * Q1 * Q2 + 680400d0 * Q1pow2 + 582600d0 * p0pow2 + &
         (47200d0 * P3 - 205200d0 * P2 + 453600d0 * P1) * p0 + 23600d0 * P3pow2 + &
         (14400d0 * P2 + 14400d0 * P1) * P3 + 81000d0 * P2pow2 - 129600d0 * P1 * P2 + 226800d0 * P1pow2) * q0 + &
         582600d0 * Q3pow3 + &
         (680400d0 * Q2 - 307800d0 * Q1) * Q3pow2 + &
         (680400d0 * Q2pow2 - 388800d0 * Q1 * Q2 + 243000d0 * Q1pow2 + 23600d0 * p0pow2 + &
         (47200d0 * P3 + 14400d0 * P2 + 14400d0 * P1) * p0 + 582600d0 * P3pow2 + &
         (453600d0 * P2 - 205200d0 * P1) * P3 + &
         226800d0 * P2pow2 - 129600d0 * P1 * P2 + 81000d0 * P1pow2) * Q3 - &
         729000d0 * Q2pow3 - 874800d0 * Q1 * Q2pow2 + (- 874800d0 * Q1pow2 - 102600d0 * p0pow2 + &
         (14400d0 * P3 + 162000d0 * P2 - 129600d0 * P1) * p0 + 226800d0 * P3pow2 + &
         (453600d0 * P2 - 129600d0 * P1) * P3 - 729000d0 * P2pow2 - &
         583200d0 * P1 * P2 - 291600d0 * P1pow2) * Q2 - 729000d0 * Q1pow3 + &
         (226800d0 * p0pow2 + (14400d0 * P3 - 129600d0 * P2 + 453600d0 * P1) * p0 - 102600d0 * P3pow2 + &
         (162000d0 * P1 - 129600d0 * P2) * P3 - 291600d0 * P2pow2 - &
         583200d0 * P1 * P2 - 729000d0 * P1pow2) * Q1) * theta - &
         33111d0 * q0pow3 + ( - 4716d0 * Q3 + 22761d0 * Q2 - 55728d0 * Q1) * q0pow2 + &
         (756d0 * Q3pow2 + (6156d0 * Q2 - 13284d0 * Q1) * Q3 - 6561d0 * Q2pow2 + 61236d0 * Q1 * Q2 - &
         78732d0 * Q1pow2 - 33111d0 * p0pow2 + &
         (- 3144d0 * P3 + 15174d0 * P2 - 37152d0 * P1) * p0 + 252d0 * P3pow2 + &
         (2052d0 * P2 - 4428d0 * P1) * P3 - 2187d0 * P2pow2 + 20412d0 * P1 * P2 - 26244d0 * P1pow2) * q0 + &
         9549d0 * Q3pow3 + (12960d0 * Q2 - 6723d0 * Q1) * &
         Q3pow2 + &
         (14580d0 * Q2pow2 + 2916d0 * Q1 * Q2 - 9477d0 * Q1pow2 - 1572d0 * p0pow2 + &
         (504d0 * P3 + 2052d0 * P2 - 4428d0 * P1) * p0 + 9549d0 * P3pow2 + &
         (8640d0 * P2 - 4482d0 * P1) * P3 + 4860d0 * P2pow2 + 972d0 * P1 * P2 - 3159d0 * P1pow2) * &
         Q3 - 6561d0 * Q2pow3 - 52488d0 * Q1 * Q2pow2 + &
         (52488d0 * Q1pow2 + 7587d0 * p0pow2 + &
         (2052d0 * P3 - 4374d0 * P2 + 20412d0 * P1) * p0 + 4320d0 * P3pow2 + (9720d0 * P2 + 972d0 * P1) * P3 - &
         6561d0 * P2pow2 - 34992d0 * P1 * P2 + 17496d0 * P1pow2) * Q2 - 137781d0 * Q1pow3 + &
         (- 18576d0 * p0pow2 + (- 4428d0 * P3 + 20412d0 * P2 - 52488d0 * P1) * p0 - 2241d0 * P3pow2 + &
         (972d0 * P2 - 6318d0 * P1) * P3 - 17496d0 * P2pow2 + 34992d0 * P1 * P2 - 137781d0 * P1pow2) * Q1
    X1 = (- matmul(A, Y11) - eps*Y12) / 665280d0

    Y21 = (1108800d0 * p0 + 1108800d0 * P3 - 1108800d0 * P2 - 1108800d0 * P1) * theta - &
         68376d0 * p0 + 12936d0 * P3 + 16632d0 * P2 - 182952d0 * P1
    Y22 = ((582600d0 * p0 + 23600d0 * P3 - 102600d0 * P2 + 226800d0 * P1) * q0pow2 + &
         ((47200d0 * p0 + 47200d0 * P3 + 14400d0 * P2 + 14400d0 * P1) * Q3 + &
         (- 205200d0 * p0 + 14400d0 * P3 + 162000d0 * P2 - 129600d0 * P1) * Q2 + &
         (453600d0 * p0 + 14400d0 * P3 - 129600d0 * P2 + 453600d0 * P1) * Q1) * q0 + &
         (23600d0 * p0 + 582600d0 * P3 + 226800d0 * P2 - 102600d0 * P1) * Q3pow2 + &
         ((14400d0 * p0 + 453600d0 * P3 + 453600d0 * P2 - 129600d0 * P1) * Q2 + &
         (14400d0 * p0 - 205200d0 * P3 - 129600d0 * P2 + 162000d0 * P1) * Q1) * Q3 + &
         (81000d0 * p0 + 226800d0 * P3 - 729000d0 * P2 - 291600d0 * P1) * Q2pow2 + &
         (- 129600d0 * p0 - 129600d0 * P3 - 583200d0 * P2 - 583200d0 * P1) * Q1 * Q2 + &
         (226800d0 * p0 + 81000d0 * P3 - 291600d0 * P2 - 729000d0 * P1) * Q1pow2 + &
         582600d0 * p0pow3 + (70800d0 * P3 - 307800d0 * P2 + 680400d0 * P1) * p0pow2 + &
         (70800d0 * P3pow2 + (43200d0 * P2 + 43200d0 * P1) * P3 + 243000d0 * P2pow2 - &
         388800d0 * P1 * P2 + 680400d0 * P1pow2) * p0 + 582600d0 * P3pow3 + &
         (680400d0 * P2 - 307800d0 * P1) * P3pow2 + &
         (680400d0 * P2pow2 - 388800d0 * P1 * P2 + 243000d0 * P1pow2) * P3 - 729000d0 * P2pow3 - &
         874800d0 * P1 * P2pow2 - 874800d0 * P1pow2 * P2 - 729000d0 * P1pow3) * theta + &
         (- 33111d0 * p0 - 1572d0 * P3 + 7587d0 * P2 - 18576d0 * P1) * q0pow2 + &
         ((- 3144d0 * p0 + 504d0 * P3 + 2052d0 * P2 - 4428d0 * P1) * Q3 + &
         (15174d0 * p0 + 2052d0 * P3 - 4374d0 * P2 + 20412d0 * P1) * Q2 + &
         (- 37152d0 * p0 - 4428d0 * P3 + 20412d0 * P2 - 52488d0 * P1) * Q1) * q0 + &
         (252d0 * p0 + 9549d0 * P3 + 4320d0 * P2 - 2241d0 * P1) * Q3pow2 + &
         ((2052d0 * p0 + 8640d0 * P3 + 9720d0 * P2 + 972d0 * P1) * Q2 + &
         (- 4428d0 * p0 - 4482d0 * P3 + 972d0 * P2 - 6318d0 * P1) * Q1) * Q3 + &
         (- 2187d0 * p0 + 4860d0 * P3 - 6561d0 * P2 - 17496d0 * P1) * Q2pow2 + &
         (20412d0 * p0 + 972d0 * P3 - 34992d0 * P2 + 34992d0 * P1) * Q1 * Q2 + &
         (- 26244d0 * p0 - 3159d0 * P3 + 17496d0 * P2 - 137781d0 * P1) * Q1pow2 - &
         33111d0 * p0pow3 + ( - 4716d0 * P3 + 22761d0 * P2 - 55728d0 * P1) * p0pow2 + &
         (756d0 * P3pow2 + (6156d0 * P2 - 13284d0 * P1) * P3 - 6561d0 * P2pow2 + 61236d0 * P1 * P2 - &
         78732d0 * P1pow2) * p0 + 9549d0 * P3pow3 + &
         (12960d0 * P2 - 6723d0 * P1) * P3pow2 + (14580d0 * P2pow2 + 2916d0 * P1 * P2 - 9477d0 * P1pow2) * P3 - &
         6561d0 * P2pow3 - 52488d0 * P1 * P2pow2 + 52488d0 * P1pow2 * P2 - 137781d0 * P1pow3
    X2 = (matmul(A, Y21) + eps * Y22) / 665280d0

    Y31 = (277200d0 * q0 + 277200d0 * Q3 - 277200d0 * Q2 - 277200d0 * Q1) * theta + 24024d0 * q0 + &
         3696d0 * Q3 + 16632d0 * Q2 + 66528d0 * Q1
    Y32 = (145650d0 * q0pow3 + &
         (17700d0 * Q3 - 76950d0 * Q2 + 170100d0 * Q1) * q0pow2 + &
         (17700d0 * Q3pow2 + (10800d0 * Q2 + 10800d0 * Q1) * Q3 + 60750d0 * Q2pow2 - 97200d0 * Q1 * Q2 + &
         170100d0 * Q1pow2 + 145650d0 * p0pow2 + &
         (11800d0 * P3 - 51300d0 * P2 + 113400d0 * P1) * p0 + 5900d0 * P3pow2 + &
         (3600d0 * P2 + 3600d0 * P1) * P3 + 20250d0 * P2pow2 - 32400d0 * P1 * P2 + 56700d0 * P1pow2) * q0 + &
         145650d0 * Q3pow3 + &
         (170100d0 * Q2 - 76950d0 * Q1) * Q3pow2 + (170100d0 * Q2pow2 - 97200d0 * Q1 * Q2 + 60750d0 * Q1pow2 + &
         5900d0 * p0pow2 + (11800d0 * P3 + 3600d0 * P2 + 3600d0 * P1) * p0 + 145650d0 * P3pow2 + &
         (113400d0 * P2 - 51300d0 * P1) * P3 + 56700d0 * P2pow2 - 32400d0 * P1 * P2 + 20250d0 * P1pow2) * Q3 - &
         182250d0 * Q2pow3 - 218700d0 * Q1 * Q2pow2 + &
         (- 218700d0 * Q1pow2 - 25650d0 * p0pow2 + (3600d0 * P3 + 40500d0 * P2 - 32400d0 * P1) * p0 + &
         56700d0 * P3pow2 + (113400d0 * P2 - 32400d0 * P1) * P3 - 182250d0 * P2pow2 - 145800d0 * P1 * P2 - &
         72900d0 * P1pow2) * Q2 - &
         182250d0 * Q1pow3 + (56700d0 * p0pow2 + (3600d0 * P3 - 32400d0 * P2 + 113400d0 * P1) * p0 - &
         25650d0 * P3pow2 + (40500d0 * P1 - 32400d0 * P2) * P3 - 72900d0 * P2pow2 - 145800d0 * P1 * P2 - &
         182250d0 * P1pow2) * Q1) * theta &
         + 11223d0 * q0pow3 + (1674d0 * Q3 - 7695d0 * Q2 + 19278d0 * Q1) * q0pow2 + &
         (306d0 * Q3pow2 + (4212d0 * Q1 - 648d0 * Q2) * Q3 + 3645d0 * Q2pow2 - 23328d0 * Q1 * Q2 + &
         27702d0 * Q1pow2 + 11223d0 * p0pow2 + &
         (1116d0 * P3 - 5130d0 * P2 + 12852d0 * P1) * p0 + 102d0 * P3pow2 + (1404d0 * P1 - 216d0 * P2) * P3 + &
         1215d0 * P2pow2 - 7776d0 * P1 * P2 + 9234d0 * P1pow2) * q0 + &
         558d0 * Q3pow3 + (2106d0 * Q2 - 324d0 * Q1) * Q3pow2 + &
         (4374d0 * Q2pow2 - 8748d0 * Q1 * Q2 + 4374d0 * Q1pow2 + 558d0 * p0pow2 + &
         (204d0 * P3 - 216d0 * P2 + 1404d0 * P1) * p0 + 558d0 * P3pow2 + (1404d0 * P2 - 216d0 * P1) * P3 + &
         1458d0 * P2pow2 - 2916d0 * P1 * P2 + 1458d0 * P1pow2) * Q3 + &
         19683d0 * Q2pow3 + 13122d0 * Q1 * Q2pow2 + &
         ( - 13122d0 * Q1pow2 - 2565d0 * p0pow2 + ( - 216d0 * P3 + 2430d0 * P2 - 7776d0 * P1) * p0 + &
         702d0 * P3pow2 + (2916d0 * P2 - 2916d0 * P1) * P3 + 19683d0 * P2pow2 + 8748d0 * P1 * P2 - &
         4374d0 * P1pow2) * Q2 + 52488d0 * Q1pow3 + &
         (6426d0 * p0pow2 + (1404d0 * P3 - 7776d0 * P2 + 18468d0 * P1) * p0 - &
         108d0 * P3pow2 + (2916d0 * P1 - 2916d0 * P2) * P3 + 4374d0 * P2pow2 - 8748d0 * P1 * P2 + &
         52488d0 * P1pow2) * Q1
    X3 = (matmul(A, Y31) + eps * Y32) / 166320d0

    Y41 = (277200d0 * p0 + 277200d0 * P3 - 277200d0 * P2 - 277200d0 * P1) * theta + 24024d0 * p0 + &
         3696d0 * P3 + 16632d0 * P2 + 66528d0 * P1
    Y42 = ((145650d0 * p0 + 5900d0 * P3 - 25650d0 * P2 + 56700d0 * P1) * q0pow2 + &
         ((11800d0 * p0 + 11800d0 * P3 + 3600d0 * P2 + 3600d0 * P1) * Q3 + &
         (- 51300d0 * p0 + 3600d0 * P3 + 40500d0 * P2 - 32400d0 * P1) * Q2 + &
         (113400d0 * p0 + 3600d0 * P3 - 32400d0 * P2 + 113400d0 * P1) * Q1) * q0 &
         + (5900d0 * p0 + 145650d0 * P3 + 56700d0 * P2 - 25650d0 * P1) * Q3pow2 + &
         ((3600d0 * p0 + 113400d0 * P3 + 113400d0 * P2 - 32400d0 * P1) * Q2 + &
         (3600d0 * p0 - 51300d0 * P3 - 32400d0 * P2 + 40500d0 * P1) * Q1) * &
         Q3 + (20250d0 * p0 + 56700d0 * P3 - 182250d0 * P2 - 72900d0 * P1) * Q2pow2 + &
         (- 32400d0 * p0 - 32400d0 * P3 - 145800d0 * P2 - 145800d0 * P1) * Q1 * Q2 + &
         (56700d0 * p0 + 20250d0 * P3 - &
         72900d0 * P2 - 182250d0 * P1) * Q1pow2 + 145650d0 * p0pow3 + &
         (17700d0 * P3 - 76950d0 * P2 + 170100d0 * P1) * p0pow2 + &
         (17700d0 * P3pow2 + (10800d0 * P2 + 10800d0 * P1) * P3 + 60750d0 * P2pow2 - &
         97200d0 * P1 * P2 + 170100d0 * P1pow2) * p0 + 145650d0 * P3pow3 + &
         (170100d0 * P2 - 76950d0 * P1) * P3pow2 + &
         (170100d0 * P2pow2 - 97200d0 * P1 * P2 + 60750d0 * P1pow2) * P3 - 182250d0 * P2pow3 - &
         218700d0 * P1 * P2pow2 - 218700d0 * P1pow2 * P2 - 182250d0 * P1pow3) * theta + &
         (11223d0 * p0 + 558d0 * P3 - 2565d0 * P2 + 6426d0 * P1) * &
         q0pow2 + ((1116d0 * p0 + 204d0 * P3 - 216d0 * P2 + 1404d0 * P1) * Q3 + &
         (- 5130d0 * p0 - 216d0 * P3 + 2430d0 * P2 - 7776d0 * P1) * Q2 + &
         (12852d0 * p0 + 1404d0 * P3 - 7776d0 * P2 + 18468d0 * P1) * Q1) * q0 + &
         (102d0 * p0 + 558d0 * P3 + 702d0 * P2 - 108d0 * P1) * Q3pow2 + &
         ((- 216d0 * p0 + 1404d0 * P3 + 2916d0 * P2 - 2916d0 * P1) * Q2 + &
         (1404d0 * p0 - 216d0 * P3 - 2916d0 * P2 + 2916d0 * P1) * Q1) * Q3 + &
         (1215d0 * p0 + 1458d0 * P3 + 19683d0 * P2 + 4374d0 * P1) * Q2pow2 + &
         (- 7776d0 * p0 - 2916d0 * P3 + 8748d0 * P2 - 8748d0 * P1) * Q1 * Q2 + &
         (9234d0 * p0 + 1458d0 * P3 - 4374d0 * P2 + 52488d0 * P1) * Q1pow2 + 11223d0 * &
         p0pow3 + (1674d0 * P3 - 7695d0 * P2 + 19278d0 * P1) * p0pow2 + &
         (306d0 * P3pow2 + (4212d0 * P1 - 648d0 * P2) * P3 + 3645d0 * P2pow2 - 23328d0 * P1 * P2 + &
         27702d0 * P1pow2) * p0 + 558d0 * P3pow3 + (2106d0 * P2 - 324d0 * P1) * P3pow2 + &
         (4374d0 * P2pow2 - 8748d0 * P1 * P2 + 4374d0 * P1pow2) * P3 + 19683d0 * P2pow3 + &
         13122d0 * P1 * P2pow2 - 13122d0 * P1pow2 * P2 + 52488d0 * P1pow3
    X4 = (-matmul(A, Y41) - eps * Y42) / 166320d0

    Y51 = 840d0 * q0 + 840d0 * Q3 + 2520d0 * Q2 + 2520d0 * Q1
    Y52 = 357d0 * q0pow3 + (60d0 * Q3 - 243d0 * Q2 + 648d0 * Q1) * q0pow2 + &
         (60d0 * Q3pow2 + (108d0 * Q2 + 108d0 * Q1) * Q3 + 243d0 * Q2pow2 - 972d0 * Q1 * Q2 + &
         972d0 * Q1pow2 + 357d0 * p0pow2 + (40d0 * P3 - 162d0 * P2 + 432d0 * P1) * p0 + 20d0 * P3pow2 + &
         (36d0 * P2 + 36d0 * P1) * P3 + 81d0 * P2pow2 - 324d0 * P1 * P2 + 324d0 * P1pow2) * q0 + &
         357d0 * Q3pow3 + (648d0 * Q2 - 243d0 * Q1) * Q3pow2 + &
         (972d0 * Q2pow2 - 972d0 * Q1 * Q2 + 243d0 * Q1pow2 + 20d0 * p0pow2 + &
         (40d0 * P3 + 36d0 * P2 + 36d0 * P1) * p0 + 357d0 * P3pow2 + &
         (432d0 * P2 - 162d0 * P1) * P3 + 324d0 * P2pow2 - 324d0 * P1 * P2 + 81d0 * P1pow2) * Q3 + &
         2187d0 * Q2pow3 + &
         (- 81d0 * p0pow2 + (36d0 * P3 + 162d0 * P2 - 324d0 * P1) * p0 + 216d0 * P3pow2 + &
         (648d0 * P2 - 324d0 * P1) * P3 + 2187d0 * P2pow2) * Q2 + 2187d0 * Q1pow3 + &
         (216d0 * p0pow2 + (36d0 * P3 - 324d0 * P2 + 648d0 * P1) * p0 - 81d0 * P3pow2 + &
         (162d0 * P1 - 324d0 * P2) * P3 + 2187d0 * P1pow2) * Q1
    X5 =  (matmul(A, Y51) + eps * Y52) / 6720d0

    Y61 = 840d0 * p0 + 840d0 * P3 + 2520d0 * P2 + 2520d0 * P1
    Y62 = (357d0 * p0 + 20d0 * P3 - 81d0 * P2 + 216d0 * P1) * q0pow2 + &
         ((40d0 * p0 + 40d0 * P3 + 36d0 * P2 + 36d0 * P1) * Q3 + &
         (- 162d0 * p0 + 36d0 * P3 + 162d0 * P2 - 324d0 * P1) * Q2 + &
         (432d0 * p0 + 36d0 * P3 - 324d0 * P2 + 648d0 * P1) * Q1) * q0 + &
         (20d0 * p0 + 357d0 * P3 + 216d0 * P2 - 81d0 * P1) * &
         Q3pow2 + ((36d0 * p0 + 432d0 * P3 + 648d0 * P2 - 324d0 * P1) * Q2 + &
         (36d0 * p0 - 162d0 * P3 - 324d0 * P2 + 162d0 * P1) * Q1) * Q3 + &
         (81d0 * p0 + 324d0 * P3 + 2187d0 * P2) * Q2pow2 + ( - 324d0 * p0 - 324d0 * P3) * Q1 * Q2 + &
         (324d0 * p0 + 81d0 * P3 + 2187d0 * P1) * Q1pow2 + 357d0 * p0pow3 + &
         (60d0 * P3 - 243d0 * P2 + 648d0 * P1) * p0pow2 + &
         (60d0 * P3pow2 + (108d0 * P2 + 108d0 * P1) * P3 + &
         243d0 * P2pow2 - 972d0 * P1 * P2 + 972d0 * P1pow2) * p0 + 357d0 * P3pow3 + &
         (648d0 * P2 - 243d0 * P1) * P3pow2 + (972d0 * P2pow2 - 972d0 * P1 * P2 + 243d0 * P1pow2) * P3 + &
         2187d0 * P2pow3 + 2187d0 * P1pow3
    X6 = (- matmul(A, Y61) - eps * Y62) / 6720d0

    Phi(1 : n) = P1 - p0 - dt * X1
    Phi(n + 1 : n * 2) = Q1 - q0 - dt * X2
    Phi(n * 2 + 1 : n * 3) = P2 - p0 - dt * X3
    Phi(n * 3 + 1 : n * 4) = Q2 - q0 - dt * X4
    Phi(n * 4 + 1 : n * 5) = P3 - p0 - dt * X5
    Phi(n * 5 + 1 : n * 6) = Q3 - q0 - dt * X6
  end subroutine bm4func

  subroutine kron_with_eye(m, X, n, K)
    ! MATLAB: K = kron(X, eye(n))
    integer, intent(in) :: m, n
    real(8), intent(in) :: X(m, m)
    real(8), intent(out) :: K(m * n, m * n)

    integer :: i, j, l, i1, j1

    do i = 1, m
      do j = 1, m
        i1 = n * (i - 1)
        j1 = n * (j - 1)
        do l = 1, n
          K(i1 + l, j1 + l) = X(i, j)
        end do
      end do
    end do
  end subroutine kron_with_eye

  subroutine reshape_to_vector(nx, ny, U, P, Q)
    ! MATLAB:
    ! P=reshape(real(U),[nx*ny,1]); %real(U)
    ! Q=reshape(imag(U),[nx*ny,1]); %imag(U)
    integer, intent(in) :: nx, ny
    complex(kind(0d0)), intent(in) :: U(nx, ny)
    real(8), intent(out) :: P(nx * ny), Q(nx * ny)

    integer :: iy

    do iy = 1, ny
      P(nx * (iy - 1) + 1 : nx * iy) = real(U(1 : nx, iy))
      Q(nx * (iy - 1) + 1 : nx * iy) = imag(U(1 : nx, iy))
    end do
  end subroutine reshape_to_vector

  subroutine gen_matE(theta, c1, c2, c3, quadtol, Emat, Tmat, Tinv, lambdas)
    ! quadtol is unused in this code.
    real(8), intent(in) :: theta, c1, c2, c3, quadtol
    real(8), intent(out) :: Emat(3, 3), Tmat(3, 3), Tinv(3, 3), lambdas(3)

    integer :: info, i, ipiv(3)
    real(8) :: Mmat(3, 3), a, a6, a36, lambdas_imag(3), work(12), mat_work(3, 3)

    a = -300d0 * theta
    a6 = a * 6d0
    a36 = a * 36d0
    Mmat(1, 1) = a + 4d0
    Mmat(2, 1) = -a6 - 6d0
    Mmat(3, 1) = a6
    Mmat(2, 2) = a36 + 12d0
    Mmat(3, 2) = -a36
    Mmat(3, 3) = a36
    Mmat(1, 2) = Mmat(2, 1)
    Mmat(1, 3) = Mmat(3, 1)
    Mmat(2, 3) = Mmat(3, 2)
    call integrate_gen_mate_integrand(c1, c2, c3, Mmat, Emat)

    mat_work(:, :) = Emat(:, :)
    call dgeev('N', 'V', 3, mat_work, 3, lambdas, lambdas_imag, &
         0d0, 1, Tmat, 3, work, 12, info)

    mat_work(:, :) = Tmat(:, :)
    Tinv(:, :) = 0d0
    do i = 1, 3
      Tinv(i, i) = 1d0
    end do
    call dgesv(3, 3, mat_work, 3, ipiv, Tinv, 3, info)

    print *, 'lambdas', lambdas
    print *, 'lambdas_imag (must be zero)', lambdas_imag
  end subroutine gen_matE

  subroutine set_jacobian(n, eps, p0, q0, A, Jf)
    ! MATLAB: Jf = [2*eps*diag(p0.*q0), A + eps * diag(p0.^2+3*q0.^2);
    !               -A - eps * diag(p0.^2+3*q0.^2), 2*eps*diag(p0.*q0)];
    integer, intent(in) :: n
    real(8), intent(in) :: eps, p0(n), q0(n), A(n, n)
    real(8), intent(out) :: Jf(n * 2, n * 2)

    integer :: i
    real(8) :: p0pow2(n), p0pow3(n), q0pow2(n)
    real(8) :: jfd1(n), jfd2(n)

    call pow(p0, 2, p0pow2)
    call pow(p0, 3, p0pow3)
    call pow(q0, 2, q0pow2)

    Jf(:, :) = 0d0
    Jf(1 : n, n + 1 : n * 2) = A(1 : n, 1 : n)
    Jf(n + 1 : n * 2, 1 : n) = -A(1 : n, 1 : n)

    jfd1 = 2d0 * eps * (p0 * q0)
    jfd2 = eps * (p0pow2 + 3d0 * q0pow2)
    do i = 1, n
      Jf(i, i) = jfd1(i)
      Jf(n + i, n + i) = jfd1(i)
      Jf(i, n + i) = Jf(i, n + i) + jfd2(i)
      Jf(n + i, i) = Jf(n + i, i) - jfd2(i)
    end do
  end subroutine set_jacobian

  subroutine set_eye(n, X)
    ! set identity matrix with size n.
    integer, intent(in) :: n
    real(8), intent(out) :: X(n, n)

    integer :: i

    X(:, :) = 0d0
    do i = 1, n
      X(i, i) = 1d0
    end do
  end subroutine set_eye

  subroutine reshape_to_matrix(nx, ny, P, Q, U)
    ! MATLAB: U = reshape(P + sqrt(-1)*Q,[nx,ny]);
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: P(nx * ny), Q(nx * ny)
    complex(kind(0d0)), intent(out) :: U(nx, ny)

    integer :: iy
    complex(kind(0d0)), parameter :: kImagUnit = dcmplx(0d0, 1d0)
    do iy = 1, ny
      U(1 : nx, iy) = P(nx * (iy - 1) + 1 : nx * iy) + kImagUnit * Q(nx * (iy - 1) + 1 : nx * iy)
    end do
  end subroutine reshape_to_matrix

  subroutine print_result(nx, ny, iter, t, U, norm, energy, iunit)
    integer, intent(in) :: nx, ny, iter, iunit
    real(8), intent(in) :: t, norm, energy
    complex(kind(0d0)), intent(in) :: U(nx, ny)

    integer :: ix, iy

    write(iunit, *) '##step'
    write(iunit, *) 'nx', nx
    write(iunit, *) 'ny', ny
    write(iunit, *) 'iter', iter
    write(iunit, '(A, E24.16e3)') 't', t
    write(iunit, '(A, E24.16e3)') 'norm', norm
    write(iunit, '(A, E24.16e3)') 'energy', energy
    write(iunit, '(A)') '# ix iy Re Im'
    do iy = 1, ny
      do ix = 1, nx
        write(iunit, '(2I6, 2E25.16e3)') ix, iy, real(U(ix, iy)), imag(U(ix, iy))
      end do
    end do
  end subroutine print_result

  subroutine bm4(dx, dy, n, nx, ny, FDx, FDy, A, eps, TolX, quadtol, dt, t_end, Uinit, norms, energies, iunit)
    ! n == nx * ny
    integer, intent(in) :: nx, ny, iunit
    real(8), intent(in) :: dx, dy, eps, TolX, quadtol, dt, t_end
    real(8), intent(in) :: FDx(nx, nx), FDy(ny, ny), A(n, n)
    complex(kind(0d0)), intent(in) :: Uinit(nx, ny)
    real(8), intent(out) :: norms(:), energies(:)

    integer :: n, iter, simple_newton_iter
    real(8), parameter :: theta = 19d0 / 8d0, c1 = 1d0 / 3d0, c2 = 2d0 / 3d0, c3 = 1d0
    real(8) :: Emat(3, 3), Tmat(3, 3), Tinv(3, 3), lambdas(3), t
    real(8) :: kronTinv(n * 6, n * 6), kronT(n * 6, n * 6)

    complex(kind(0d0)) :: U(nx, ny)
    real(8) :: P(n), Q(n), p0(n), q0(n)
    real(8) :: P1(n), Q1(n), P2(n), Q2(n), P3(n), Q3(n)
    real(8) :: Jf(n * 2, n * 2), Phi(n * 6), modBM4(n * 6), rho(n * 6)
    real(8) :: X1(n * 2, n * 2), X2(n * 2, n * 2), X3(n * 2, n * 2)
    integer :: ipiv1(n * 2), ipiv2(n * 2), ipiv3(n * 2), info1, info2, info3

    call gen_matE(theta, c1, c2, c3, quadtol, Emat, Tmat, Tinv, lambdas)
    call kron_with_eye(3, Tinv, 2 * n, kronTinv)
    call kron_with_eye(3, Tmat, 2 * n, kronT)

    iter = 1
    t = 0d0
    U(:, :) = Uinit(:, :)
    call reshape_to_vector(nx, ny, U, P, Q)
    do
      call get_norm(nx, ny, dx, dy, U, norms(iter))
      call get_energy(nx, ny, dx, dy, eps, FDx, FDy, U, energies(iter))
      print *, 'iteration ', iter, ', time ', t, ', (norm, energy) = (', norms(iter), energies(iter), ')'
      call print_result(nx, ny, iter, t, U, norms(iter), energies(iter), iunit)
      if (t >= t_end) then
        exit
      end if
      p0(:) = P(:)
      q0(:) = Q(:)
      call set_jacobian(n, eps, p0, q0, A, Jf)

      ! Initial value for iteration.
      P1(:) = p0(:)
      Q1(:) = q0(:)
      P2(:) = p0(:)
      Q2(:) = q0(:)
      P3(:) = p0(:)
      Q3(:) = q0(:)

      call set_eye(n * 2, X1)
      call set_eye(n * 2, X2)
      call set_eye(n * 2, X3)
      X1 = X1 - dt * lambdas(1) * Jf
      X2 = X2 - dt * lambdas(2) * Jf
      X3 = X3 - dt * lambdas(3) * Jf
      call dgetrf(n * 2, n * 2, X1, n * 2, ipiv1, info1)
      call dgetrf(n * 2, n * 2, X2, n * 2, ipiv2, info2)
      call dgetrf(n * 2, n * 2, X3, n * 2, ipiv3, info3)
      if (info1 /= 0 .or. info2 /= 0 .or. info3 /= 0) then
        print *, 'result of dgetrf: ', info1, info2, info3
      end if

      simple_newton_iter = 0
      rho(:) = 1d0
      do while (maxval(abs(rho)) > TolX)  ! Simplified Newton method.
        simple_newton_iter = simple_newton_iter + 1
        !print *, 'simple Newton method iter ', simple_newton_iter, ' start'
        call bm4func(n, A, eps, dt, theta, p0, q0, P1, Q1, P2, Q2, P3, Q3, Phi)
        modBM4 = matmul(kronTinv, Phi)
        call dgetrs('No', n * 2, 1, X1, n * 2, ipiv1, modBM4(1 : n * 2), n * 6, info1)
        call dgetrs('No', n * 2, 1, X2, n * 2, ipiv2, modBM4(n * 2 + 1 : n * 4), n * 6, info2)
        call dgetrs('No', n * 2, 1, X3, n * 2, ipiv3, modBM4(n * 4 + 1 : n * 6), n * 6, info3)
        if (info1 /= 0 .or. info2 /= 0 .or. info3 /= 0) then
          print *, 'result of dgetrs: ', info1, info2, info3
        end if

        rho = matmul(kronT, - modBM4)
        P1 = P1 + rho(1 : n)
        Q1 = Q1 + rho(n + 1 : n * 2)
        P2 = P2 + rho(n * 2 + 1 : n * 3)
        Q2 = Q2 + rho(n * 3 + 1 : n * 4)
        P3 = P3 + rho(n * 4 + 1 : n * 5)
        Q3 = Q3 + rho(n * 5 + 1 : n * 6)

        if (simple_newton_iter == 20) then
          print *, '[Warn] simple Newton method max iter ', simple_newton_iter, ' reached'
          exit
        end if
      end do

      P = P3
      Q = Q3
      call reshape_to_matrix(nx, ny, P, Q, U)

      iter = iter + 1
      t = dt * (iter - 1)
    end do
  end subroutine bm4

  subroutine get_norm(nx, ny, dx, dy, U, norm)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: dx, dy
    complex(kind(0d0)), intent(in) :: U(nx, ny)
    real(8), intent(out) :: norm

    real(8) :: Uabs(nx, ny)

    Uabs = abs(U)
    norm = sum(abs(U) * abs(U)) / real(nx * ny)
  end subroutine get_norm

  subroutine get_energy(nx, ny, dx, dy, eps, FDx, FDy, U, energy)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: dx, dy, eps
    real(8), intent(in) :: FDx(nx, nx), FDy(ny, ny)
    complex(kind(0d0)), intent(in) :: U(nx, ny)
    real(8), intent(out) :: energy

    real(8) :: Uabs(nx, ny)
    real(8) :: FDxUabs(nx, ny), UFDyabs(nx, ny)

    Uabs = abs(U)
    FDxUabs = abs(matmul(FDx, U))
    UFDyabs = abs(matmul(U, FDy))
    energy = sum(FDxUabs * FDxUabs + UFDyabs * UFDyabs + eps / 2d0 * Uabs ** 4d0) / real(nx * ny)
  end subroutine get_energy
end module bm4_m

program main
  use bm4_m
  implicit none
  integer :: nx, ny, nt, n, iunit
  real(8) :: lx, ly, kPi, dx, dy, dt, t_end, eps, quadtol, TolX
  real(8), allocatable :: FDx(:, :), FDy(:, :), D2x(:, :), D2y(:, :), A(:, :), energies(:), norms(:)
  complex(kind(0d0)), allocatable :: U(:, :)

  kPi = 3.141592653589793d0
  lx = 2d0 * kPi  ! Spatial region size for x-axis.
  ly = 2d0 * kPi  ! Spatial region size for y-axis.
  nx = 10  ! Mesh size for x-axis.
  ny = 10  ! Mesh size for y-axis.
  dx = lx / nx
  dy = ly / ny

  t_end = 4d0  ! Total simulation time.
  dt = 0.01d0  ! Simulation time step width.
  nt = ceiling(t_end / dt)
  eps = 0.1  ! Factor for the nonlinear term |U|^2 .* U.
  quadtol = 1d-12  ! Unused in this code.
  TolX = 1d-13  ! Convergence criterion for the simplified Newton method.

  n = nx * ny  ! Total number of mesh points.
  allocate(FDx(nx, nx), FDy(ny, ny), D2x(nx, nx), D2y(ny, ny), A(n, n), energies(nt), norms(nt), U(nx, ny))

  call set_forward_difference_matrix(nx, dx, FDx)
  call set_forward_difference_matrix(ny, dy, FDy)
  call set_second_difference_matrix(nx, dx, D2x)
  call set_second_difference_matrix(ny, dy, D2y)
  call set_difference_xy_matrix(nx, ny, D2x, D2y, A)

  call set_initial(nx, ny, dx, dy, U)
  iunit = 7
  open(iunit, file='result_bm4_fortran.txt')
  call bm4(dx, dy, n, nx, ny, FDx, FDy, A, eps, TolX, quadtol, dt, t_end, U, norms, energies, iunit)
  close(iunit)
end program main
