Program Main

  character(*), parameter :: InputFile = 'input.txt', &
  OutputFile = 'data.plt'  ! имена входного и выходного файлов
  character MeshFile*30        ! имя файла с расчетной сеткой
  integer, parameter :: IO = 12 ! единица ввода-вывода
  real, allocatable, dimension(:,:) :: X, Y, CellVolume          ! скалярные массивы
  real, allocatable, dimension(:,:,:) :: CellCenter, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector  ! векторные массивы
  real, allocatable, dimension(:,:,:) :: GradP
  real, allocatable, dimension(:,:) :: Ro, U, V, T, P, M
  real :: Pin, Uin, Vin, Tin, gamma, cp, Rm, cfl, rtmp, dt
  integer :: NI, NJ, Nit, I, J, scheme, order, limiter_type
 
  ! === ЧТЕНИЕ ВХОДНОГО ФАЙЛА ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO, FILE = InputFile)
  READ(IO,*) MeshFile  ! чтение имени файла с расчетной сеткой
  READ(IO,*) Pin
  READ(IO,*) Uin
  READ(IO,*) Vin
  READ(IO,*) Tin 
  READ(IO,*) gamma
  READ(IO,*) cp
  READ(IO,*) Rm
  READ(IO,*) cfl
  READ(IO,*) Nit
  READ(IO,*) scheme
  READ(IO,*) order
  READ(IO,*) limiter_type
  CLOSE(IO)

  ! === ЧТЕНИЕ КОЛИЧЕСТВА УЗЛОВ (NI,NJ) ИЗ ФАЙЛА С СЕТКОЙ ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO, FILE = MeshFile)
  READ(IO,*) NI, NJ
  WRITE(*,*) 'NI, NJ = ', NI, NJ

  ! === ВЫВОД ИСХОДНЫХ ДАННЫХ ===
  WRITE(*,*) 'Pin = ', Pin
  WRITE(*,*) 'Uin = ', Uin
  WRITE(*,*) 'Vin = ', Vin
  WRITE(*,*) 'Tin = ', Tin
  WRITE(*,*) 'gamma = ', gamma, 'cp = ', cp, 'Rm = ', Rm
  WRITE(*,*) 'scheme = ', scheme
  WRITE(*,*) 'order = ', order
  WRITE(*,*) 'scheme = ', limiter_type
  WRITE(*,*) 'CFL = ', cfl
  WRITE(*,*) 'Nit = ', Nit

  ! === ВЫДЕЛЕНИЕ ПАМЯТИ ДЛЯ ВСЕХ МАССИВОВ ===
  WRITE(*,*) 'Allocate arrays'
  allocate(X(NI,NJ))               ! mesh nodes X-coordinates
  allocate(Y(NI,NJ))               ! mesh nodes Y-coordinates
  allocate(P(0:NI, 0:NJ))          ! Давление
  allocate(U(0:NI, 0:NJ))          ! X-компонента скорости
  allocate(V(0:NI, 0:NJ))          ! Y-компонента скорости
  allocate(T(0:NI, 0:NJ))          ! Температура
  allocate(Ro(0:NI, 0:NJ))         ! Плотность
  allocate(M(0:NI, 0:NJ))          ! Число Маха
  allocate(CellVolume(NI-1, NJ-1)) ! Объемы ячеек
  allocate(CellCenter(0:NI, 0:NJ, 2))     ! Центры ячеек
  allocate(IFaceCenter(NI, NJ-1, 2))      ! Центры I-граней
  allocate(IFaceVector(NI, NJ-1, 2))      ! Векторы I-граней
  allocate(JFaceCenter(NI-1, NJ, 2))      ! Центры J-граней
  allocate(JFaceVector(NI-1, NJ, 2))      ! Векторы J-граней

  ! === ЧТЕНИЕ СЕТКИ ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J), Y(I,J), rtmp, I = 1, NI), J = 1, NJ)
  CLOSE(IO)

print *, 'Xmin = ', minval(X), ' Xmax = ', maxval(X)
print *, 'Ymin = ', minval(Y), ' Ymax = ', maxval(Y)

  ! === ВЫЧИСЛЕНИЕ МЕТРИКИ ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI, NJ, X, Y, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector) 


  ! ===== ИНИЦИАЛИЗАЦИЯ ПОЛЕЙ =====
  do J = 1, NJ-1      ! внутри расчётной области
    do I = 1, NI-1
      P = Pin
      U = Uin
      V = Vin
      T = Tin
      Ro = Pin / (Rm * Tin)
    enddo
  enddo

  ! ===== ВЫЗОВ РАСЧЕТНОЙ ПРОЦЕДУРЫ =====
  WRITE(*,*) 'Start Euler solver'       
  call B_Euler(NI, NJ, P, U, V, T, Ro, Pin, Uin, Vin, Tin, gamma, cp, Rm, cfl, Nit, &
    CellVolume, IFaceVector, JFaceVector, scheme, order, limiter_type)

  ! ===== ВЫЧИСЛЕНИЕ Ro И M (ПОСЛЕ ИНИЦИАЛИЗАЦИИ ПОЛЕЙ!) =====
  WRITE(*,*) 'Calculate Ro and M'       
  Ro(:,:) = P(:,:) / (Rm * T(:,:))
  M(:,:) = sqrt((U(:,:) * U(:,:) + V(:,:) * V(:,:)) / (gamma * Rm * T(:,:)))

  ! ===== ВЫВОД ПОЛЕЙ =====
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO, FILE = OutputFile)
  Call B_OutputFields(IO, NI, NJ, X, Y, P, U, V, T, Ro, M)
  Close(IO)

  open (IO, File = 'pres_on_wall.dat')
  write (IO,*) 'x', 'p'
  do i = 1, NI
    write (IO,*) x(i,NJ-1), p(i,NJ-1)
  enddo
  close(IO)


END PROGRAM Main