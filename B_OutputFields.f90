Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,U,V,T,Ro,M)
  Integer:: IO, NI, NJ, I, J
  Real,Dimension(NI,NJ):: X, Y
  Real,Dimension(0:NI,0:NJ):: P, U, V, T, Ro, M

  Write(IO,*) 'VARIABLES = "X","Y","U","V","T","Ro","M","P"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ, 'DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F20.8)') X(1:NI,1:NJ)
  Write(IO,'(100F20.8)') Y(1:NI,1:NJ)
  Write(IO,'(100F20.8)') U(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') V(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') Ro(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') M(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') P(1:NI-1,1:NJ-1)
  
  ! Отладочный вывод в консоль
  WRITE(*,*) 'OutputFields: M(1,1) = ', M(1,1)
  WRITE(*,*) 'OutputFields: Max M = ', MAXVAL(M(1:NI-1,1:NJ-1))
  WRITE(*,*) 'OutputFields: Min M = ', MINVAL(M(1:NI-1,1:NJ-1))
  

End Subroutine B_OutputFields
