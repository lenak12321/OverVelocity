SUBROUTINE B_CalcMetric(NI, NJ, X, Y, CellCenter, CellVolume, IFaceCenter, IFaceVector, &
                        JFaceCenter, JFaceVector) 

  REAL X(NI,NJ), Y(NI,NJ)                                ! Координаты узлов сетки
  REAL CellCenter(0:NI,0:NJ,2), CellVolume(NI-1,NJ-1)    ! Центры ячеек и их объемы
  REAL IFaceCenter(NI,NJ-1,2), IFaceVector(NI,NJ-1,2)    ! Центры и векторы I-граней
  REAL JFaceCenter(NI-1,NJ,2), JFaceVector(NI-1,NJ,2)    ! Центры и векторы J-граней
  REAL r(2)                                              ! Вспомогательный вектор
  INTEGER I, J, NBOUND, IBOUND, IOUT, JBOUND, JOUT

  ! =============================================
  ! ЦЕНТРЫ И ВЕКТОРЫ ГРАНЕЙ
  ! =============================================

  ! I-НАПРАВЛЕНИЕ (вертикальные грани)
  DO J = 1, NJ-1
    DO I = 1, NI
      ! Вектор между узлами вдоль J-направления
      r(1) = X(I,J+1) - X(I,J)
      r(2) = Y(I,J+1) - Y(I,J)
      
      ! Вектор грани - поворот на 90 градусов
      IFaceVector(I,J,1) =  r(2)  ! x-компонента
      IFaceVector(I,J,2) = -r(1)  ! y-компонента
      
      ! Центр грани - середина между узлами
      IFaceCenter(I,J,1) = 0.5 * (X(I,J) + X(I,J+1))
      IFaceCenter(I,J,2) = 0.5 * (Y(I,J) + Y(I,J+1))
    ENDDO
  ENDDO

  ! J-НАПРАВЛЕНИЕ (горизонтальные грани)
  DO J = 1, NJ
    DO I = 1, NI-1
      ! Вектор между узлами вдоль I-направления
      r(1) = X(I+1,J) - X(I,J)
      r(2) = Y(I+1,J) - Y(I,J)
      
      ! Вектор грани - поворот на -90 градусов
      JFaceVector(I,J,1) = -r(2)  ! x-компонента
      JFaceVector(I,J,2) =  r(1)  ! y-компонента
      
      ! Центр грани - середина между узлами
      JFaceCenter(I,J,1) = 0.5 * (X(I,J) + X(I+1,J))
      JFaceCenter(I,J,2) = 0.5 * (Y(I,J) + Y(I+1,J))
    ENDDO
  ENDDO

  ! =============================================
  ! ОБЪЕМЫ ЯЧЕЕК
  ! =============================================
  DO J = 1, NJ-1
    DO I = 1, NI-1
      ! Вектор диагонали ячейки
      r(1) = X(I+1,J+1) - X(I,J)
      r(2) = Y(I+1,J+1) - Y(I,J)
      
      ! Вычисление объема как суммы площадей двух треугольников
      CellVolume(I,J) = 0.5 * DOT_PRODUCT(IFaceVector(I,J,:), r) + &
                        0.5 * DOT_PRODUCT(JFaceVector(I,J,:), r)
    ENDDO
  ENDDO

  ! =============================================
  ! ЦЕНТРЫ ЯЧЕЕК
  ! =============================================

  ! ВНУТРЕННИЕ ЯЧЕЙКИ: центр как средневзвешенное центров граней
  DO J = 1, NJ-1
    DO I = 1, NI-1
      CellCenter(I,J,:) = (IFaceCenter(I,   J, :) * Norm2(IFaceVector(I,   J, :)) + &
                           IFaceCenter(I+1, J, :) * Norm2(IFaceVector(I+1, J, :)) + &
                           JFaceCenter(I, J,   :) * Norm2(JFaceVector(I, J,   :)) + &
                           JFaceCenter(I, J+1, :) * Norm2(JFaceVector(I, J+1, :)) ) / &
                         ( Norm2(IFaceVector(I,   J, :)) + &
                           Norm2(IFaceVector(I+1, J, :)) + &
                           Norm2(JFaceVector(I, J,   :)) + &
                           Norm2(JFaceVector(I, J+1, :)) )
    ENDDO
  ENDDO

  ! ГРАНИЧНЫЕ ЯЧЕЙКИ: центр совпадает с центром грани
  ! I-ГРАНИЦЫ (левая и правая)
  DO NBOUND = 1, 2
    IF (NBOUND .EQ. 1) THEN
      IBOUND = 1   ! Левая граница
      IOUT   = 0   ! Внешняя ячейка слева
    ELSE 
      IBOUND = NI  ! Правая граница
      IOUT   = NI  ! Внешняя ячейка справа
    ENDIF
    
    DO J = 1, NJ-1
      CellCenter(IOUT, J, :) = IFaceCenter(IBOUND, J, :)
    ENDDO
  ENDDO

  ! J-ГРАНИЦЫ (нижняя и верхняя)
  DO NBOUND = 1, 2
    IF (NBOUND .EQ. 1) THEN
      JBOUND = 1   ! Нижняя граница
      JOUT   = 0   ! Внешняя ячейка снизу
    ELSE 
      JBOUND = NJ  ! Верхняя граница
      JOUT   = NJ  ! Внешняя ячейка сверху
    ENDIF
    
    DO I = 1, NI-1
      CellCenter(I, JOUT, :) = JFaceCenter(I, JBOUND, :) 
    ENDDO
  ENDDO

END SUBROUTINE B_CalcMetric