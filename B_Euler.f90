subroutine B_Euler(NI, NJ, P, U, V, T, Ro, Pin, Uin, Vin, Tin, &
	gamma, cp, Rm, cfl, Nit, CellVolume, IFaceVector, JFaceVector, &
    scheme, order, limiter_type)

implicit none

integer :: I, J, NI, NJ, Nit, k, scheme, order, limiter_type

real :: Pin, Uin, Vin, Tin, Roin
real :: gamma, cp, Rm, dt

real, dimension(0:NI,0:NJ) :: P, U, V, T

real, dimension(0:NI,0:NJ) :: ro, roU, roV, roE
real, dimension(0:NI,0:NJ) :: ro1, roU1, roV1, roE1
real, dimension(0:NI,0:NJ) :: RES1, RES2, RES3, RES4

! Массивы для реконструированных значений (только для 2-го порядка)
real, dimension(0:NI,0:NJ) :: pL, pR, uL, uR, vL, vR, TL, TR

real, dimension(4) :: qR_face, qL_face
real, dimension(2) :: SF
real, dimension(4) :: FLUX
real :: a, volume, dt_local, cfl, Vnorm, Vnorm_in, Vnorm_out
real :: a_local, U_conv

real :: CellVolume(NI-1,NJ-1), IFaceVector(NI,NJ-1,2), JFaceVector(NI-1,NJ,2)

! Объявление функций
real :: TVD_LIMITER

do I = 1, NI-1
	do J = 1, NJ-1
		ro(I,J)	= P(I,J) / (Rm * T(I,J)) 									! плотность
		roU(I,J) = ro(I,J) * U(I,J)												! импульс по X
		roV(I,J) = ro(I,J) * V(I,J)												! импульс по Y
		! полная энергия, [gamma = cp / cv => cv = cp / gamma]
		roE(I,J) = ro(I,J) * ((cp / gamma) * T(I,J) + 0.5 * (U(I,J)**2.0 + V(I,J)**2.0))
	enddo
enddo

	open (12, file='Residuals.dat')

do k = 1, Nit

	RES1(:,:) = 0.0
	RES2(:,:) = 0.0	
	RES3(:,:) = 0.0
	RES4(:,:) = 0.0

    ! ==================================================================
    ! TVD-реконструкция примитивных переменных
    ! ==================================================================
    if (order == 2) then
        call B_TVD(NI, NJ, P, U, V, T, &
                               pL, pR, uL, uR, vL, vR, TL, TR, limiter_type)
    endif

    ! ================================================================
    ! I-НАПРАВЛЕНИЕ (вход/выход)
    ! ================================================================

	do J = 1, NJ-1
		do I = 2, NI-1
            ! Выбор значений на гранях в зависимости от порядка точности
            if (order == 1) then
                ! 1-й порядок: центральные значения ячеек
                qR_face(1) = P(I,J)      ! P справа
                qR_face(2) = U(I,J)      ! U справа
                qR_face(3) = V(I,J)      ! V справа
                qR_face(4) = T(I,J)      ! T справа
                
                qL_face(1) = P(I-1,J)    ! P слева
                qL_face(2) = U(I-1,J)    ! U слева
                qL_face(3) = V(I-1,J)    ! V слева
                qL_face(4) = T(I-1,J)    ! T слева
                
            else if (order == 2) then
                ! 2-й порядок: реконструированные значения
                ! ПРАВАЯ грань между ячейками I-1 и I:
                ! - для ячейки I-1: берем значение на её ПРАВОЙ грани
                ! - для ячейки I: берем значение на её ЛЕВОЙ грани
                qL_face(1) = pR(I-1,J)  ! p на правой грани ячейки i-1
                qL_face(2) = uR(I-1,J)  ! u на правой грани ячейки i-1
                qL_face(3) = vR(I-1,J)  ! v на правой грани ячейки i-1
                qL_face(4) = TR(I-1,J)  ! T на правой грани ячейки i-1
                
                qR_face(1) = pL(I,J)    ! p на левой грани ячейки i
                qR_face(2) = uL(I,J)    ! u на левой грани ячейки i
                qR_face(3) = vL(I,J)    ! v на левой грани ячейки i
                qR_face(4) = TL(I,J)    ! T на левой грани ячейки i
            endif											

			SF(:) = IFaceVector(I,J,:) !face vector
			call CALC_FLUX(SF, qR_face, qL_face, gamma, cp, Rm, FLUX, scheme)

			RES1(I,J) = RES1(I,J) - FLUX(1)
			RES2(I,J) = RES2(I,J) - FLUX(2)
			RES3(I,J) = RES3(I,J) - FLUX(3)
			RES4(I,J) = RES4(I,J) - FLUX(4)

			RES1(I-1,J) = RES1(I-1,J) + FLUX(1)
			RES2(I-1,J) = RES2(I-1,J) + FLUX(2)
			RES3(I-1,J) = RES3(I-1,J) + FLUX(3)
			RES4(I-1,J) = RES4(I-1,J) + FLUX(4)

		enddo
	enddo

    ! ================================================================
    ! J-НАПРАВЛЕНИЕ (стенка/симметрия)
    ! ================================================================
    do I = 1, NI-1
        do J = 2, NJ-1
            ! Выбор значений на гранях в зависимости от порядка точности
            if (order == 1) then
                ! 1-й порядок: центральные значения ячеек
                qR_face(1) = P(I,J)      ! P сверху
                qR_face(2) = U(I,J)      ! U сверху
                qR_face(3) = V(I,J)      ! V сверху
                qR_face(4) = T(I,J)      ! T сверху
                
                qL_face(1) = P(I,J-1)    ! P снизу
                qL_face(2) = U(I,J-1)    ! U снизу
                qL_face(3) = V(I,J-1)    ! V снизу
                qL_face(4) = T(I,J-1)    ! T снизу
                
            else if (order == 2) then
                ! 2-й порядок: реконструированные значения
                ! ВЕРХНЯЯ грань между ячейками J-1 и J:
                ! - для ячейки J-1: берем значение на её ВЕРХНЕЙ грани
                ! - для ячейки J: берем значение на её НИЖНЕЙ грани
                qL_face(1) = pR(I,J-1)  ! p на верхней грани ячейки j-1
                qL_face(2) = uR(I,J-1)  ! u на верхней грани ячейки j-1
                qL_face(3) = vR(I,J-1)  ! v на верхней грани ячейки j-1
                qL_face(4) = TR(I,J-1)  ! T на верхней грани ячейки j-1
                
                qR_face(1) = pL(I,J)    ! p на нижней грани ячейки j
                qR_face(2) = uL(I,J)    ! u на нижней грани ячейки j
                qR_face(3) = vL(I,J)    ! v на нижней грани ячейки j
                qR_face(4) = TL(I,J)    ! T на нижней грани ячейки j
            endif

            SF(:) = JFaceVector(I,J,:)
            call CALC_FLUX(SF, qR_face, qL_face, gamma, cp, Rm, FLUX, scheme)

			RES1(I,J) = RES1(I,J) - FLUX(1)
			RES2(I,J) = RES2(I,J) - FLUX(2)
			RES3(I,J) = RES3(I,J) - FLUX(3)
			RES4(I,J) = RES4(I,J) - FLUX(4)

			RES1(I,J-1) = RES1(I,J-1) + FLUX(1)
			RES2(I,J-1) = RES2(I,J-1) + FLUX(2)
			RES3(I,J-1) = RES3(I,J-1) + FLUX(3)
			RES4(I,J-1) = RES4(I,J-1) + FLUX(4)

		enddo
	enddo

    ! ================================================================
    ! ГРАНИЧНЫЕ УСЛОВИЯ
    ! ================================================================

	do J = 1, NJ-1

	! ===== supersonic inlet =====

		SF(:) = IFaceVector(1,J,:)
		!modSF = sqrt(SF(1)**2 + SF(2)**2)

        ! На входе задаются параметры набегающего потока
        Vnorm_in = (Uin * SF(1) + Vin * SF(2))

        FLUX(1) = (Pin / Rm / Tin) * Vnorm_in
        FLUX(2) = (Pin / Rm / Tin) * Vnorm_in * Uin + Pin * SF(1)
        FLUX(3) = (Pin / Rm / Tin) * Vnorm_in * Vin + Pin * SF(2)
        FLUX(4) = (Pin / Rm / Tin) * Vnorm_in * (cp * Tin + (Uin**2 + Vin**2) / 2.0)

				Res1(1,J) = Res1(1,J) - FLUX(1)
				Res2(1,J) = Res2(1,J) - FLUX(2)
				Res3(1,J) = Res3(1,J) - FLUX(3)
				Res4(1,J) = Res4(1,J) - FLUX(4)

	! ===== supersonic outlet =====

		SF(:) = IFaceVector(NI,J,:)

		Vnorm_out = (U(NI-1,J) * SF(1) + V(NI-1,J) * SF(2))
        
        FLUX(1) = ro(NI-1,J) * Vnorm_out
        FLUX(2) = ro(NI-1,J) * U(NI-1,J) * Vnorm_out + P(NI-1,J) * SF(1)
        FLUX(3) = ro(NI-1,J) * V(NI-1,J) * Vnorm_out + P(NI-1,J) * SF(2)
        FLUX(4) = ro(NI-1,J) * Vnorm_out * (Cp*T(NI-1,J) + (U(NI-1,J)**2 + V(NI-1,J)**2)/2.0)

				Res1(NI-1,J) = Res1(NI-1,J) + FLUX(1)
				Res2(NI-1,J) = Res2(NI-1,J) + FLUX(2)
				Res3(NI-1,J) = Res3(NI-1,J) + FLUX(3)
				Res4(NI-1,J) = Res4(NI-1,J) + FLUX(4)

	enddo


	do I = 1, NI-1

	! ===== down symmetry =====

		SF(:) = JFaceVector(I,1,:)

        ! Для симметрии: нормальная скорость = 0, а на нее всё домножается
        
        FLUX(1) = 0.0  ! Нет потока массы через стенку
        FLUX(2) = P(I,1) * SF(1)  ! Только давление
        FLUX(3) = P(I,1) * SF(2)  ! Только давление
        FLUX(4) = 0.0  ! Нет потока энергии через стенку

				Res1(I,1) = Res1(I,1) - FLUX(1)
				Res2(I,1) = Res2(I,1) - FLUX(2)
				Res3(I,1) = Res3(I,1) - FLUX(3)
				Res4(I,1) = Res4(I,1) - FLUX(4)

	! ===== up wall =====

		SF(:) = JFaceVector(I,NJ,:)

        ! Для стенки: вся скорость = 0 (непротекание)

        FLUX(1) = 0.0 								! Нет потока массы через стенку
        FLUX(2) = P(I,NJ-1) * SF(1)		! Только давление
        FLUX(3) = P(I,NJ-1) * SF(2)		! Только давление
        FLUX(4) = 0.0									! Нет потока энергии через стенку

				Res1(I,NJ-1) = Res1(I,NJ-1) + FLUX(1)
				Res2(I,NJ-1) = Res2(I,NJ-1) + FLUX(2)
				Res3(I,NJ-1) = Res3(I,NJ-1) + FLUX(3)
				Res4(I,NJ-1) = Res4(I,NJ-1) + FLUX(4)

	enddo

    ! ================================================================
    ! ОБНОВЛЕНИЕ КОНСЕРВАТИВНЫХ ПЕРЕМЕННЫХ
    ! ================================================================
	do I = 1, NI-1
		do J = 1, NJ-1

			dt = cfl*sqrt(CellVolume(I,J))/sqrt(gamma*Rm*Tin)

			!a_local = sqrt(gamma * Rm * T(I,J))
			!U_conv = sqrt(U(I,J)**2 + V(I,J)**2) + a_local
			!dt = cfl * CellVolume(I,J) / U_conv

            ro1(I,J) = ro(I,J) - dt * RES1(I,J) / CellVolume(I,J)
            roU1(I,J) = roU(I,J) - dt * RES2(I,J) / CellVolume(I,J)
            roV1(I,J) = roV(I,J) - dt * RES3(I,J) / CellVolume(I,J)
            roE1(I,J) = roE(I,J) - dt * RES4(I,J) / CellVolume(I,J)
		enddo
	enddo

	ro(1:NI-1, 1:NJ-1) = ro1(1:NI-1, 1:NJ-1)
	roU(1:NI-1, 1:NJ-1) = roU1(1:NI-1, 1:NJ-1)
	roV(1:NI-1, 1:NJ-1) = roV1(1:NI-1, 1:NJ-1)
	roE(1:NI-1, 1:NJ-1) = roE1(1:NI-1, 1:NJ-1)


    ! ================================================================
    ! ВЫЧИСЛЕНИЕ ПРИМИТИВНЫХ ПЕРЕМЕННЫХ
    ! ================================================================
	do I = 1, NI-1
		do J = 1, NJ-1
			U(I,J) = roU(I,J) / ro(I,J)
			V(I,J) = roV(I,J) / ro(I,J)
			P(I,J) = (gamma - 1.0) * (roE(I,J) - 0.5 * ro(I,J) * (U(I,J)**2 + V(I,J)**2))
			T(I,J) = P(I,J) / (ro(I,J) * Rm)
		enddo
	enddo

    ! ================================================================
    ! ВЫВОД НЕВЯЗОК
    ! ================================================================

if ((k==1).or.(k==10).or.(k==20).or.(k==30).or.(k==40).or.(k==50).or.&
	(k==(Nit/5)).or.(k==(Nit/5*2)).or.(k==(Nit/5*3)).or.(k==(Nit/5*4)).or.(k==Nit)) then
write (*,*)		k, maxval(Res1(1:NI-1,1:NJ-1)), maxval(Res2(1:NI-1,1:NJ-1)), &
				maxval(Res3(1:NI-1,1:NJ-1)), maxval(Res4(1:NI-1,1:NJ-1))                    !, point_U(i,j)
endif

write (12,*)	k, maxval(Res1(1:NI-1,1:NJ-1)), maxval(Res2(1:NI-1,1:NJ-1)), &
				maxval(Res3(1:NI-1,1:NJ-1)), maxval(Res4(1:NI-1,1:NJ-1))

enddo

close(12)

end subroutine B_Euler

! ==================================================================
! Подпрограмма TVD-реконструкции
! ==================================================================
subroutine B_TVD(NI, NJ, P, U, V, T, &
                             pL, pR, uL, uR, vL, vR, TL, TR, limiter_type)
implicit none
integer :: NI, NJ, limiter_type
real, dimension(0:NI,0:NJ) :: P, U, V, T
real, dimension(0:NI,0:NJ) :: pL, pR, uL, uR, vL, vR, TL, TR
integer :: i, j
real :: r, psi, delta
real :: delta_right, delta_left
real :: TVD_LIMITER

do i = 1, NI-1
    do j = 1, NJ-1
        ! ------------------------------------------------------------
        ! Реконструкция в I-направлении
        ! ------------------------------------------------------------
        
        ! Давление P
        if (i > 1 .and. i < NI-1) then
            delta_right = P(i+1,j) - P(i,j)
            delta_left  = P(i,j) - P(i-1,j)
            
            if (abs(delta_right) > 1e-10) then
                r = delta_left / delta_right
            else
                r = 1.0
            endif
            
            psi = TVD_LIMITER(r, limiter_type)
            delta = 0.5 * psi * delta_right
            
            pR(i,j) = P(i,j) + delta   ! На правой грани
            pL(i,j) = P(i,j) - delta   ! На левой грани
        else
            ! На границах - первый порядок
            pR(i,j) = P(i,j)
            pL(i,j) = P(i,j)
        endif
        
        ! Скорость U
        if (i > 1 .and. i < NI-1) then
            delta_right = U(i+1,j) - U(i,j)
            delta_left  = U(i,j) - U(i-1,j)
            
            if (abs(delta_right) > 1e-10) then
                r = delta_left / delta_right
            else
                r = 1.0
            endif
            
            psi = TVD_LIMITER(r, limiter_type)
            delta = 0.5 * psi * delta_right
            
            uR(i,j) = U(i,j) + delta
            uL(i,j) = U(i,j) - delta
        else
            uR(i,j) = U(i,j)
            uL(i,j) = U(i,j)
        endif
        
        ! Скорость V
        if (i > 1 .and. i < NI-1) then
            delta_right = V(i+1,j) - V(i,j)
            delta_left  = V(i,j) - V(i-1,j)
            
            if (abs(delta_right) > 1e-10) then
                r = delta_left / delta_right
            else
                r = 1.0
            endif
            
            psi = TVD_LIMITER(r, limiter_type)
            delta = 0.5 * psi * delta_right
            
            vR(i,j) = V(i,j) + delta
            vL(i,j) = V(i,j) - delta
        else
            vR(i,j) = V(i,j)
            vL(i,j) = V(i,j)
        endif
        
        ! Температура T
        if (i > 1 .and. i < NI-1) then
            delta_right = T(i+1,j) - T(i,j)
            delta_left  = T(i,j) - T(i-1,j)
            
            if (abs(delta_right) > 1e-10) then
                r = delta_left / delta_right
            else
                r = 1.0
            endif
            
            psi = TVD_LIMITER(r, limiter_type)
            delta = 0.5 * psi * delta_right
            
            TR(i,j) = T(i,j) + delta
            TL(i,j) = T(i,j) - delta
        else
            TR(i,j) = T(i,j)
            TL(i,j) = T(i,j)
        endif
    enddo
enddo

end subroutine B_TVD


! ==================================================================
! Функция TVD-ограничителя
! ==================================================================
real function TVD_LIMITER(r, limiter_type)
implicit none
real :: r
integer :: limiter_type

if (r <= 0.0) then
    TVD_LIMITER = 0.0
    return
endif

select case(limiter_type)
case(1)  ! van Leer
    TVD_LIMITER = 2.0 * r / (r + 1.0)
case(2)  ! Superbee
    TVD_LIMITER = max(0.0, min(2.0*r, 1.0), min(r, 2.0))
case(3)  ! minmod
    TVD_LIMITER = max(0.0, min(1.0, r))
case(4)  ! van Albada
    TVD_LIMITER = (r**2 + r) / (r**2 + 1.0)
case default
    TVD_LIMITER = 1.0  ! Без ограничителя (чисто второй порядок)
end select

end function TVD_LIMITER