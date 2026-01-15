subroutine CALC_FLUX(SF, qR, qL, gamma, cp, Rm, FLUX, scheme)

    implicit none
    
    real, dimension(2) :: SF, nf
    real, dimension(4) :: qR, qL, FLUX, wL, wR, FL, FR
    real :: gamma, cp, Rm, modSF
    
    real :: uR, vR, pR, TR
    real :: uL, vL, pL, TL
    real :: Uf, Vf, Pf, Tf
    real :: VnL, VnR, VnF
    real :: Rof, Ef, Hf
    real :: mass_flux

    INTEGER :: scheme, i
    real :: rhoL, cL, HL, SL
    real :: rhoR, cR, HR, SR
    real :: rho_roe, u_roe, v_roe, H_roe
    real :: Vn_roe, csr
    real :: M_plus, M_minus, p_plus, p_minus, M_face, p_face, ML, MR, c_face

! ===== Расчёт примитивных переменных =====
pR = qR(1)
uR = qR(2)
vR = qR(3)
TR = qR(4)

pL = qL(1)
uL = qL(2)
vL = qL(3)
TL = qL(4)

modSF = sqrt(SF(1)**2 + SF(2)**2)
nf(1) = SF(1)/modSF
nf(2) = SF(2)/modSF

! Нормальные скорости
VnL = uL * nf(1) + vL * nf(2)   !left
VnR = uR * nf(1) + vR * nf(2)   !right
VnF = 0.5 * (VnL + VnR)         !face

! Плотности
rhoL = pL / (Rm * TL)
rhoR = pR / (Rm * TR)

! Полные энтальпии
HL = cp * TL + 0.5 * (uL**2 + vL**2)
HR = cp * TR + 0.5 * (uR**2 + vR**2)

! Скорости звука  
cL = sqrt(gamma * Rm * TL)
cR = sqrt(gamma * Rm * TR)

! ======================================================================================================================================================

SELECT CASE(scheme)

CASE(1) ! противопоточная схема [SUPER]

    if (VnF < 0.0) then

        FLUX(1) = rhoR * VnR * modSF
        FLUX(2) = (rhoR * VnR * UR + PL * nF(1)) * modSF
        FLUX(3) = (rhoR * VnR * VR + PL * nF(2)) * modSF
        FLUX(4) = rhoR * VnR * HR * modSF

    elseif (VnF > 0.0) then

        FLUX(1) = rhoL * VnL * modSF
        FLUX(2) = (rhoL * VnL * UL + PR * nF(1)) * modSF
        FLUX(3) = (rhoL * VnL * VL + PR * nF(2)) * modSF
        FLUX(4) = rhoL * VnL * HL * modSF

    endif

! ======================================================================================================================================================

CASE(2) ! схема HLL (решение задачи Римана) [SUPER]

! ===== Осреднение по Роу =====
rho_roe = sqrt(rhoL * rhoR)
Vn_roe = (sqrt(rhoL) * VnL + sqrt(rhoR) * VnR) / (sqrt(rhoL) + sqrt(rhoR))
H_roe = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR))

! Средняя скорость звука
csr = sqrt((gamma - 1.0) * (H_roe - 0.5 * (Vn_roe**2)))

! ===== Расчёт скоростей волн =====
SL = minval([VnL - cL, Vn_roe - csr])
SR = maxval([VnR + cR, Vn_roe + csr])

! ===== Векторы консервативных переменных =====
wL(1) = rhoL
wL(2) = rhoL * uL
wL(3) = rhoL * vL  
wL(4) = rhoL * (cp * TL - Rm * TL + 0.5 * (uL**2 + vL**2))

wR(1) = rhoR
wR(2) = rhoR * uR
wR(3) = rhoR * vR
wR(4) = rhoR * (cp * TR - Rm * TR + 0.5 * (uR**2 + vR**2))

! ===== Векторы потоков =====
FL(1) = rhoL * VnL
FL(2) = rhoL * uL * VnL + pL * nF(1)
FL(3) = rhoL * vL * VnL + pL * nF(2)
FL(4) = rhoL * HL * VnL

FR(1) = rhoR * VnR
FR(2) = rhoR * uR * VnR + pR * nF(1)
FR(3) = rhoR * vR * VnR + pR * nF(2)
FR(4) = rhoR * HR * VnR

! ===== Расчёт потока HLL =====
IF (SL >= 0.0) THEN
    FLUX = FL * modSF
ELSE IF (SR <= 0.0) THEN
    FLUX = FR * modSF
ELSE
    FLUX = (SR * FL - SL * FR + SL * SR * (wR - wL)) / (SR - SL) * modSF
ENDIF

! ======================================================================================================================================================

CASE(3) !AUSM

ML = VnL / cL
MR = VnR / cR

if (abs(ML) <= 1) then
    M_plus = 0.25 * (ML + 1.0)**2
    p_plus = 0.25 * pL * (ML + 1.0)**2 * (2.0 - ML)
elseif (abs(ML) > 1) then
    M_plus = 0.5 * ( ML + abs(ML) )
    p_plus = 0.5 * pL * ( ML + abs(ML) ) / ML
endif

if (abs(MR) <= 1) then
    M_minus = -0.25 * (MR - 1.0)**2
    p_minus = 0.25 * pR * (MR - 1.0)**2 * (2.0 + MR)    ! убрала минус
elseif (abs(MR) > 1) then
    M_minus = 0.5 * ( MR - abs(MR) )
    p_minus = 0.5 * pR * ( MR - abs(MR) ) / MR
endif

M_face = M_plus + M_minus
p_face = p_plus + p_minus

FL(1) = rhoL * cL
FL(2) = rhoL * uL * cL
FL(3) = rhoL * vL * cL
FL(4) = rhoL * HL * cL

FR(1) = rhoR * cR
FR(2) = rhoR * uR * cR
FR(3) = rhoR * vR * cR
FR(4) = rhoR * HR * cR

if ( M_face >= 0 ) then
    FLUX = M_face * FL
elseif ( M_face < 0 ) then
    FLUX = M_face * FR
endif

FLUX(2) = FLUX(2) + p_face * nF(1)
FLUX(3) = FLUX(3) + p_face * nF(2)

FLUX = FLUX * modSF

! ======================================================================================================================================================

CASE(4) !AUSM+

cL = sqrt(gamma * Rm * TL)
cR = sqrt(gamma * Rm * TR)

c_face = 0.5 * (cL + cR)

ML = VnL / c_face
MR = VnR / c_face

if (abs(ML) <= 1) then
    M_plus = 0.25 * (ML + 1.0)**2 + 1/8 * (ML**2 - 1.0)**2
    p_plus = 0.25 * pL * (ML + 1.0)**2 * (2.0 - ML) + 3/16 * ML * (ML + 1.0)**2
elseif (abs(ML) > 1) then
    M_plus = 0.5 * ( ML + abs(ML) )
    p_plus = 0.5 * pL * ( ML + abs(ML) ) / ML
endif

if (abs(MR) <= 1) then
    M_minus = -0.25 * (MR - 1.0)**2 - 1/8 * (MR**2 - 1.0)**2
    p_minus = 0.25 * pR * (MR - 1.0)**2 * (2.0 + MR) - 3/16 * MR * (MR - 1.0)**2
elseif (abs(MR) > 1) then
    M_minus = 0.5 * ( MR - abs(MR) )
    p_minus = 0.5 * pR * ( MR - abs(MR) ) / MR
endif

M_face = M_plus + M_minus
p_face = p_plus + p_minus

FL(1) = rhoL * cL
FL(2) = rhoL * uL * cL
FL(3) = rhoL * vL * cL
FL(4) = rhoL * HL * cL

FR(1) = rhoR * cR
FR(2) = rhoR * uR * cR
FR(3) = rhoR * vR * cR
FR(4) = rhoR * HR * cR

if ( M_face >= 0 ) then
    FLUX = M_face * FL
elseif ( M_face < 0 ) then
    FLUX = M_face * FR
endif

FLUX(2) = FLUX(2) + p_face * nF(1)
FLUX(3) = FLUX(3) + p_face * nF(2)

FLUX = FLUX * modSF

END SELECT

end subroutine CALC_FLUX