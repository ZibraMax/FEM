program Triangle
    IMPLICIT NONE
    REAL :: a, b, c, arear, test

    test = 1.45

    PRINT *, 'Welcome, please enter the&
                &lengths of the 3 sides.'
    READ *, a, b, c

    arear = Area(a, b, c)

    PRINT *, 'Triangle''s area:  ', arear
    PRINT *, 'Goodbye'

contains
    function Area(x, y, z) result(res)
        IMPLICIT NONE
        REAL :: res
        REAL, INTENT(IN) :: x, y, z
        REAL :: theta, height
        theta = ACOS((x**2 + y**2 - z**2)/(2.0*x*y))
        height = x*SIN(theta)
        res = 0.5*y*height
    end function Area

end program Triangle
