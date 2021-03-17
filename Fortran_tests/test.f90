PROGRAM Triangle
    IMPLICIT NONE
    REAL :: a, b, c, Area, arear, test

    test = 1.45
    
    PRINT *, 'Welcome, please enter the&
                &lengths of the 3 sides.'
    READ *, a, b, c
    arear = Area(a,b,c)
    PRINT *, 'Triangle''s area:  ', arear
    PRINT *, 'Goodbye'
END PROGRAM Triangle

FUNCTION Area(x,y,z)
    IMPLICIT NONE
    REAL :: Area
    REAL, INTENT( IN ) :: x, y, z
    REAL :: theta, height
    theta = ACOS((x**2+y**2-z**2)/(2.0*x*y))
    height = x*SIN(theta)
    Area = 0.5*y*height
END FUNCTION Area