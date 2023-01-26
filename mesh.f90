program maillage

    implicit none

    doubleprecision, dimension (:,:), allocatable :: X,Y  !les matrices des les cordonnes
    doubleprecision, dimension (:), allocatable :: Sav,Scord,Sap   !Sav pour la pndulation en avale d'aube Scord pour la pondulation dans l'aube  Sap pour la pondulation en aval d'aube
    integer :: sxav,sxcord,sxac,sx,sxap,sy,i,j
    real :: L

    print*, "Enter the cord length: "
    read*, L
    print*, "Enter the cell number of cells in x axes infront of the airfoil:"
    read*, sxav
    print*, "Enter the cell number of cells in x axes inside the airfoil:"
    read*, sxcord
    print*, "Enter the cell number of cells in x axes behind the airfoil:"
    read*, sxap
    print*, "Enter the number of cells along the y axes:"
    read*, sy

    sxac= sxav+sxcord
    sx= sxav+sxcord+sxap

    allocate(X(0:sx,0:sy))
    allocate(Y(0:sx,0:sy))

    allocate(Sav(sxav+1))
    allocate(Scord(sxcord/2+1))
    allocate(Sap(Sxap+1))

    !stretching (p,q,n,s)
    call STRETCHING(3.0,1.0,sxav+1,Sav)
    call STRETCHING(3.0,1.0,sxcord/2+1,Scord)
    call STRETCHING(3.0,1.0,sxap+1,Sap)

    j=0
!*---------------------------------la distébution des x
    do i=0,sxav
       X(i,j)=Sav(i+1)*L
    end do
    do i=sxav+sxcord/2,sxav,-1
        X(2*sxav+sxcord/2-i,j)=3.*L/2.-Scord(i-sxav+1)*L/2.
    end do
    do i=sxav+sxcord/2,sxac
        X(i,j)=3.*L/2.+Scord(i-sxav-sxcord/2+1)*L/2.
    end do
    do i=sx,sxac,-1
        X(2*sx-sxap-i,j)=L*5-Sap(i-sxac+1)*L*3.
    end do

!*----------------------------------------- la premiere ligne
    do i=0,sx
        if (i.gt.(sxav).and.i.lt.(sxac)) then
                Y(i,j)= 5*0.12*L*(0.2969*sqrt(X(i,j)/L-1)-0.1260*(X(i,j)/L-1)&
                        &-0.3516*(X(i,j)/L-1)**2+0.2843*(X(i,j)/L-1)**3-0.1015*(X(i,j)/L-1)**4)
                else
                    Y(i,j)=0
        end if
    end do
!*---------------------------------------------les autres lignes
    do j=1,sy
        do i=0,sx
            X(i,j)= X(i,0)
            Y(i,j)=Y(i,j-1)+(L-Y(i,j-1))/real(sy-j+1)
        end do
    end do
!*---------------------------------------------------
!*ramplisage des donnes
!*----------------------------------------------------
    open(unit=10,file="Data.dat")
        do j=0,sy
            do i=0,sx
                write(10,*) X(i,j),Y(i,j)
            end do
            write(10,*)
        end do
        do i=0,sx
            do j=0,sy
                write(10,*) X(i,j),Y(i,j)
            end do
            write(10,*)
        end do

        do j=0,sy
            do i=0,sx
                write(10,*) X(i,j),-Y(i,j)
            end do
            write(10,*)
        end do
        do i=0,sx
            do j=0,sy
                write(10,*) X(i,j),-Y(i,j)
            end do
            write(10,*)
        end do
    close(10)
!*-----------------------------------------------------------------
    contains

    subroutine STRETCHING (p,q,n,S)
!*----------------------------------------------------------------------
!* COMPUTE ONE-DIMENSIONAL STRETCHING FUNCTION,
!*  S = P*ETA + (1.-P)*(1.-TANH(Q*(1.-ETA))/TANH(Q)) ,
!* FOR GIVEN PARAMETERS, P and Q.
!*----------------------------------------------------------------------
       integer :: i,an,al,n
       real :: deta,tqi,eta,dum,p,q
       doubleprecision, DIMENSION (n) :: S
!
!        p=0.3
!        q=1.0
!        p=3.0
!        q=1.0

       AN = N-1
       DETA = 1./AN
       TQI = 1./TANH(Q)
       s(1)=0
       DO i = 2, n
            AL = i-1
            ETA = AL*DETA
            DUM = Q*(1.-ETA)
            DUM = 1.-TANH(DUM)*TQI
            S(i) = P*ETA+(1.-P)*DUM
        enddo

    end subroutine STRETCHING

end program
