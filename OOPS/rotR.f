      PROGRAM ROT                                                               
CC    History:                                                                  
CC    The program rot was originated at the Universidad Iberoamericana and      
CC    is a modification of the program 3Dbox which in turn                      
CC    is the last heir of a series of codes written by Professor K.T. Yang      
CC    and his students over the last 15 years. The author of 3Dbox is Dr        
CC    Devadata Mukutmoni.                                                      
CC    In September 1991 the code was transferred to Solar Energy Laboratory     
CC    of the University of Mexico.                                              
CC    Comments added in Mexico are denoted by double C.                         
C                                                                               
C     version 1.4, 5/8/91, Uses relative rather than an absolute toleran        
C                          Does not provide an initial guess for the cou        
C                          constants du,dv, and dw in the momentum equat        
C                 7/13/91, this version can only be employed for uniform        
C                          grids                                                
C                                                                               
C     THIS PROGRAM SOLVES THREE DIMENSIONAL BUOYANCY DRIVEN FLOWS IN            
C     A RECTANGULAR CAVITY.                                                     
CC    THE CAVITY IS ASSUMED TO BE ROTATING ARROUND THE VECTOR OMEGA             
CC                                                                              
C     IMPORTANT PARAMETERS ARE THE FOLLOWING:                                   
C     RA :   THE RAYLEIGH NUMBER                                                
C     PR :   PRANDTL NUMBER                                                     
CC    RAR :  ROTATIONAL RAYLEIGH NUMBER                                         
CC    TA :   TAYLOR NUMBER                                                      
C     AX :   ASPECT RATIO IN THE X DIRECTION                                    
C     AZ :   ASPECT RATIO IN THE Z DIRECTION                                    
C     DT :   NON DIMENSIONAL TIME STEP                                          
C     NMAX:  MAXIMUM NUMBER OF TIME STEPS                                       
C     EPS:   MINIMUM RESIDUAL                                                   
C     IMAX:  MAXIMUM NUMBER OF ITERATIONS IN THE ITERATIVE SOLVER               
C                                                                               
C***********************************************************************        
C                                                                               
CC   THE PROGRAM SOLVES THE MASS, MOMENTUM AND ENERGY BALANCE EQUATIONS         
CC   FOR A BENARD-LIKE NATURAL CONVECTIVE FLOW INSIDE A ROTATING CAVITY         
CC   THERE ARE THREE INDEPENDENT DIRECTIONS THAT MUST BE SPECIFIED TO           
CC   UNAMBIGUOSLY DEFINE THE PROBLEM:                                          
CC   A) THE ORIENTATION OF THE WALL TEMPERATURE GRADIENT                        
CC   B) THE DIRECTION OF THE GRAVITY ACCELERATION VECTOR                        
CC   C) THE ROTATION VECTOR ORIENTATION.                                        
CC   IN THE PRESENT CODE, THE RELATIVE ORIENTATION BETWEEN THE THREE            
CC   VECTORS IS ARBITRARY AND SMALL (AND INDEED, LARGE) DEVIATIONS              
CC   FROM THE CLASSICAL ARRANGEMENTS ARE PERFECTLY POSSIBLE.                    
C    THE PROGRAM SOLVES THE EQUATIONS USING THE CONTROL( FINITE) VOLUME         
C    METHOD. THE GOVERNING NAVIER STOKES EQUATIONS , CONTINUITY EQUATION        
C    AND THE ENERGY EQUATIONS ARE DISCRETISED BY FINITE DIFFERENCES.            
C    THE INTERPOLATION FOR THE CONVECTIVE TERMS IS DONE USING THE QUICK         
C    SCHEME. THE ALGORITHM IS A VARIANT OF SIMPLE. A BOUSSINESQ FLUID IS        
C    ASSUMED. PRESSURE WORK AND VISCOUS DISSIPATION ARE BEING NEGLECTED.        
C                                                                               
C***********************************************************************        
CC
CC   Al final del programa se anexa una subrutina que calcula la trayectoria
CC   de una particula, usando Series de Taylor a primer orden.
CC   El nombre de la subrutina es TRACE.
CC**********************************************************************

      common/bl1/ dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum             
      COMMON/BL7/ NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                       
     &            NIP2,NJP2,NKP2,ITER,NNMAX                                     
CC
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/ smp(24,24,24),resorm(93),thot,tcool,                         
     &             du(24,24,24),dv(24,24,24),                                   
     &             dw(24,24,24),pp(24,24,24)                                    
      common/bl36/ ap(24,24,24),
     &             ae(24,24,24),aw(24,24,24),                                   
     &             an(24,24,24),as(24,24,24),                                   
     &             af(24,24,24),ab(24,24,24),
     &             sp(24,24,24),su(24,24,24)                                             
CC                                                                              
      COMMON/BL101/GX1U                                                         
      COMMON/BL102/RAR,OX,OY,OZ,OMEGAX0,OMEGAY0,OMAGS                           
      COMMON/BL103/TA                                                           
      COMMON/BL201/GY1U                                                         
      COMMON/BL301/GZ1U                                                         
CC                                                                              
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
      COMMON/ASP/AX,AZ,EPS                                                      
      COMMON/ERROR/IMAX                                                         
      common/count/nt,nmax
CC
      DIMENSION UU(40000),VV(40000),WW(40000)                                   
      DIMENSION VECTOR(3,24,24,24)
CC    DIMENSION VORT1(3,24,24,24)
CC    DIMENSION VORT(3,24,24,24)
CC    DIMENSION FCOR(3,24,24,24)
CC
      DATA sormax,itmax/2.00,5/                                                 
      DATA pi,roll/3.14159265359,2.0/   
CC
      CHARACTER*16 ARCH8, ARCH9, ARCH10, TIEMPO
      CHARACTER*16 ARCH32, ARCH33, ARCH37
      ARCH8 = 'TEMP.'
      ARCH9 = 'VELC.'
CC    ARCH10 = 'FCOR.'
      ARCH32 = 'JTEMP.'
      ARCH33 = 'JVELC.'
CC    ARCH37 = 'VORT.'

CC  Inicializamos VECTOR y VORT.      
      DO I=1,NIP1
      DO J=1,NJP1
      DO K=1,NKP1
       VECTOR(1,I,J,K) = 0.0
       VECTOR(2,I,J,K) = 0.0
       VECTOR(3,I,J,K) = 0.0
CC     VORT(1,I,J,K) = 0.0
CC     VORT(2,I,J,K) = 0.0
CC     VORT(3,I,J,K) = 0.0
      ENDDO
      ENDDO
      ENDDO 

CC  
CC    MNij : Monitores de la velocidad en (IMONj,JMONj,KMONj). i=U,V,W
CC   	     j = 1,2,3,4.  
CC    NUH : Numero de Nusselt en la pared inferior (TH).
CC    NUC : Numero de Nusselt en la pared superior (TC).
CC    Data : Datos de entrada.
CC    CINI : Temperatura, Velocidad y Presion finales, para continuar corridas.
CC    TINI : Tiempo final, para continuar corridas.
CC    xyz : Posiciones de una particula que inicia en (xa,ya,za).
CC    LOOPi : Proyeccion del espacio fase de velocidad en diferentes planos.
CC            i=1 en UV, i=2 en UW y i=3 en VW. Todos en (ILM,JLM,KLM).
CC    LOOP : Espacio fase de velocidad en 3D. En (ILM,JLM,KLM).
CC
      OPEN(12,FILE='MNV1',STATUS='UNKNOWN')                 
      OPEN(13,FILE='MNV2',STATUS='UNKNOWN')                 
      OPEN(34,FILE='MNU2',STATUS='UNKNOWN')
      OPEN(35,FILE='MNW2',STATUS='UNKNOWN')
      OPEN(14,FILE='MNU3',STATUS='UNKNOWN')                 
      OPEN(15,FILE='MNU4',STATUS='UNKNOWN')                 
      OPEN(16,FILE='MNV4',STATUS='UNKNOWN')                 
      OPEN(17,FILE='MNW4',STATUS='UNKNOWN')                 
      OPEN(18,FILE='NUH',STATUS='UNKNOWN')                 
      OPEN(19,FILE='NUC',STATUS='UNKNOWN')                 
      OPEN(20,FILE='Data',STATUS='UNKNOWN')
      OPEN(30,FILE='CINI',STATUS='UNKNOWN')
      OPEN(31,FILE='TINI',STATUS='UNKNOWN')
      OPEN(36,FILE='xyz',STATUS='UNKNOWN')
      OPEN(41,FILE='LOOP1',STATUS='UNKNOWN')
      OPEN(42,FILE='LOOP2',STATUS='UNKNOWN')
      OPEN(43,FILE='LOOP3',STATUS='UNKNOWN')
      OPEN(44,FILE='LOOP',STATUS='UNKNOWN')


CC                                                                              
CC    KRUN = 1: Continuing run. KRUN >< 1: New run.                 
CC    XPER = X-dimension of perturbation roll                                   
CC                                                                              
      READ(20,*) 
      READ(20,*) KRUN
      READ(20,*) 
      READ(20,*) RA                                          
      READ(20,*) 
      READ(20,*) PR                                           
      READ(20,*) 
      READ(20,*) IMAX                                           
      READ(20,*) 
      READ(20,*) XPER                                           
      READ(20,*) 
      READ(20,*) NMAX                                           
      READ(20,*) 
      READ(20,*) DTIME                                           
      READ(20,*) 
      READ(20,*) AX                                           
      READ(20,*) 
      READ(20,*) AZ                                           
      READ(20,*) 
      READ(20,*) EPS                                           
      READ(20,*) 
CC                                                                              
      READ(20,*) GRAVIX,GRAVIY,GRAVIZ                                           
CC    GRAVI ARE THE COMPONENTS OF THE GRAVITATIONAL ACCELERATION VECTOR         
CC    AT TIME ZERO                                                              
      READ(20,*) 
      READ(20,*) OMEGAX,OMEGAY,OMEGAZ                                           
CC    OMEGA ARE THE COMPONENTS OF THE ROTATION VECTOR                           
      READ(20,*) 
      READ(20,*) OMEGAX0,OMEGAY0                                                
CC    OMEGA0 ARE THE COORDINATES OF INTERSECTION OF THE ROTATION VECTOR         
CC    WITH THE Z=0 PLANE.                                                       
      READ(20,*) 
      READ(20,*) ALPHA
CC    ALPHA IS THE ANGLE BETWEEN THE ROTATION VECTOR AND
CC    THE GRAVITY ACCELERATION VECTOR AT TIME EQUALS ZERO.
      READ(20,*)
      READ(20,*) NPR
      READ(20,*)
      READ(20,*) NPR1
      READ(20,*)
      READ(20,*) NPRT
      READ(20,*)
      READ(20,*) IMON1
      READ(20,*) JMON1
      READ(20,*) KMON1
      READ(20,*) IMON2
      READ(20,*) JMON2
      READ(20,*) KMON2
      READ(20,*) IMON3
      READ(20,*) JMON3
      READ(20,*) KMON3
      READ(20,*) IMON4
      READ(20,*) JMON4
      READ(20,*) KMON4
      READ(20,*)
      READ(20,*) JMONO
      READ(20,*) 
      READ(20,*) TA
      READ(20,*) 
      READ(20,*) xa,ya,za
      READ(20,*)
      READ(20,*) ILM,JLM,KLM
      print *,'THE RAYLEIGH NUMBER IS ',ra                                       
      print *,'THE PRANDTL NUMBER IS ',pr                                        
      print *,'THE ASPECT RATIOS ARE ',ax,az                                    
      print *,'THE TIME STEP IS ',dtime                                          
      print *,'THE NUMBER OF TIME STEPS ',nmax                                       
CC                                                                              
CC    CALCULATE THE PARAMETERS REQUIRED IN THE ROTATION TERMS                   
CC    A) CALCULATE GRAVITY VECTOR MAGNITUDE                                     
      GMAG = SQRT(GRAVIX*GRAVIX+GRAVIY*GRAVIY+GRAVIZ*GRAVIZ)                    
CC    B) CALCULATE THE UNITY GRAVITY VECTOR                                     
C
      GX = GRAVIX/GMAG                                                          
      GY = GRAVIY/GMAG                                                          
      GZ = GRAVIZ/GMAG                                                          
C
CC    C) CALCULATE ROTATION VECTOR MAGNITUDE                                    
      OMAGS=OMEGAX*OMEGAX+OMEGAY*OMEGAY+OMEGAZ*OMEGAZ                           
      OMAG=SQRT(OMAGS)                                                          
CC    D) CALCULATE THE UNITY ROTATION VECTOR                                    
C
      OX=OMEGAX/OMAG                                                            
      OY=OMEGAY/OMAG                                                            
      OZ=OMEGAZ/OMAG                                                            
CC                                                                              
CC    DEFINE THE NONDIMENSIONAL PARAMETERS                                      
CC    ROTATIONAL REYNOLDS NUMBER                                                
      RER=0.0                                                                   
CC    ROTATIONAL RAYLEIGH NUMBER                                                
      RAR=0.0                                                                   
CC    TAYLOR NUMBER                                                             
CC      TA=300000000.0 
CC    TAYLOR NUMBER WITHOUT ROTATION
CC      TA=0.0
CC                                                                              
      PRINT*, 'THE GRAVITY VECTOR COMPONENTS ARE',GRAVIX,GRAVIY,GRAVIZ          
      PRINT*, 'THE GRAVITY VECTOR MAGNIITUDE IS',GMAG                           
      PRINT*, 'THE ROTATION VECTOR COMPONENTS ARE',OMEGAX,OMEGAY,OMEGAZ         
      PRINT*, 'THE INTERSECTION POINT OF THE ROTATION'                          
      PRINT*, 'VECTOR WITH THE Z=0 PLANE IS',OMEGAX0,OMEGAY0                    
                                                                                
         call grid                                                              
c                                                                               
c     FOR CONTINUING RUN, READ FROM DATASET    
c 
      if(krun .eq. 1) then
        do 1370 k=1,nkp1
        do 1370 j=1,njp1
        do 1370 i=1,nip1
         read(30,3030) tod(i,j,k),uod(i,j,k),vod(i,j,k),
     &                            wod(i,j,k),pod(i,j,k)
 1370   continue
      else
         do 220 j=1,njp1                                                        
           yy = (float(j) - 1.5) / float(nim1)                                  
           do 220 i=1,nip1                                                      
              do 220 k=1,nkp1                                                   
                 uod(i,j,k) = 0.                                                
                 u(i,j,k) = 0.                                                  
                 upd(i,j,k) = 0.                                                
                 vod(i,j,k) = 0.                                                
                 v(i,j,k) = 0.                                                  
                 vpd(i,j,k) = 0.                                                
                 w(i,j,k) = 0.                                                  
                 wpd(i,j,k) = 0.                                                
                 wod(i,j,k) = 0.                                                
                 pod(i,j,k) = 0.                                                
                 p(i,j,k) = 0.                                                  
                 ppd(i,j,k) = 0.                                                
                 du(i,j,k) = 0.                                                 
                 dv(i,j,k) = 0.                                                 
                 dw(i,j,k) = 0.                                                 
                 su(i,j,k) = 0.                                                 
                 sp(i,j,k) = 0.                                                 
                 pp(i,j,k) = 0.                                                 
                 ap(i,j,k) = 0.                                                 
                 aw(i,j,k) = 0.                                                 
                 ae(i,j,k) = 0.                                                 
                 an(i,j,k) = 0.                                                 
                 as(i,j,k) = 0.                                                 
                 smp(i,j,k) = 0.                                                
CC               tod(i,j,k) = 0.5 - yy                                          
                 tod(i,j,k) = 0.                                          
                 t(i,j,k) = tod(i,j,k)                                          
                 tpd(i,j,k) = tod(i,j,k)                                        
                 T_MEAN(I,J,K) = 0.                                             
                 U_MEAN(I,J,K) = 0.                                             
                 V_MEAN(I,J,K) = 0.                                             
                 W_MEAN(I,J,K) = 0.                                             
                 P_MEAN(I,J,K) = 0.                                             
  220         continue                                                          
c                                                                               
c        AMPLITUDE OF PERTURBATIONS                                             
c                                                                               
       A1 = 1. / AX                                                      
       ROLL = 1.0
       YPER = ROLL * XPER / AX                                           
c
c	 Defasamiento en las componentes de la velocidad
c	 de la perturbacion para obtener un rollo simetrico.
c
       DF1 = 1.5
       DF2 = 1.5
c                                                                               
c        U VELOCITY PERTURBATIONS                                               
c                                                                               
       do 230 k=2,nk                                                            
          DO 230 J=2,NJ                                                         
             DO 230 I=2,NI                                                      
                U(I,J,K) = - XPER * COS(PI * DY * (FLOAT(J) - DF1))             
     &                * SIN(ROLL * PI * DX * A1 * (FLOAT(I) - DF2))           
                UOD(I,J,K) = U(I,J,K)                                           
                UPD(I,J,K) = U(I,J,K)                                           
 230  CONTINUE                                                                  
c                                                                               
c        V VELOCITY PERTURBATIONS                                               
c                                                                              
       do 235 k=2,nk                                                            
          DO 235 J=2,NJ                                                         
             DO 235 I=2,NI                                                      
                V(I,J,K) = YPER * SIN(PI * DY * (FLOAT(J) - DF2))               
     &              * COS(ROLL * PI * DX * A1 * (FLOAT(I) - DF1))               
                VOD(I,J,K) = V(I,J,K)                                           
                VPD(I,J,K) = V(I,J,K)                                           
 235   continue                                                                 
      endif                                                                     
c                                                                               
c          SET THE HOT AND COLD WALL TEMPERATURE                                
c                                                                               
      thot = 0.5                                                                
      TCOOL = - 0.5                                                             
c                                                                               
c          SET THE VELOCITIES OUTSIDE THE COMPUTATIONAL                         
c                DOMAIN TO ZERO                                                 
c                                                                               
      do 505 i=1,nip1                                                           
         DO 505 J=1,NJP1                                                        
            UOD(I,J,1) = 0.                                                     
            UOD(I,J,NKP1) = 0.                                                  
            VOD(I,J,1) = 0.                                                     
            VOD(I,J,NKP1) = 0.                                                  
            WOD(I,J,1) = 0.                                                     
            WOD(I,J,2) = 0.                                                     
            WOD(I,J,NKP1) = 0.                                                  
  505 continue                                                                  
      do 510 i=1,nip1                                                           
         DO 510 K=1,NKP1                                                        
            UOD(I,1,K) = 0.                                                     
            UOD(I,NJP1,K) = 0.                                                  
            VOD(I,1,K) = 0.                                                     
            VOD(I,2,K) = 0.                                                     
            VOD(I,NJP1,K) = 0.                                                  
            WOD(I,1,K) = 0.                                                     
            WOD(I,NJP1,K) = 0.                                                  
  510 CONTINUE                                                                  
      do 520 j=1,njp1                                                           
         DO 520 K=1,NKP1                                                        
            UOD(1,J,K) = 0.                                                     
            UOD(2,J,K) = 0.                                                     
            UOD(NIP1,J,K) = 0.                                                  
            VOD(1,J,K) = 0.                                                     
            VOD(NIP1,J,K) = 0.                                                  
            WOD(1,J,K) = 0.                                                     
            WOD(NIP1,J,K) = 0.                                                  
  520 continue                                                                  
c                                                                               
c        INITIALISE u,v,w,t,p                                                   
c                                                                               
      do 210 k=1,nkp1                                                           
         do 210 j=1,njp1                                                        
            do 210 i=1,nip1                                                     
               t(i,j,k) = tod(i,j,k)                                            
               u(i,j,k) = uod(i,j,k)                                            
               v(i,j,k) = vod(i,j,k)                                            
               w(i,j,k) = wod(i,j,k)                                            
               p(i,j,k) = pod(i,j,k)                                            
  210 continue                                                                  
c 
CC                                                                              
      ICOUNTER =1
CC
      if(krun .eq. 1) then
	read(31,3001) nt,xtime
      else
        xtime = 0.                     
        nt = 0                                     
      endif
c                                                                               
c    THE MAIN TIME MARCHING LOOP BEGINS HERE                                    
c                                                                               
CC                                                                              
CC   ALPHA IS THE ROTATED ANGLE                                                 
CC   EXPRESS THELAR VELOCITY IN RADIANS PER SECOND                         
CC
CC    ALPHA=OMAGSR*XTIME                                                        
                                                        
CC   THE VARIABLES GX1,GY1,GZ1 ARE THE GRAVITY ACCELERATION                     
CC   COMPONENTS IN  A SYSTEM ROTATING AROUND THE VECTOR OMEGA                   
      GDOTOM=GX*OX+GY*OY+GZ*OZ                          
      GX1=GDOTOM*OX+COS(ALPHA)*(GX-GDOTOM*OX)                       
     &   -SIN(ALPHA)*(OY*GZ-OZ*GY)                              
      GY1=GDOTOM*OY+COS(ALPHA)*(GY-GDOTOM*OY)                       
     &   +SIN(ALPHA)*(OX*GZ-OZ*GX)                              
      GZ1=GDOTOM*OZ+COS(ALPHA)*(GZ-GDOTOM*OZ)                       
     &   -SIN(ALPHA)*(OX*GY-OY*GX)                             
CC                                                                              
CC    THE VARIABLES GX1U,GY1U,GZ1U ARE THE UNIT VECTORS CORRESPONDING           
CC    TO GX1,GY1 AND GZ1                                                        
CC    G1M= SQRT(GX1*GX1+GY1*GY1+GZ1*GZ1)                                        
CC    GX1U=GX1/G1M                                                              
CC    GY1U=GY1/G1M                                                 
CC    GZ1U=GZ1/G1M                                                 
CC	
CC	The next components are the corrects (31-Aug-94)
CC
      GX1U = 0
      GY1U = 1
      GZ1U = 0
CC                                                                              
  300 continue                                                                  
c                                                                               
      nt = nt + 1                                                               
CC                                                                              
      PRINT *, 'TIME STEP =',NT                                                 
CC                                                                              
c                                                                               
c    UPDATE THE MEAN VELOCITIES,TEMPERATURE AND                                 
c               PRESSURE                                                        
c                                                                               
      C1 = 1. / FLOAT(NT)                                                       
      do 310 k=1,nkp1                                                           
         DO 310 J=1,NJP1                                                        
            DO 310 I=1,NIP1                                        
               U_MEAN(I,J,K) = (1. - C1) * U_MEAN(I,J,K)                        
     &                       + c1 * u(i,j,k)                                    
               V_MEAN(I,J,K) = (1. - C1) * V_MEAN(I,J,K)                        
     &                       + c1 * v(i,j,k)                                    
               W_MEAN(I,J,K) = (1. - C1) * W_MEAN(I,J,K)                        
     &                       + c1 * w(i,j,k)                                    
               T_MEAN(I,J,K) = (1. - C1) * T_MEAN(I,J,K)                        
     &                       + c1 * t(i,j,k)                                    
  310 continue

CC  Si ya llego al maximo paso en el tiempo entonces termina
 
      if(nt.gt.nmax) then                                                       
         stop                                                                   
      endif

CC  Se escriben los resultados en los archivos TEMP.***** y VELC.*****. 
CC  El numero de archivo se calcula de acuerdo al tiempo en que se calcularon
CC  los datos.
                                                                     
      IF (MOD((NT-0),NPR) .EQ. 0) THEN
           IU = MOD(NT,10)
           ITOT = (NT - IU)/10
	   ID = MOD(ITOT,10)
	   ITOT = (ITOT - ID)/10
	   IC = MOD(ITOT,10)
           ITOT = (ITOT - IC)/10
           IM = MOD(ITOT,10) 
           ITOT = (ITOT - IM)/10
           IX = MOD(ITOT,10) 
           WRITE(TIEMPO(1:1),3000) IX
           WRITE(TIEMPO(2:2),3000) IM
	   WRITE(TIEMPO(3:3),3000) IC
	   WRITE(TIEMPO(4:4),3000) ID
	   WRITE(TIEMPO(5:5),3000) IU
 3000 FORMAT(I1)
	   ARCH8(6:10) = TIEMPO
	   ARCH9(6:10) = TIEMPO
	   ARCH10(6:10) = TIEMPO
 	   ARCH37(6:10) = TIEMPO
           OPEN(8,FILE=ARCH8,STATUS='UNKNOWN')    
           OPEN(9,FILE=ARCH9,STATUS='UNKNOWN')     
CC         OPEN(10,FILE=ARCH10,STATUS='UNKNOWN')
CC	   OPEN(37,FILE=ARCH37,STATUS='UNKNOWN')

CC  Al asignar las velocidades al arreglo VECTOR, se hace una interpolacion
CC  a los centros de los volumenes usados en la resolucion de las ecuaciones.
CC  Se usan NIP1 X NJP1 X NKP1 volumenes finitos (prismas rectangulares).
CC  Tales centros corresponden a los nodos de la malla en MPGS (que tambien
CC  es de NIP1 X NJP1 X NKP1 ).

           DO J=2,NJ
           DO I=2,NI
           DO K=2,NK
            VECTOR(1,I,J,K)=(U(I,J,K) + U(I+1,J,K))/2.
            VECTOR(2,I,J,K)=(V(I,J,K) + V(I,J+1,K))/2.
            VECTOR(3,I,J,K)=(W(I,J,K) + W(I,J,K+1))/2.

CC Se calcula la fuerza de coriolis a partir de la velocidad y del vector de
CC de rotacion.

CC	    FCOR(1,I,J,K)=OY*VECTOR(3,I,J,K)-OZ*VECTOR(2,I,J,K)
CC	    FCOR(2,I,J,K)=OZ*VECTOR(1,I,J,K)-OX*VECTOR(3,I,J,K)
CC	    FCOR(3,I,J,K)=OX*VECTOR(2,I,J,K)-OY*VECTOR(1,I,J,K)
           ENDDO
           ENDDO
           ENDDO
	
CC
CC Interpolacion de la temperatura
CC
        DO J=1,NJP1
        DO K=1,NKP1
        DO I=1,NIP1
        if((J.NE.1).AND.(J.NE.NJP1)) then
         if((I.EQ.1).AND.(K.EQ.1))       T(I,J,K) = (T(I,J,K) +
     &                                               T(2,J,2))/2
         if((I.EQ.NIP1).AND.(K.EQ.NKP1)) T(I,J,K) = (T(I,J,K) +
     &                                               T(NI,J,NK))/2
         if((I.EQ.1).AND.(K.EQ.NKP1)) T(I,J,K) = (T(I,J,K) +
     &                                           T(2,J,NK))/2
         if((I.EQ.NIP1).AND.(K.EQ.1)) T(I,J,K) = (T(I,J,K) +
     &                                           T(NI,J,2))/2
        endif
        if(J.EQ.1) then
          T(I,1,K)= (T(I,1,K) + T(I,2,K))/2
          if((I.EQ.1).OR.(I.EQ.NIP1)) T(I,J,K) = 0.5
          if((K.EQ.1).OR.(K.EQ.NKP1)) T(I,J,K) = 0.5
        endif
        If(J.EQ.NJP1) then
          T(I,NJP1,K)= (T(I,NJP1,K) + T(I,NJ,K))/2
          if((I.EQ.NIP1).OR.(I.EQ.1)) T(I,J,K) = -0.5
          if((K.EQ.NKP1).OR.(K.EQ.1)) T(I,J,K) = -0.5
        endif
        ENDDO
        ENDDO
        ENDDO

CC
CC Se escriben los resultados con el formato MPGS.
CC
           WRITE(8,3010) "TFIELD"                                 
           WRITE(8,3030) (((T(I,J,K),I=1,NIP1),J=1,NJP1),K=1,NKP1)
           WRITE(9,3010) "VECTORFIELD"                                 
           WRITE(9,3030) ((((VECTOR(IJK,I,J,K),IJK=1,3),I=1,NIP1)
     &                   ,J=1,NJP1),K=1,NKP1)
CC	   WRITE(10,3010) "CORIOLIS FORCE"
CC	   WRITE(10,3030) ((((FCOR(IJK,I,J,K),IJK=1,3),I=1,NIP1)
CC   &                   ,J=1,NJP1),K=1,NKP1)

CC
CC Se calcula la vorticidad.
CC
CC         DO I=1,NI
CC         DO J=1,NJ
CC         DO K=1,NK
CC          DWDY=(W(I+1,J+1,K+1)-W(I+1,J,K+1))/DY
CC          DVDZ=(V(I+1,J+1,K+1)-V(I+1,J+1,K))/DZ 
CC          DUDZ=(U(I+1,J+1,K+1)-U(I+1,J+1,K))/DZ 
CC          DWDX=(W(I+1,J+1,K+1)-W(I,J+1,K+1))/DX
CC          DUDY=(U(I+1,J+1,K+1)-U(I+1,J,K+1))/DY
CC          DVDX=(V(I+1,J+1,K+1)-V(I,J+1,K+1))/DX
CC          VORT1(1,I,J,K)=  DWDY-DVDZ
CC          VORT1(2,I,J,K)=-(DWDX-DUDZ)
CC          VORT1(3,I,J,K)=  DVDX-DUDY
CC         ENDDO
CC         ENDDO
CC         ENDDO

CC         DO I=2,NI
CC         DO J=2,NJ
CC         DO K=2,NK
CC          VORT(1,I,J,K)=(VORT1(1,I,J,K)+VORT1(1,I,J,K-1)+
CC   &                    VORT1(1,I,J-1,K-1)+VORT1(1,I,J-1,K))/4.
CC          VORT(2,I,J,K)=(VORT1(2,I,J,K)+VORT1(2,I,J,K-1)+
CC   &                    VORT1(2,I-1,J,K-1)+VORT1(2,I-1,J,K))/4.
CC          VORT(3,I,J,K)=(VORT1(3,I,J,K)+VORT1(3,I-1,J,K)+
CC   &                    VORT1(3,I-1,J-1,K)+VORT1(3,I,J-1,K))/4.
CC
CC         ENDDO
CC         ENDDO
CC         ENDDO

CC	   WRITE(37,3010) "VORTICITYFIELD"
CC         WRITE(37,3030) ((((VORT(IJK,I,J,K),IJK=1,3),I=1,NIP1)
CC   &                   ,J=1,NJP1),K=1,NKP1)

	   ARCH32(7:11) = TIEMPO
	   ARCH33(7:11) = TIEMPO
	   OPEN(32,FILE=ARCH32,STATUS='UNKNOWN')
	   OPEN(33,FILE=ARCH33,STATUS='UNKNOWN')
	  
CC Se escriben los resultados en un corte paralelo al plano horizontal 
CC a una altura de JMONO.
 
	   do 1372 k=1,nkp1
            do 1372 i=1,nip1
             x = AX*float(i-1)/(float(ni))
             z = AZ*float(k-1)/(float(nk))
             write(33,3002) x, z, v(i,JMONO,k)
             write(32,3002) x, z, t(i,JMONO,k)
 1372  continue

	   CLOSE(8)
	   CLOSE(9)
CC	   CLOSE(10)
	   CLOSE(32)
	   CLOSE(33)
CC	   CLOSE(37)
CC
      ENDIF
 3001 FORMAT(I10,1X,F12.8)
 3002 FORMAT(3E12.5)
 3010 FORMAT(A80)
 3030 FORMAT(6E12.5)
 3040 FORMAT(2(F10.4,1X))
 3050 FORMAT(5(F10.4,1X))

      xtime = xtime + dtime                                                     
c                                                                               
c        START CALCULATIONS                                                     
c                                                                               
      iter = 0                                                                  
      jterm = 0                                                                 
      jjterm = 0                                                                
c                                                                               
c     define the updated tpd(i,j,k), upd(i,j,k) and vpd(i,j,k)                  
c                             for  su(i,j,k)                                    
c                                                                               
      do 48 k=1,nkp1                                                            
         do 48 j=1,njp1                                                         
            do 48 i=1,nip1                                                      
               tpd(i,j,k) = t(i,j,k)                                            
               upd(i,j,k) = u(i,j,k)                                            
               vpd(i,j,k) = v(i,j,k)                                            
               wpd(i,j,k) = w(i,j,k)                                            
   48 continue                                                                  
c                                                                               
   29 continue                                                                  
      jterm = jterm + 1                                                         
c                                                                               
      call calt                                                                 
c                                                                               
      do 2220 i=1,nip1                                                          
         do 2220 k=1,nkp1                                                       
            do 2220 j=2,nj                                                      
               IF (T(I,J,K) .LT. TCOOL) T(I,J,K) = TCOOL                        
               IF (T(I,J,K). GT. THOT) T(I,J,K) = THOT                          
 2220 continue                                                                  
c                                                                               
c        PRESSURE CORRECTION LOOP                                               
c                                                                               
  301 continue                                                                  
      ITER = ITER + 1                                                           
c                                                                               
      call calu                                                                 
c                                                                               
      call calv                                                                 
c                                                                               
      call calw                                                                 
c                                                                               
      call calp                                                                 
c                                                                               
      if(resorm(iter) .le. sormax) go to 49                                     
      if(iter .eq. 1) go to 302                                                 
      if(resorm(iter) .le. resorm(iter-1)) go to 302                            
      go to 304                                                                 
  302 if(jterm .ge. 2) go to 37                                                 
      source=resorm(iter)                                                       
      go to 39                                                                  
   37 if(resorm(iter) .le. source) go to 38                                     
      go to 304                                                                 
   38 source=resorm(iter)                                                       
   39 continue                                                                  
      do 23 k=1,nkp1                                                            
         do 23 j=1,njp1                                                         
            do 23 i=1,nip1                                                      
               tpd(i,j,k) = t(i,j,k)                                            
               upd(i,j,k) = u(i,j,k)                                            
               vpd(i,j,k) = v(i,j,k)                                            
               wpd(i,j,k) = w(i,j,k)                                            
               ppd(i,j,k) = p(i,j,k)                                            
   23 continue                                                                  
      jjterm=0                                                                  
      if(iter .eq. itmax) go to 49                                              
      if(jterm .eq. 2) go to 35                                                 
      if(iter .eq. 4) go to 29                                                  
   35 continue                                                                  
      if(jterm .eq. 3) go to 58                                                 
      if(iter .eq. 7) go to 29                                                  
   58 continue                                                                  
      jjterm=0                                                                  
      go to 301                                                                 
  304 continue                                                                  
      jjterm=jjterm+1                                                           
      if(jterm .eq. 1) go to 41                                                 
      if(jterm .eq. 2 .and. jjterm .eq. 1 .and. iter .ne. 5) go to 41           
      go to 82                                                                  
   41 continue                                                                  
      do 40 k=1,nkp1                                                            
         do 40 j=1,njp1                                                         
            do 40 i=1,nip1                                                      
               u(i,j,k) = upd(i,j,k)                                            
               v(i,j,k) = vpd(i,j,k)                                            
               w(i,j,k) = wpd(i,j,k)                                            
               p(i,j,k) = ppd(i,j,k)                                            
   40 continue                                                                  
      if(iter .eq. itmax) go to 49                                              
      go to 29                                                                  
   82 continue                                                                  
      do 43 k=1,nkp1                                                            
         do 43 j=1,njp1                                                         
            do 43 i=1,nip1                                                      
               t(i,j,k) = tpd(i,j,k)                                            
               u(i,j,k) = upd(i,j,k)                                            
               v(i,j,k) = vpd(i,j,k)                                            
               w(i,j,k) = wpd(i,j,k)                                            
               p(i,j,k) = ppd(i,j,k)                                            
   43 continue                                                                  
      if(iter .eq. itmax) go to 49                                              
      if((jterm .eq. 3 .and. iter .ne. 8) .or. jjterm .eq. 2) go to 49          
      go to 301                                                                 
   49 continue                                                                  
c                                                                               
      II = NIP1 / 3                                                             
      JJ = NJP1 / 3                                                             
      KK = NKP1 / 3                                                             
      uu(nt) = u(ii,jj,kk)                                                      
      vv(nt) = v(ii,jj,kk)                                                      
      ww(nt) = w(ii,jj,kk)                                                      

c                                                                               
c                        PRINT OUTPUT                                           
c                                                                               
c        SORSUM IS THE SUM OF "SMP'S" FROM ALL CONTROL VOLUMES                  
c                                                                               
      IF (MOD(NT,1) .EQ. 0) THEN                                                
         call nu(rnuc,rnuh)                                                     
         write(6,500) nt,xtime,iter,resorm(iter),sorsum,                        
     &                rnuc,rnuh                                                 
  500    format(1x, 'nt =',i9,5x,'time=',f6.4,/,                                
     &          1x, 'iter=',i2,  5x,'source=',f9.6,5x,'sorsum=',f9.6,/,         
     &          1x, 'nuc=',f10.6, 5x,'nuh=',f10.6,/)    

c TRACE calcula la trayectoria de una particula en el tiempo
         call TRACE(xa,ya,za) 

      ENDIF 

CC    
CC  Monitores puntuales de la velocidad
CC                                                                        
      IF (MOD(NT,NPRT).EQ.0) THEN 

C       WRITE(11,3070) XTIME,(T(NIP1,NJP1/2,NKP1/2)
C     &                      +T(NIP2,NJP1/2,NKP1/2))/2.
C     &                     -(T(1,NJP1/2,NKP1/2)
C     &                      +T(2,NJP1/2,NKP1/2))/2.

       WRITE(12,3070) XTIME,V(IMON1,JMON1,KMON1)
       WRITE(13,3070) XTIME,V(IMON2,JMON2,KMON2)
       WRITE(34,3070) XTIME,U(IMON2,JMON2,KMON2)
       WRITE(35,3070) XTIME,W(IMON2,JMON2,KMON2)
       WRITE(14,3070) XTIME,U(IMON3,JMON3,KMON3)
       WRITE(15,3070) XTIME,U(IMON4,JMON4,KMON4)
       WRITE(16,3070) XTIME,V(IMON4,JMON4,KMON4)
       WRITE(17,3070) XTIME,W(IMON4,JMON4,KMON4)
       WRITE(18,3070) XTIME,RNUH
       WRITE(19,3070) XTIME,RNUC
       WRITE(41,3070) U(ILM,JLM,KLM),V(ILM,JLM,KLM)
       WRITE(42,3070) U(ILM,JLM,KLM),W(ILM,JLM,KLM)
       WRITE(43,3070) V(ILM,JLM,KLM),W(ILM,JLM,KLM)
       WRITE(44,3002) U(ILM,JLM,KLM),V(ILM,JLM,KLM),W(ILM,JLM,KLM)
      ENDIF
 3070 FORMAT(2(F12.6,1X))
CC
c        UPDATE VARIABLES FOR THE NEXT TIME STEP                                
c                                
      do 305 k=1,nkp1                                                           
         do 305 j=1,njp1                                                        
            do 305 i=1,nip1                                                     
               tod(i,j,k) = t(i,j,k)                                            
               uod(i,j,k) = u(i,j,k)                                            
               vod(i,j,k) = v(i,j,k)                                            
               wod(i,j,k) = w(i,j,k) 
  305 continue      
      
      if(MOD(NT,NPR) .eq. 0) then
       rewind 30
       do 1371 k=1,nkp1
	 do 1371 j=1,njp1
	   do 1371 i=1,nip1
	    write(30,3030) t(i,j,k),u(i,j,k),v(i,j,k),
     &				    w(i,j,k),p(i,j,k)
 1371  continue
       rewind 31
       write(31,3001) nt,xtime
      endif

      go to 300

CC  Hay que cerrar los archivos abiertos

      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(20)
      CLOSE(30)
      CLOSE(31)
      CLOSE(34)
      CLOSE(35)
      end                                                                       

c *********************************************************************         
      subroutine calt                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      common/bl7/ni,nip1,nim1,nj,njp1,njm1,nk,nkp1,nkm1                         
     &  ,NIP2,NJP2,NKP2,ITER,NNMAX                                              
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/smp(24,24,24),resorm(93),thot,tcool,                          
     &     du(24,24,24),dv(24,24,24)                                            
     &    ,dw(24,24,24),pp(24,24,24)                                            
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &       aw(24,24,24),an(24,24,24),                                         
     &        as(24,24,24),af(24,24,24),ab(24,24,24)                            
     &   ,sp(24,24,24),su(24,24,24)                                             
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
c                                                                               
c         CALCULATE COEFFICIENTS                                                
c                                                                               
      do 100 k=2,nk                                                             
      do 100 j=2,nj                                                             
      do 100 i=2,ni                                                             
c                                                                               
      condn1 = dzx / dy                                                         
      conds1 = dzx / dy                                                         
      condf1 = dxy / dz                                                         
      condb1 = dxy / dz                                                         
      conde1 = dyz / dx                                                         
      condw1 = dyz / dx                                                         
c                                                                               
      ce = u(i+1,j,k) * dyz                                                     
      cw = u(i  ,j,k) * dyz                                                     
      cn = v(i,j+1,k) * dzx                                                     
      cs = v(i,j  ,k) * dzx                                                     
      cf = w(i,j,k+1) * dxy                                                     
      cb = w(i,j,k  ) * dxy                                                     
c                                                                               
      if (ce .gt. 0.) then                                                      
          cep = ce * 0.125                                                      
          cem = 0.                                                              
      else                                                                      
          cep = 0.                                                              
          cem = - ce * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cw .gt. 0.) then                                                      
          cwp = cw * 0.125                                                      
          cwm = 0.                                                              
      else                                                                      
          cwp = 0.                                                              
          cwm = - cw * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cn .gt. 0.) then                                                      
          cnp = cn * 0.125                                                      
          cnm = 0.                                                              
      else                                                                      
          cnp = 0.                                                              
          cnm = - cn * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cs .gt. 0.) then                                                      
          csp = cs * 0.125                                                      
          csm = 0.                                                              
      else                                                                      
          csp = 0.                                                              
          csm = - cs * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cf .gt. 0.) then                                                      
          cfp = cf * 0.125                                                      
          cfm = 0.                                                              
      else                                                                      
          cfp = 0.                                                              
          cfm = - cf * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cb .gt. 0.) then                                                      
          cbp = cb * 0.125                                                      
          cbm = 0.                                                              
      else                                                                      
          cbp = 0.                                                              
          cbm = - cb * 0.125                                                    
      endif                                                                     
c                                                                               
      ae(i,j,k) = - 0.5 * ce + cep + conde1                                     
     &          + 2. * cem + cwm                                                
      aw(i,j,k) = 0.5 * cw + cwm + condw1                                       
     &          + 2. * cwp + cep                                                
      an(i,j,k) = - 0.5 * cn + cnp + condn1                                     
     &          + 2. * cnm + csm                                                
      as(i,j,k) = 0.5 * cs + csm + conds1                                       
     &          + 2. * csp + cnp                                                
      af(i,j,k) = - 0.5 * cf + cfp + condf1                                     
     &          + 2. * cfm + cbm                                                
      ab(i,j,k) = 0.5 * cb + cbm + condb1                                       
     &          + 2. * cbp + cfp                                                
c                                                                               
      if (i .lt. nim1) then                                                     
          aee = - cem                                                           
          aeer = aee * tpd(i+2,j,k)                                             
      else if (i .eq. nim1) then                                                
          aee = - cem                                                           
          aeer = aee * tpd(i+1,j,k)                                             
      else                                                                      
          aee = 0.                                                              
          aeer = 0.                                                             
      endif                                                                     
c                                                                               
      if (i .gt. 3) then                                                        
          aww = - cwp                                                           
          awwr = aww * tpd(i-2,j,k)                                             
      else if (i.eq.3) then                                                     
          aww = - cwp                                                           
          awwr = aww * tpd(i-1,j,k)                                             
      else                                                                      
          aww = 0.                                                              
          awwr = 0.                                                             
      endif                                                                     
c                                                                               
      if (j .lt. njm1) then                                                     
          ann = - cnm                                                           
          annr = ann * tpd(i,j+2,k)                                             
      else if (j .eq. njm1) then                                                
          ann = - cnm                                                           
          annr = ann * (2. * tcool - tpd(i,j+1,k))                              
      else                                                                      
          ann = 0.                                                              
          annr = 0.                                                             
      endif                                                                     
c                                                                               
      if (j .gt. 3) then                                                        
          ass = - csp                                                           
          assr = ass * tpd(i,j-2,k)                                             
      else if (j .eq. 3) then                                                   
          ass = - csp                                                           
          assr = ass * (2. * thot - tpd(i,j-1,k))                               
      else                                                                      
          ass = 0.                                                              
          assr = 0.                                                             
      endif                                                                     
c                                                                               
      if (k .lt. nkm1) then                                                     
          aff = - cfm                                                           
          affr = aff * tpd(i,j,k+2)                                             
      else if (k .eq. nkm1) then                                                
          aff = - cfm                                                           
          affr = aff * tpd(i,j,k+1)                                             
      else                                                                      
          aff = 0.                                                              
          affr = 0.                                                             
      endif                                                                     
c                                                                               
      if (k .gt. 3) then                                                        
          abb = - cbp                                                           
          abbr = abb * tpd(i,j,k-2)                                             
      else if (k.eq.3) then                                                     
          abb = - cbp                                                           
          abbr = abb * tpd(i,j,k-1)                                             
      else                                                                      
          abb = 0.                                                              
          abbr = 0.                                                             
      endif                                                                     
c                                                                               
      ap(i,j,k) = ae(i,j,k) + aw(i,j,k) + an(i,j,k) + as(i,j,k)                 
     &          + af(i,j,k) + ab(i,j,k) + aee + aww + ann                       
     &          + ass + aff + abb                                               
      sp(i,j,k) = - voldt                                                       
      su(i,j,k) = voldt * tod(i,j,k)                                            
      su(i,j,k) = su(i,j,k) + aeer + awwr + annr + assr                         
     &          + affr + abbr                                                   
  100 continue                                                                  
c                                                                               
c        ISOTHERMAL WALLS                                                       
c                                                                               
c        BOTTOM AND TOP WALLS                                                   
c                                                                               
      do 101 k=2,nk                                                             
         DO 101 I=2,NI                                                          
            sp(i,2,k) = sp(i,2,k) - as(i,2,k)                                   
            su(i,2,k) = su(i,2,k) + 2. * as(i,2,k) * thot                       
            AS(I,2,K) = 0.                                                      
            sp(i,nj,k) = sp(i,nj,k) - an(i,nj,k)                                
            su(i,nj,k) = su(i,nj,k) + 2. * an(i,nj,k) * tcool                   
            an(i,nj,k) = 0.                                                     
  101 continue                                                                  
c                                                                               
c        ADIABATIC WALLS (LEFT AND RIGHT)                                       
c                                                                               
      do 102 k=2,nk                                                             
         DO 102 J=2,NJ                                                          
            sp(2,j,k) = sp(2,j,k) + aw(2,j,k)                                   
            aw(2,j,k) = 0.                                                      
            sp(ni,j,k) = sp(ni,j,k) + ae(ni,j,k)                                
            ae(ni,j,k) = 0.                                                     
  102 continue                                                                  
c                                                                               
c        ADIABATIC WALLS (BACK AND FRONT)                                       
c                                                                               
      do 103 j=2,nj                                                             
         DO 103 I=2,NI                                                          
            sp(i,j,2) = sp(i,j,2) + ab(i,j,2)                                   
            ab(i,j,2) = 0.                                                      
            sp(i,j,nk) = sp(i,j,nk) + af(i,j,nk)                                
            af(i,j,nk) = 0.                                                     
  103 continue                                                                  
c                                                                               
C       ASSEMBLE COEFFICIENTS AND SOLVE DIFFERENCE EQUATION                     
c                                                                               
      do 300 k=2,nk                                                             
         do 300 j=2,nj                                                          
            do 300 i=2,ni                                                       
               ap(i,j,k) = ap(i,j,k) - sp(i,j,k)                                
  300 continue                                                                  
      call sip (2,2,2,ni,nj,nk,t)                                               
c                                                                               
c     PRESCRIBE WALL TEMPERATURES CONSISTENT WITH THE                           
c                 BOUNDARY CONDITIONS                                           
c                                                                               
c             RIGHT AND LEFT WALL                                               
c                                                                               
      do 310 k=2,nk                                                             
         do 310 j=2,nj                                                          
            t(nip1,j,k) = t(ni,j,k)                                             
            t(1,j,k) = t(2,j,k)                                                 
  310 continue                                                                  
c                                                                               
c             FRONT AND BACK WALL                                               
c                                                                               
      do 320 j=2,nj                                                             
         do 320 i=2,ni                                                          
            t(i,j,nkp1) = t(i,j,nk)                                             
            t(i,j,1) = t(i,j,2)                                                 
  320 continue                                                                  
c                                                                               
c            TOP AND BOTTOM WALL                                                
c                                                                               
      do 330 k=2,nk                                                             
         DO 330 I=2,NI
            T(I,1,K) = 2. * THOT - T(I,2,K)                                     
            T(I,NJP1,K) = 2. * TCOOL - T(I,NJ,K)                                
  330 continue                                                                  
      end                                                                       
c**********************************************************************         
      subroutine calu                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/smp(24,24,24),resorm(93),thot,tcool,                          
     &     du(24,24,24),dv(24,24,24)                                            
     &    ,dw(24,24,24),pp(24,24,24)                                            
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &        aw(24,24,24),an(24,24,24),                                        
     &        as(24,24,24),af(24,24,24),ab(24,24,24)                            
     &   ,sp(24,24,24),su(24,24,24)                                             
      COMMON/BL101/GX1U                                                         
      COMMON/BL102/RAR,OX,OY,OZ,OMEGAX0,OMEGAY0,OMAGS                           
      COMMON/BL103/TA                                                           
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
c                                                                               
c ***    calculate coefficients                                                 
c                                                                               
      do 100 k=2,nk                                                             
      do 100 j=2,nj                                                             
      do 100 i=3,ni                                                             
c                                                                               
CC    Calculate parameters for centrifugal acceleration
CC
      FI=FLOAT(I)
      FJ=FLOAT(J)
      FK=FLOAT(K)
CC
      XR=(FI-2.0)*DX
      YR=(FJ-1.5)*DY
      ZR=(FK-1.5)*DZ

CC
      SK=((XR-OMEGAX0)*OX+(YR-OMEGAY0)*OY+ZR*OZ)/OMAGS
CC
      XR0=SK*OX+OMEGAX0
      YR0=SK*OY+OMEGAY0
      ZR0=SK*OZ
CC
      CENTX=(YR-YR0)*OX*OY-(XR-XR0)*OY*OY
     &     +(ZR-ZR0)*OX*OZ-(XR-XR0)*OZ*OZ
CC   
      zxoyn = dzx / dy                                                          
      zxoys = dzx / dy                                                          
      xyozf = dxy / dz                                                          
      xyozb = dxy / dz                                                          
      yzoxe = dyz / dx                                                          
      yzoxw = dyz / dx                                                          
c                                                                               
      GN = V(I,J+1,K)                                                           
      gnw = v(i-1,j+1,k)                                                        
      GS = V(I,J,K)                                                             
      GSW = V(I-1,J,K)                                                          
c                                                                               
      ge = u(i+1,j,k)                                                           
      GP = U(I,J,K)                                                             
      gw = u(i-1,j,k)                                                           
c                                                                               
      GF = W(I,J,K+1)                                                           
      gfw = w(i-1,j,k+1)                                                        
      GB = W(I,J,K)                                                             
      GBW = W(I-1,J,K)                                                          
c                                                                               
CC    Calculate parameters for Coriolis acceleration
CC
      VP=(GN+GNW+GS+GSW)/4.
      WP=(GF+GFW+GB+GBW)/4.
CC
CC  Sign correct in the next expresion
CC
      CORIX=OY*WP-OZ*VP
CC
      cn = 0.5 * (gn + gnw) * dzx                                               
      cs = 0.5 * (gs + gsw) * dzx                                               
      ce = 0.5 * (ge + gp) * dyz                                                
      cw = 0.5 * (gp + gw) * dyz                                                
      cf = 0.5 * (gf + gfw) * dxy                                               
      CB = 0.5 * (GB + GBW) * DXY                                               
c                                                                               
      vise1 = yzoxe * pr                                                        
      visw1 = yzoxw * pr                                                        
      visn1 = zxoyn * pr                                                        
      viss1 = zxoys * pr                                                        
      visf1 = xyozf * pr                                                        
      visb1 = xyozb * pr                                                        
c                                                                               
      IF (CE .GT. 0.) THEN                                                      
        cep = ce * .125                                                         
        cem = 0.                                                                
      else                                                                      
        cep = 0.                                                                
        CEM = - CE * .125                                                       
      endif                                                                     
c                                                                               
      IF (CW .GT. 0.) THEN                                                      
        cwp = cw * .125                                                         
        cwm = 0.                                                                
      else                                                                      
        cwp = 0.                                                                
        cwm = - cw * .125                                                       
      endif                                                                     
c                                                                               
      if (cn .gt. 0.) then                                                      
          cnp = cn * 0.125                                                      
          cnm = 0.                                                              
      else                                                                      
          cnp = 0.                                                              
          cnm = - cn * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cs .gt. 0.) then                                                      
          csp = cs * .125                                                       
          csm = 0.                                                              
      else                                                                      
          csp = 0.                                                              
          csm = - cs * .125                                                     
      endif                                                                     
c                                                                               
      if (cf .gt. 0.) then                                                      
          cfp = cf * .125                                                       
          cfm = 0.                                                              
      else                                                                      
        cfp = 0.                                                                
        cfm = - cf * .125                                                       
      endif                                                                     
c                                                                               
      if (cb .gt. 0.) then                                                      
          cbp = cb * .125                                                       
          cbm = 0.                                                              
      else                                                                      
          cbp = 0.                                                              
          cbm = - cb * .125                                                     
      endif                                                                     
c                                                                               
c     BOUNDARY CONDITIONS                                                       
c                                                                               
      aee = - cem                                                               
      if (i .lt. ni) then                                                       
          aeer = aee * upd(i+2,j,k)                                             
      else                                                                      
          aeer = - aee * upd(i-1,j,k)                                           
      endif                                                                     
c                                                                               
      aww = - cwp                                                               
      if (i .gt. 3) then                                                        
          awwr = aww * upd(i-2,j,k)                                             
      else                                                                      
          awwr = - aww * upd(i+1,j,k)                                           
      endif                                                                     
c                                                                               
      if (j .lt. njm1) then                                                     
          ann = - cnm                                                           
          annr = ann * upd(i,j+2,k)                                             
      else if (j .eq. njm1) then                                                
          ann = - cnm                                                           
          annr = - ann * upd(i,j+1,k)                                           
      else                                                                      
          ANN = 0.                                                              
          ANNR = 0.                                                             
      endif                                                                     
c                                                                               
      if (j .gt. 3) then                                                        
          ass = - csp                                                           
          assr = ass * upd(i,j-2,k)                                             
      else if (j.eq.3) then                                                     
          ass = - csp                                                           
          assr = - ass * upd(i,j-1,k)                                           
      else                                                                      
          ass = 0.                                                              
          assr = 0.                                                             
      endif                                                                     
c                                                                               
      if (k .lt. nkm1) then                                                     
          aff = - cfm                                                           
          affr = aff * upd(i,j,k+2)                                             
      else if (k .eq. nkm1) then                                                
          aff = - cfm                                                           
          affr = - aff * upd(i,j,k+1)                                           
      else                                                                      
          aff = 0.                                                              
          affr = 0.                                                             
      endif                                                                     
c                                                                               
      if (k .gt. 3) then                                                        
          abb = - cbp                                                           
          abbr = abb * upd(i,j,k-2)                                             
      else if (k .eq. 3) then                                                   
          abb = - cbp                                                           
          abbr = - abb * upd(i,j,k-1)                                           
      else                                                                      
          abb = 0.                                                              
          abbr = 0.                                                             
      endif                                                                     
c                                                                               
      ae(i,j,k) = - 0.5 * ce + cep + 2. * cem + cwm + vise1                     
      aw(i,j,k) = 0.5 * cw + 2. * cwp + cwm + cep + visw1                       
      an(i,j,k) = - 0.5 * cn + cnp + 2. * cnm + csm + visn1                     
      as(i,j,k) = 0.5 * cs + csm + 2. * csp + cnp + viss1                       
      af(i,j,k) = - 0.5 * cf + cfp + 2. * cfm + cbm + visf1                     
      ab(i,j,k) = 0.5 * cb + cbm + 2. * cbp + cfp + visb1                       
      ap(i,j,k) = ae(i,j,k) + aw(i,j,k) + an(i,j,k) + as(i,j,k)                 
     &          + af(i,j,k) + ab(i,j,k) + aee + aww                             
     &          + ann + ass + aff + abb                                         
      sp(i,j,k) = - voldt                                                       
      su(i,j,k) = voldt * uod(i,j,k)                                            
      su(i,j,k) = su(i,j,k) + dyz * ( p(i-1,j,k) - p(i,j,k) )                   
     &          + aeer + awwr + annr + assr + affr + abbr                       
     &          + 0.5 * RA * PR * (T(I,J,K) + T(I-1,J,K)) * VOL               
     &          * GX1U                                                          
     &          + 0.5 * RAR * PR * (T(I,J,K) + T(I-1,J,K)) 
     &          * CENTX * VOL
     &          + SQRT(TA) * PR * CORIX * VOL
  100 continue                                                                  
c                                                                               
c        BOUNDARY CONDITIONS                                                    
c                                                                               
c        TOP AND BOTTOM WALL                                                    
c                                                                               
      do 101 k=2,nk                                                             
         DO 101 I=3,NI                                                          
            sp(i,2,k) = sp(i,2,k) - as(i,2,k)                                   
            as(i,2,k) = 0.                                                      
            sp(i,nj,k) = sp(i,nj,k) - an(i,nj,k)                                
            an(i,nj,k) = 0.                                                     
  101 continue                                                                  
c                                                                               
c        LEFT AND RIGHT WALL                                                    
c                                                                               
      do 102 k=2,nk                                                             
         DO 102 J=2,NJ                                                          
            aw(3,j,k) = 0.                                                      
            ae(ni,j,k) = 0.                                                     
  102 continue                                                                  
c                                                                               
c       FRONT AND BACK WALL                                                     
c                                                                               
      do 103 j=2,nj                                                             
         DO 103 I=3,NI                                                          
            sp(i,j,2) = sp(i,j,2) - ab(i,j,2)                                   
            ab(i,j,2) = 0.                                                      
            sp(i,j,nk) = sp(i,j,nk) - af(i,j,nk)                                
            af(i,j,nk) = 0.                                                     
  103 continue                                                                  
c                                                                               
c        CALCULATE  AP                                                          
c                                                                               
      do 300 k=2,nk                                                             
         do 300 j=2,nj                                                          
            do 300 i=3,ni                                                       
               ap(i,j,k) = ap(i,j,k) - sp(i,j,k)                                
  300 continue                                                                  
c                                                                               
c       SOLVE THE LINEARISED U MOMENTUM EQUATIONS                               
c                                                                               
      call sip (3,2,2,ni,nj,nk,u)                                               
c                                                                               
c       CALCULATE DU NEEDED FOR PRESSURE CORRECTION                             
c                                                                               
      do 301 k=2,nk                                                             
         DO 301 J=2,NJ                                                          
            DO 301 I=3,NI                                                       
               SU(I,J,K) = DYZ                                                  
  301 continue                                                                  
c                                                                               
c      SOLVE EQUATION FOR DU                                                    
c                                                                               
      call sip (3,2,2,ni,nj,nk,du)                                              
      end                                                                       
c**********************************************************************         
      subroutine calv                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/smp(24,24,24),resorm(93),thot,tcool,                          
     &     du(24,24,24),dv(24,24,24)                                            
     &    ,dw(24,24,24),pp(24,24,24)                                            
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &        aw(24,24,24),an(24,24,24),                                        
     &        as(24,24,24),af(24,24,24),ab(24,24,24)                            
     &   ,sp(24,24,24),su(24,24,24)                                             
      COMMON/BL201/GY1U                                                         
      COMMON/BL102/RAR,OX,OY,OZ,OMEGAX0,OMEGAY0,OMAGS                           
      COMMON/BL103/TA                                                           
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
c                                                                               
c ***    calculate coefficients                                                 
c                                                                               
      do 100 k=2,nk                                                             
      do 100 j=3,nj                                                             
      do 100 i=2,ni                                                             
c                                                                               
CC    Calculate parameters for centrifugal acceleration
CC
      FI=FLOAT(I)
      FJ=FLOAT(J)
      FK=FLOAT(K)
CC
      XR=(FI-1.5)*DX
      YR=(FJ-2.0)*DY
      ZR=(FK-1.5)*DZ
CC
      SK=((XR-OMEGAX0)*OX+(YR-OMEGAY0)*OY+ZR*OZ)/OMAGS
CC
      XR0=SK*OX+OMEGAX0
      YR0=SK*OY+OMEGAY0
      ZR0=SK*OZ
CC
      CENTY=-((YR-YR0)*OX*OX-(XR-XR0)*OX*OY
     &       -(ZR-ZR0)*OY*OZ+(YR-YR0)*OZ*OZ)
CC   
      zxoyn = dzx / dy                                                          
      zxoys = dzx / dy                                                          
      xyozf = dxy / dz                                                          
      xyozb = dxy / dz                                                          
      yzoxe = dyz / dx                                                          
      yzoxw = dyz / dx                                                          
c                                                                               
      gn = v(i,j+1,k)                                                           
      gp = v(i,j,k)                                                             
      gs = v(i,j-1,k)                                                           
      ge = u(i+1,j,k)                                                           
      gse = u(i+1,j-1,k)                                                        
      gw = u(i,j,k)                                                             
      gsw = u(i,j-1,k)                                                          
      gf = w(i,j,k+1)                                                           
      gsf = w(i,j-1,k+1)                                                        
      gb = w(i,j,k)                                                             
      gsb = w(i,j-1,k)                                                          
c                                                                               
CC    Calculate parameters for Coriolis acceleration
CC
      UP=(GE+GSE+GW+GSW)/4.
      WP=(GF+GSF+GB+GSB)/4.
CC
CC Sign correct in the next expresion
CC
      CORIY=OZ*UP-OX*WP
CC
      cn = 0.5 * (gn + gp) * dzx                                                
      cs = 0.5 * (gp + gs) * dzx                                                
      ce = 0.5 * (ge + gse) * dyz                                               
      cw = 0.5 * (gw + gsw) * dyz                                               
      cf = 0.5 * (gf + gsf) * dxy                                               
      cb = 0.5 * (gb + gsb) * dxy                                               
c                                                                               
      vise1 = yzoxe * pr                                                        
      visw1 = yzoxw * pr                                                        
      visn1 = zxoyn * pr                                                        
      viss1 = zxoys * pr                                                        
      visf1 = xyozf * pr                                                        
      visb1 = xyozb * pr                                                        
c                                                                               
      if (ce .gt. 0.) then                                                      
          cep = ce * 0.125                                                      
          cem = 0.                                                              
      else                                                                      
          cep = 0.                                                              
          cem = - ce * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cw .gt. 0.) then                                                      
          cwp = cw * 0.125                                                      
          cwm = 0.                                                              
      else                                                                      
          cwp = 0.                                                              
          cwm = - cw * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cn .gt. 0.) then                                                      
          cnp = cn * 0.125                                                      
          cnm = 0.                                                              
      else                                                                      
          cnp = 0.                                                              
          cnm = - cn * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cs .gt. 0.) then                                                      
          csp = cs * 0.125                                                      
          csm = 0.                                                              
      else                                                                      
          csp = 0.                                                              
          csm = - cs * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cf .gt. 0.) then                                                      
          cfp = cf * 0.125                                                      
          cfm = 0.                                                              
      else                                                                      
          cfp = 0.                                                              
          cfm = - cf * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cb .gt. 0.) then                                                      
          cbp = cb * 0.125                                                      
          cbm = 0.                                                              
      else                                                                      
          cbp = 0.                                                              
          cbm = - cb * 0.125                                                    
      endif                                                                     
c                                                                               
      if (i .lt. nim1) then                                                     
          aee = - cem                                                           
          aeer = aee * vpd(i+2,j,k)                                             
      else if (i .eq. nim1) then                                                
          aee = - cem                                                           
          aeer = - aee * vpd(i+1,j,k)                                           
      else                                                                      
          aee = 0.                                                              
          aeer = 0.                                                             
      endif                                                                     
c                                                                               
      if (i .gt. 3) then                                                        
          aww = - cwp                                                           
          awwr = aww * vpd(i-2,j,k)                                             
      else if (i.eq.3) then                                                     
          aww = -cwp                                                            
          awwr = - aww * vpd(i-1,j,k)                                           
      else                                                                      
          aww = 0.                                                              
          awwr = 0.                                                             
      endif                                                                     
c                                                                               
      ann = - cnm                                                               
      if (j .lt. nj) then                                                       
          annr = ann * vpd(i,j+2,k)                                             
      else                                                                      
          annr = - ann * vpd(i,j-1,k)                                           
      endif                                                                     
c                                                                               
      ass = - csp                                                               
      if (j .gt. 3) then                                                        
          assr = ass * vpd(i,j-2,k)                                             
      else                                                                      
          assr = - ass * vpd(i,j+1,k)                                           
      endif                                                                     
c                                                                               
      if (k .lt. nkm1) then                                                     
          aff = - cfm                                                           
          affr = aff * vpd(i,j,k+2)                                             
      else if (k .eq. nkm1) then                                                
          aff = - cfm                                                           
          affr = - aff * vpd(i,j,k+1)                                           
      else                                                                      
          aff = 0.                                                              
          affr = 0.                                                             
      endif                                                                     
c                                                                               
      if (k .gt. 3) then                                                        
          abb = - cbp                                                           
          abbr = abb * vpd(i,j,k-2)                                             
      else if (k .eq. 3) then                                                   
          abb = - cbp                                                           
          abbr = - abb * vpd(i,j,k-1)                                           
      else                                                                      
          abb = 0.                                                              
          abbr = 0.                                                             
      endif                                                                     
c                                                                               
      ae(i,j,k) = - 0.5 * ce + cep + 2. * cem + cwm + vise1                     
      aw(i,j,k) = 0.5 * cw + cwm +2. * cwp + cep + visw1                        
      an(i,j,k) = - 0.5 * cn + cnp + 2. * cnm + csm + visn1                     
      as(i,j,k) = 0.5 * cs + 2. * csp + csm + cnp + viss1                       
      af(i,j,k) = - 0.5 * cf + cfp + 2. * cfm + cbm + visf1                     
      ab(i,j,k) = 0.5 * cb + cbm + 2. * cbp + cfp + visb1                       
      ap(i,j,k) = ae(i,j,k) + aw(i,j,k) + an(i,j,k) + as(i,j,k)                 
     &          + af(i,j,k) + ab(i,j,k) + aee + aww + ann                       
     &          + ass + aff + abb                                               
      sp(i,j,k) = - voldt                                                       
      su(i,j,k) = voldt * vod(i,j,k)                                            
      su(i,j,k) = su(i,j,k) + dzx * (p(i,j-1,k) - p(i,j,k))                     
     &          + aeer + awwr + annr + assr + affr + abbr                       
     &          + 0.5 * ra * pr * (t(i,j,k) + t(i,j-1,k)) * vol                 
     &          * GY1U                                                          
     &          + 0.5 * RAR * PR * (T(I,J,K) + T(I,J-1,K)) 
     &          * CENTY * VOL
     &          + SQRT(TA) * PR * CORIY * VOL
  100 continue                                                                  
c                                                                               
c        BOUNDARY CONDITIONS                                                    
c                                                                               
c        TOP AND BOTTOM WALL                                                    
c                                                                               
      do 101 k=2,nk                                                             
         DO 101 I=2,NI                                                          
            as(i,3,k) = 0.                                                      
            an(i,nj,k) = 0.                                                     
  101 continue                                                                  
c                                                                               
c         LEFT AND RIGHT WALL                                                   
c                                                                               
      do 102 k=2,nk                                                             
         DO 102 J=3,NJ                                                          
            sp(2,j,k) = sp(2,j,k) - aw(2,j,k)                                   
            aw(2,j,k) = 0.                                                      
            sp(ni,j,k) = sp(ni,j,k) - ae(ni,j,k)                                
            ae(ni,j,k) = 0.                                                     
  102 continue                                                                  
c                                                                               
c        FRONT AND BACK WALL                                                    
c                                                                               
      do 103 j=3,nj                                                             
         DO 103 I=2,NI                                                          
            sp(i,j,2) = sp(i,j,2) - ab(i,j,2)                                   
            ab(i,j,2) = 0.                                                      
            sp(i,j,nk) = sp(i,j,nk) - af(i,j,nk)                                
            af(i,j,nk) = 0.                                                     
  103 continue                                                                  
c                                                                               
c        CALCULATE AP                                                           
c                                                                               
      do 300 k=2,nk                                                             
         do 300 j=3,nj                                                          
            do 300 i=2,ni                                                       
               ap(i,j,k) = ap(i,j,k) - sp(i,j,k)                                
  300 continue                                                                  
c                                                                               
c       SOLVE THE LINEARISED V MOMENTUM EQUATION                                
c                                                                               
      call sip (2,3,2,ni,nj,nk,v)                                               
c                                                                               
c       SET UP SU FOR CALCULATING DV                                            
c                                                                               
      do 301 k=2,nk                                                             
         DO 301 J=3,NJ                                                          
            DO 301 I=2,NI                                                       
               SU(I,J,K) = DZX                                                  
  301 continue                                                                  
c                                                                               
c     SOLVE EQUATION FOR CALCULATING DV FOR PRESSURE CORRECTION                 
c                                                                               
      call sip (2,3,2,ni,nj,nk,dv)                                              
      end                                                                       
c***********************************************************************        
      subroutine calw                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/smp(24,24,24),resorm(93),thot,tcool,                          
     &     du(24,24,24),dv(24,24,24)                                            
     &    ,dw(24,24,24),pp(24,24,24)                                            
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &        aw(24,24,24),an(24,24,24),                                        
     &        as(24,24,24),af(24,24,24),ab(24,24,24)                            
     &   ,sp(24,24,24),su(24,24,24)                                             
      COMMON/BL301/GZ1U                                                         
      COMMON/BL102/RAR,OX,OY,OZ,OMEGAX0,OMEGAY0,OMAGS                           
      COMMON/BL103/TA                                                           
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
c                                                                               
c ***        calculate coefficients                                             
c                                                                               
      do 100 k=3,nk                                                             
      do 100 j=2,nj                                                             
      do 100 i=2,ni                                                             
c                                                                               
CC    Calculate parameters for centrifugal acceleration
CC
      FI=FLOAT(I)
      FJ=FLOAT(J)
      FK=FLOAT(K)
CC
      XR=(FI-1.5)*DX
      YR=(FJ-1.5)*DY
      ZR=(FK-2.0)*DZ
CC
      SK=((XR-OMEGAX0)*OX+(YR-OMEGAY0)*OY+ZR*OZ)/OMAGS
CC
      XR0=SK*OX+OMEGAX0
      YR0=SK*OY+OMEGAY0
      ZR0=SK*OZ
CC
      CENTZ=(ZR-ZR0)*OX*OX-(XR-XR0)*OX*OZ
     &     -(ZR-ZR0)*OY*OY+(YR-YR0)*OY*OZ
CC   
      zxoyn = dzx / dy                                                          
      zxoys = dzx / dy                                                          
      xyozf = dxy / dz                                                          
      xyozb = dxy / dz                                                          
      yzoxe = dyz / dx                                                          
      yzoxw = dyz / dx                                                          
c                                                                               
      gn = v(i,j+1,k)                                                           
      gnb = v(i,j+1,k-1)                                                        
      gs = v(i,j,k)                                                             
      gsb = v(i,j,k-1)                                                          
      ge = u(i+1,j,k)                                                           
      geb = u(i+1,j,k-1)                                                        
      gw = u(i,j,k)                                                             
      gwb = u(i,j,k-1)                                                          
      gf = w(i,j,k+1)                                                           
      gp = w(i,j,k)                                                             
      gb = w(i,j,k-1)                                                           
c                                                                               
CC    Calculate parameters for Coriolis acceleration
CC
      UP=(GE+GEB+GW+GWB)/4.
      VP=(GN+GNB+GS+GSB)/4.
CC
      CORIZ=OX*VP-OY*UP
CC
      cn = 0.5 * (gn + gnb) * dzx                                               
      cs = 0.5 * (gs + gsb) * dzx                                               
      ce = 0.5 * (ge + geb) * dyz                                               
      cw = 0.5 * (gw + gwb) * dyz                                               
      cf = 0.5 * (gf + gp) * dxy                                                
      cb = 0.5 * (gb + gp) * dxy                                                
c                                                                               
      vise1 = yzoxe * pr                                                        
      visw1 = yzoxw * pr                                                        
      visn1 = zxoyn * pr                                                        
      viss1 = zxoys * pr                                                        
      visf1 = xyozf * pr                                                        
      visb1 = xyozb * pr                                                        
c                                                                               
      if (ce .gt. 0.) then                                                      
          cep = ce * 0.125                                                      
          cem = 0.                                                              
      else                                                                      
          cep = 0.                                                              
          cem = - ce * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cw .gt. 0.) then                                                      
          cwp = cw * 0.125                                                      
          cwm = 0.                                                              
      else                                                                      
          cwp = 0.                                                              
          cwm = - cw * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cn .gt. 0.) then                                                      
          cnp = cn * 0.125                                                      
          cnm = 0.                                                              
      else                                                                      
          cnp = 0.                                                              
          cnm = - cn * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cs .gt. 0.) then                                                      
          csp = cs * 0.125                                                      
          csm = 0.                                                              
      else                                                                      
          csp = 0.                                                              
          csm = - cs * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cf .gt. 0.) then                                                      
          cfp = cf * 0.125                                                      
          cfm = 0.                                                              
      else                                                                      
          cfp = 0.                                                              
          cfm = - cf * 0.125                                                    
      endif                                                                     
c                                                                               
      if (cb .gt. 0.) then                                                      
          cbp = cb * 0.125                                                      
          cbm = 0.                                                              
      else                                                                      
          cbp = 0.                                                              
          cbm = - cb * 0.125                                                    
      endif                                                                     
c                                                                               
      if (i .lt. nim1) then                                                     
          aee = - cem                                                           
          aeer = aee * wpd(i+2,j,k)                                             
      else if (i .eq. nim1) then                                                
          aee = - cem                                                           
          aeer = - aee * wpd(i+1,j,k)                                           
      else                                                                      
          aee = 0.                                                              
          aeer = 0.                                                             
      endif                                                                     
c                                                                               
      if (i .gt. 3) then                                                        
          aww = - cwp                                                           
          awwr = aww * wpd(i-2,j,k)                                             
      else if (i.eq.3) then                                                     
          aww = - cwp                                                           
          awwr = - aww * wpd(i-1,j,k)                                           
      else                                                                      
          aww = 0.                                                              
          awwr = 0.                                                             
      endif                                                                     
c                                                                               
      if (j .lt. njm1) then                                                     
          ann = - cnm                                                           
          annr = ann * wpd(i,j+2,k)                                             
      else if (j.eq.njm1) then                                                  
          ann = - cnm                                                           
          annr = - ann * wpd(i,j+1,k)                                           
      else                                                                      
          ann = 0.                                                              
          annr = 0.                                                             
      endif                                                                     
c                                                                               
      if (j .gt. 3) then                                                        
          ass = - csp                                                           
          assr = ass * wpd(i,j-2,k)                                             
      else if (j .eq. 3) then                                                   
          ass = - csp                                                           
          assr = - ass * wpd(i,j-2,k)                                           
      else                                                                      
          ass = 0.                                                              
          assr = 0.                                                             
      endif                                                                     
c                                                                               
      aff = - cfm                                                               
      if (k .lt. nk) then                                                       
          affr = aff * wpd(i,j,k+2)                                             
      else                                                                      
          affr = - aff * wpd(i,j,k-1)                                           
      endif                                                                     
c                                                                               
      abb = - cbp                                                               
      if (k .gt. 3) then                                                        
          abbr = abb * wpd(i,j,k-2)                                             
      else                                                                      
          abbr = - abb * wpd(i,j,k+1)                                           
      endif                                                                     
c                                                                               
      ae(i,j,k) = - 0.5 * ce + cep + 2. * cem + cwm + vise1                     
      aw(i,j,k) = 0.5 * cw + cwm + 2. * cwp + cep + visw1                       
      an(i,j,k) = - 0.5 * cn + cnp + 2. * cnm + csm + visn1                     
      as(i,j,k) = 0.5 * cs + csm + 2. * csp + cnp + viss1                       
      af(i,j,k) = - 0.5 * cf + cfp + 2. * cfm + cbm + visf1                     
      ab(i,j,k) = 0.5 * cb + 2. * cbp + cbm + cfp + visb1                       
      ap(i,j,k) = ae(i,j,k) + aw(i,j,k) + an(i,j,k) + as(i,j,k)                 
     &          + af(i,j,k) + ab(i,j,k) + aee + aww + ann                       
     &          + ass + aff + abb                                               
      sp(i,j,k) = - voldt                                                       
      su(i,j,k) = voldt * wod(i,j,k)                                            
      su(i,j,k) = su(i,j,k) + dxy * (p(i,j,k-1) - p(i,j,k))                     
     &          + aeer + awwr + annr + assr + affr + abbr                       
     &          + 0.5 * RA * PR * (T(I,J,K) + T(I,J,K-1)) * VOL                 
     &          * GZ1U                                                          
     &          + 0.5 * RAR * PR *(T(I,J,K) + T(I,J,K-1)) 
     &          * CENTZ * VOL
     &          + SQRT(TA) * PR * CORIZ * VOL
  100 continue                                                                  
c                                                                               
c         BOUNDARY CONDITIONS                                                   
c                                                                               
c           TOP AND BOTTOM WALL                                                 
c                                                                               
      do 101 k=3,nk                                                             
         DO 101 I=2,NI                                                          
            sp(i,2,k) = sp(i,2,k) - as(i,2,k)                                   
            as(i,2,k) = 0.                                                      
            sp(i,nj,k) = sp(i,nj,k) - an(i,nj,k)                                
            an(i,nj,k) = 0.                                                     
  101 continue                                                                  
c                                                                               
c          LEFT AND RIGHT WALL                                                  
c                                                                               
      do 102 k=3,nk                                                             
         DO 102 J=2,NJ                                                          
            sp(2,j,k) = sp(2,j,k) - aw(2,j,k)                                   
            aw(2,j,k) = 0.                                                      
            sp(ni,j,k) = sp(ni,j,k) - ae(ni,j,k)                                
            ae(ni,j,k) = 0.                                                     
  102 continue                                                                  
c                                                                               
c          FRONT AND BACK WALL                                                  
c                                                                               
      do 103 j=2,nj                                                             
         DO 103 I=2,NI                                                          
            ab(i,j,3) = 0.                                                      
            af(i,j,nk) = 0.                                                     
 103  continue                                                                  
c                                                                               
c     SET UP AP AND ESTIMATE DW FOR PRESSURE CORRECTION                         
c                                                                               
      do 300 k=3,nk                                                             
         do 300 j=2,nj                                                          
            do 300 i=2,ni                                                       
               ap(i,j,k) = ap(i,j,k) - sp(i,j,k)                                
  300 continue                                                                  
c                                                                               
c     SOLVE THE LINEARISED W MOMENTUM EQUATIONS                                 
c                                                                               
      call sip (2,2,3,ni,nj,nk,w)                                               
c                                                                               
c     SET UP SU FOR CALCULATING DW                                              
c                                                                               
      do 301 k=3,nk                                                             
         do 301 j=2,nj                                                          
            do 301 i=2,ni                                                       
               SU(I,J,K) = DXY                                                  
  301 continue                                                                  
c                                                                               
c     SOLVE EQUATIONS FOR DW REQUIRED IN PRESSURE CORRECTION                    
c                                                                               
      call sip (2,2,3,ni,nj,nk,dw)                                              
      end                                                                       
c***********************************************************************        
      subroutine calp                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl31/ tod(24,24,24),uod(24,24,24),vod(24,24,24),                   
     &             wod(24,24,24),pod(24,24,24)                                  
      common/bl32/ t(24,24,24)  ,u(24,24,24)  ,v(24,24,24),                     
     &             w(24,24,24)  ,p(24,24,24)                                    
      common/bl33/ tpd(24,24,24),upd(24,24,24),vpd(24,24,24),                   
     &             wpd(24,24,24),ppd(24,24,24)                                  
      common/bl34/smp(24,24,24),resorm(93),thot,tcool,                          
     &     du(24,24,24),dv(24,24,24)                                            
     &    ,dw(24,24,24),pp(24,24,24)                                            
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &        aw(24,24,24),an(24,24,24),                                        
     &        as(24,24,24),af(24,24,24),ab(24,24,24)                            
     &   ,sp(24,24,24),su(24,24,24)                                             
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
c                                                                               
c          SET UP COEFFICIENTS                                                  
c                                                                               
      do 100 k=2,nk                                                             
         do 100 j=2,nj                                                          
            do 100 i=2,ni                                                       
               an(i,j,k) = dzx * dv(i,j+1,k)                                    
               as(i,j,k) = dzx * dv(i,j,k)                                      
               ae(i,j,k) = dyz * du(i+1,j,k)                                    
               aw(i,j,k) = dyz * du(i,j,k)                                      
               af(i,j,k) = dxy * dw(i,j,k+1)                                    
               ab(i,j,k) = dxy * dw(i,j,k)                                      
               cn = v(i,j+1,k) * dzx                                            
               cs = v(i,j,k) * dzx                                              
               ce = u(i+1,j,k) * dyz                                            
               cw = u(i,j,k) * dyz                                              
               cf = w(i,j,k+1) * dxy                                            
               cb = w(i,j,k) * dxy                                              
               smp(i,j,k) = - ce + cw - cn + cs - cf + cb                       
               su(i,j,k) = smp(i,j,k)                                           
               sp(i,j,k) = 0.                                                   
  100 continue                                                                  
c                                                                               
c        BOUNDARY CONDITIONS                                                    
c                                                                               
c          FLOOR AND CEILING                                                    
c                                                                               
      do 101 k=2,nk                                                             
         do 101 i=2,ni                                                          
            as(i,2,k) = 0.                                                      
            an(i,nj,k) = 0.                                                     
  101 continue                                                                  
c                                                                               
c         LEFT AND RIGHT WALL                                                   
c                                                                               
      do 102 k=2,nk                                                             
         do 102 j=2,nj                                                          
            aw(2,j,k) = 0.                                                      
            ae(ni,j,k) = 0.                                                     
  102 continue                                                                  
c                                                                               
c         FRONT AND BACK WALL                                                   
c                                                                               
      do 103 i=2,ni                                                             
         do 103 j=2,nj                                                          
            ab(i,j,2) = 0.                                                      
            af(i,j,nk) = 0.                                                     
 103  continue                                                                  
c                                                                               
c       CALCULATE AP                                                            
c                                                                               
      do 200 j=2,nj                                                             
         do 200 i=2,ni                                                          
            do 200 k=2,nk                                                       
            AP(I,J,K) = AN(I,J,K) + AS(I,J,K)+ AE(I,J,K) + AW(I,J,K)            
     &                + AF(I,J,K) + AB(I,J,K) - SP(I,J,K)                       
  200 continue                                                                  
c                                                                               
c       SOLVE THE PRESSURE CORRECTION EQUATION                                  
c                                                                               
      call sip (2,2,2,ni,nj,nk,pp)                                              
c                                                                               
c       U VELOCITY CORRECTION                                                   
c                                                                               
      do 201 i=3,ni                                                             
         do 201 j=2,nj                                                          
            do 201 k=2,nk                                                       
         U(I,J,K) = U(I,J,K) + DU(I,J,K) * (PP(I-1,J,K) - PP(I,J,K))            
  201 continue                                                                  
c                                                                               
c       V VELOCITY CORRECTION                                                   
c                                                                               
      do 202 j=3,nj                                                             
         do 202 k=2,nk                                                          
            do 202 i=2,ni                                                       
         V(I,J,K) = V(I,J,K) + DV(I,J,K) * (PP(I,J-1,K) - PP(I,J,K))            
  202 continue                                                                  
c                                                                               
c      W VELOCITY CORRECTION                                                    
c                                                                               
      do 203 k=3,nk                                                             
         do 203 i=2,ni                                                          
            do 203 j=2,nj                                                       
         W(I,J,K) = W(I,J,K) + DW(I,J,K) * (PP(I,J,K-1) - PP(I,J,K))            
 203  continue                                                                  
c                                                                               
c      PRESSURE CORRECTION                                                      
c                                                                               
      do 204 j=2,nj                                                             
         do 204 i=2,ni                                                          
            do 204 k=1,nk                                                       
               p(i,j,k) = p(i,j,k) + pp(i,j,k)                                  
               PP(I,J,K) = 0.                                                   
  204 continue                                                                  
c                                                                               
c      RECALCULATE MASS FLUX ERRORS AFTER U,V,W,P,CORRECTIONS                   
c                                                                               
      sorsum = 0.                                                               
      resorm(iter) = 0.                                                         
      do 205 j=2,nj                                                             
         do 205 i=2,ni                                                          
            do 205 k=2,nk                                                       
               cn = v(i,j+1,k) * dzx                                            
               cs = v(i,j,k) * dzx                                              
               ce = u(i+1,j,k) * dyz                                            
               cw = u(i,j,k) * dyz                                              
               cf = w(i,j,k+1) * dxy                                            
               cb = w(i,j,k) * dxy                                              
               smp(i,j,k) = - ce + cw - cn + cs - cf + cb                       
c                                                                               
c        SORSUM IS THE ALGEBRAIC SUM OF THE MASS FLUX ERROR FROM                
c               ALL THE CONTROL VOLUMES                                         
c                                                                               
               sorsum = sorsum + smp(i,j,k)                                     
c                                                                               
c        RESORM IS THE SUM OF THE ABSOLUTE VALUE OF THE MASS FLUX               
c               ERROR FROM ALL CONTROL VOLUMES                                  
c                                                                               
               resorm(iter) = resorm(iter) + abs(smp(i,j,k))                    
  205 continue                                                                  
      end                                                                       
c **********************************************************************        
      subroutine nu(rnuc,rnuh)                                                  
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl32/ t(24,24,24)                                                  
     &           ,u(24,24,24),v(24,24,24)                                       
     &          ,w(24,24,24),p(24,24,24)                                        
      COMMON/ASP/AX,AZ,EPS                                                      
c                                                                               
      rnuh = 0.                                                                 
      rnuc = 0.                                                                 
      a = 1. / dy                                                               
      do 20 i=2,ni                                                              
        do 20 k=2,nk                                                            
          dtdyc = a * (t(i,2,k) - t(i,1,k))                                     
          dtdyh = a * (t(i,njp1,k) - t(i,nj,k))                                 
          rnuc = rnuc - dtdyc * dz * dx                                         
          rnuh = rnuh - dtdyh * dz * dx                                         
   20 continue                                                                  
      aa = 1. / (ax * az)                                                       
      rnuh = rnuh * aa                                                          
      rnuc = rnuc * aa                                                          
      end                                                                       
c***********************************************************************        
      block data                                                                
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      data nip2,nip1,ni,nim1/25,24,23,22/
      data njp2,njp1,nj,njm1/25,24,23,22/   
      data nkp2,nkp1,nk,nkm1/25,24,23,22/  
      DATA NNMAX/40000/                   
      end                                                                       
c***********************************************************************        
      subroutine sip(ist,jst,kst,isp,jsp,ksp,phi)                               
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/bl36/ap(24,24,24),ae(24,24,24),                                    
     &        aw(24,24,24),an(24,24,24),                                        
     &        as(24,24,24),af(24,24,24),ab(24,24,24),                           
     &        sp(24,24,24),su(24,24,24)                                         
      common/mean/t_mean(24,24,24),u_mean(24,24,24),                            
     &            v_mean(24,24,24),w_mean(24,24,24),                            
     &            p_mean(24,24,24)                                              
      COMMON/ASP/AX,AZ,EPS                                                      
      COMMON/ERROR/IMAX                                                         
      common/abc/kmax,n1,k1                                                     
      common/count/nt,nmax                                                      
      common/ mat/a1(40000),a2(40000),a3(40000),a4(40000),a5(40000),                 
     &            a6(40000),a7(40000),b(40000),d(40000)                             
      dimension x1(40000),x2(40000),xr1(40000),                                    
     &          s(24,24,24),phi(24,24,24)                                       
c                                                                               
      n1  = isp - ist + 1                                                       
      k1  = (isp - ist + 1) * (jsp - jst + 1)                                   
      kmax = (isp - ist + 1) * (jsp - jst + 1) * (ksp - kst + 1)                
      eps2 = kmax * eps**2                                                      
c                                                                               
c   DEFINING A SCRATCH ARRAY FOR FURTHER MODIFICATIONS                          
c                                                                               
      do 9 k=kst,ksp                                                            
         do 9 j=jst,jsp                                                         
            do 9 i=ist,isp                                                      
                 s(i,j,k) = su(i,j,k)                                           
    9 continue                                                                  
c                                                                               
c   EXTRACT MATRIX A, RHS VECTOR B AND VECTOR XU-AN INITIAL                     
c                     GUESS FROM GIVEN DATA                                     
c                                                                               
      do 2 k=kst,ksp                                                            
         do 2 j=jst,jsp                                                         
            do 2 i=ist,isp                                                      
               m  = i - ist + 1 + (j - jst) * n1 + (k - kst) * k1               
               a1(m) = ap(i,j,k)                                                
               a2(m) = -ae(i,j,k)                                               
               a3(m) = -aw(i,j,k)                                               
               a4(m) = -an(i,j,k)                                               
               a5(m) = -as(i,j,k)                                               
               a6(m) = -af(i,j,k)                                               
               a7(m) = -ab(i,j,k)                                               
               b(m) = s(i,j,k)                                                  
               x1(m) = phi(i,j,k)                                               
    2 continue                                                                  
c                                                                               
c       THE DIAGONAL ELEMENT OF L IN THE INCOMPLETE LU                          
c                     DECOMPOSITION                                             
c                                                                               
      d(1) = a1(1)                                                              
      do 113 i=2,n1                                                             
         D(I) = A1(I) - A2(I-1) * A3(I) / D(I-1)                                
  113 continue                                                                  
      do 114 i=n1+1,k1                                                          
         D(I) = A1(I) - A2(I-1) * A3(I) / D(I-1) - A4(I-N1)                     
     &        * A5(I) / D(I-N1)                                                 
  114 continue                                                                  
      do 115 i=k1+1,kmax                                                        
         D(I) = A1(I) - A2(I-1) * A3(I) / D(I-1)                                
     &                - A4(I-N1) * A5(I) / D(I-N1)                              
     &                - A6(I-K1) * A7(I) / D(I-K1)                              
  115 continue                                                                  
      do 116 i=1,kmax                                                           
         D(I) = 1. / D(I)                                                       
  116 continue                                                                  
      ITR = 0                                                                   
c                                                                               
c              THE MAIN ITERATION LOOP BEGINS                                   
c                                                                               
   7  CONTINUE                                                                  
      ITR = ITR + 1                                                             
      call res(x1,xr1)                                                          
      call xl(xr1,x2)                                                           
      do 117 i=1,kmax                                                           
         X2(I) = X2(I) + X1(I)                                                  
  117 continue                                                                  
      v2 = 0.                                                                   
      do 55 i=1,kmax                                                            
         V2 = V2 + XR1(I)**2                                                    
   55 continue                                                                  
      x2mag = 0.                                                                
      do 56 i=1,kmax                                                            
         X2MAG = X2MAG + X2(I)**2                                               
   56 continue                                                                  
      if (x2mag .gt. 0.0001) v2 = v2 / x2mag                                    
                                                                                
c                                                                               
c            RELABELING THE VECTORS BEFORE NEXT ITERATION                       
c                                                                               
      do 50 i=1,kmax                                                            
         X1(I) = X2(I)                                                          
  50  continue                                                                  
      IF (V2 .GT. EPS2 .AND. ITR .LT. IMAX) GO TO 7                             
c                                                                               
c    RECOVER THE SOLUTION VECTOR AND RETURN THE AVERAGE                         
c            PHI AS WELL AS THE RMS ERROR                                       
c                                                                               
      IF (MOD(NT,100) .EQ. 0) PRINT *, ITR,V2                                   
      do 3 k=kst,ksp                                                            
         do 3 j=jst,jsp                                                         
            do 3 i=ist,isp                                                      
               m  = i - ist + 1 + (j - jst) * n1 + (k - kst) * k1               
               phi(i,j,k) = x1(m)                                               
    3 continue                                                                  
      end                                                                       
c***********************************************************************        
      subroutine res(x1,x2)                                                     
      common/abc/kmax,n1,k1                                                     
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      common/ mat/a1(40000),a2(40000),a3(40000),a4(40000),a5(40000),                 
     &            a6(40000),a7(40000),b(40000),d(40000)                             
      dimension x1(40000),x2(40000)                                               
C                                                                               
c                                                                               
c     CALCULATE Ax1                                                             
c                                                                               
      x2(1) = a1(1) * x1(1) + a2(1) * x1(2)                                     
     &      + a4(1) * x1(1+n1) + a6(1) * x1(1+k1)                               
      do 10 i=2,n1                                                              
         x2(i) = a1(i) * x1(i) + a2(i) * x1(i+1)                                
     &         + a4(i) * x1(i+n1) + a6(i) * x1(i+k1)                            
     &         + a3(i) * x1(i-1)                                                
  10  continue                                                                  
      do 20 i=n1+1,k1                                                           
         x2(i) = a1(i) * x1(i) + a2(i) * x1(i+1)                                
     &         + a4(i) * x1(i+n1) + a5(i) * x1(i-n1)                            
     &         + a3(i) * x1(i-1) + a6(i) * x1(i+k1)                             
  20  continue                                                                  
      DO 30 I=K1+1,KMAX-K1                                                      
         x2(i) = a1(i) * x1(i) + a2(i) * x1(i+1)                                
     &         + a4(i) * x1(i+n1) + a5(i) * x1(i-n1)                            
     &         + a3(i) * x1(i-1) + a6(i) * x1(i+k1)                             
     &         + a7(i) * x1(i-k1)                                               
  30  continue                                                                  
      DO 32 I=KMAX-K1+1,KMAX-N1                                                 
         x2(i) = a1(i) * x1(i) + a2(i) * x1(i+1)                                
     &         + a4(i) * x1(i+n1) + a5(i) * x1(i-n1)                            
     &         + A3(I) * X1(I-1) + A7(I) * X1(I-K1)                             
  32  CONTINUE                                                                  
      DO 34 I=KMAX-N1+1,KMAX-1                                                  
         x2(i) = a1(i) * x1(i) + a2(i) * x1(i+1)                                
     &         + A5(I) * X1(I-N1) + A3(I) * X1(I-1)                             
     &         + A7(I) * X1(I-K1)                                               
  34  CONTINUE                                                                  
      X2(KMAX) = A1(KMAX) * X1(KMAX) + A5(KMAX) * X1(KMAX-N1)                   
     &         + A3(KMAX) * X1(KMAX-1) + A7(KMAX) * X1(KMAX-K1)                 
      do 40 i=1,kmax                                                            
         X2(I) = B(I) - X2(I)                                                   
  40  continue                                                                  
      end                                                                       
c***********************************************************************        
      subroutine xl(x1,x3)                                                      
      common/abc/kmax,n1,k1                                                     
      common/ mat/a1(40000),a2(40000),a3(40000),a4(40000),a5(40000),                 
     &            a6(40000),a7(40000),b(40000),d(40000)                             
      dimension x1(40000),x2(40000),x3(40000)                                      
      x2(1) = x1(1) * d(1)                                                      
      do 31 i=2,n1                                                              
         X2(I) = D(I) * (X1(I) - A3(I) * X2(I-1))                               
  31  continue                                                                  
      do 32 i=n1+1,k1                                                           
         X2(I) = D(I) * (X1(I) - A3(I) * X2(I-1)                                
     &         - A5(I) * X2(I-N1))                                              
  32  continue                                                                  
      do 33 i=k1+1,kmax                                                         
         X2(I) = D(I) * (X1(I) - A3(I) * X2(I-1)                                
     &         - A5(I) * X2(I-N1) - A7(I) * X2(I-K1))                           
  33  continue                                                                  
      x3(kmax) = x2(kmax)                                                       
      do 67 i=kmax-1,kmax-n1+1,-1                                               
         X3(I) = X2(I) - D(I) * A2(I) * X3(I+1)                                 
  67  continue                                                                  
      do 68 i=kmax-n1,kmax-k1+1,-1                                              
         X3(I) = X2(I) - D(I)*(A2(I)*X3(I+1) + A4(I)*X3(I+N1))                  
  68  continue                                                                  
      do 69 i=kmax-k1,1,-1                                                      
         X3(I) = X2(I) - D(I)*(A2(I)*X3(I+1) + A4(I)*X3(I+N1)                   
     &                      - a6(i)*x3(i+k1))                                   
  69  continue                                                                  
      end                                                                       
c***********************************************************************        
      subroutine grid                                                           
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum              
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,                        
     &           NIP2,NJP2,NKP2,ITER,NNMAX                                      
      COMMON/ASP/AX,AZ,EPS                                                      
c                                                                               
c       GENERATION OF GRID                                                      
c                                                                               
      dx = ax / float(nim1)                                                     
      dy = 1. / float(njm1)                                                     
      dz = az / float(nkm1)                                                     
c                                                                               
      vol = dx * dy * dz                                                        
      voldt = vol / dtime                                                       
c                                                                               
      dxy = dx * dy                                                             
      dyz = dy * dz                                                             
      dzx = dz * dx                                                             
CC      print *, 'X GRID SPACING IS',dx                                         
CC      print *, 'Y GRID SPACING IS',dy                                         
CC      print *, 'Z GRID SPACING IS',dz                                         
CC      print *, 'XY SURFACE ELEMENT IS',dxy                                    
CC      print *, 'YZ SURFACE ELEMENT IS',dyz                                    
CC      print *, 'ZX SURFACE ELEMENT IS',dzx                                    
CC      print *, 'VOLUME ELEMENT IS',vol                                        
c                                                                               
      end
c***********************************************************************
      subroutine TRACE(xa,ya,za)
      common/bl1/dx,dy,dz,dxy,dyz,dzx,vol,dtime,voldt,ra,pr,sorsum
      COMMON/ASP/AX,AZ,EPS
      COMMON/BL7/NI,NIP1,NIM1,NJ,NJP1,NJM1,NK,NKP1,NKM1,
     &           NIP2,NJP2,NKP2,ITER,NNMAX
      common/bl32/ t(24,24,24),u(24,24,24),v(24,24,24)
     &            ,w(24,24,24),p(24,24,24)

      DIMENSION VELIN(3)

      write(36,5010) xa,ya,za
      deltax = AX*2./NI
      deltay = 2./NJ
      deltaz = AZ*2./NK
    
      DDT = dtime 
      
      i = NINT(xa/deltax)
      j = NINT(ya/deltay)
      k = NINT(za/deltaz)
      ddx = MOD(xa,deltax)
      ddy = MOD(ya,deltay)
      ddz = MOD(za,deltaz)
      f1 = ddx/deltax
      f2 = ddy/deltay
      f3 = ddz/deltaz
 
      if(ddx .LT. deltax/2.) then
          iii = i + 1
      else
          iii = i - 1
      endif
      if(ddy .LT. deltay/2.) then
          jjj = j + 1
      else
          jjj = j - 1
      endif
      if(ddz .LT. deltaz/2.) then
          kkk = k + 1 
      else
          kkk = k - 1
      endif

      VELIN(1) = u(i,j,k)*(1-(f1 + f2 + f3))
     &           + u(iii,j,k)*f1 + u(i,jjj,k)*f2
     &           + u(i,j,kkk)*f3

      VELIN(2) = v(i,j,k)*(1-(f1 + f2 + f3))
     &           + v(iii,j,k)*f1 + v(i,jjj,k)*f2
     &           + v(i,j,kkk)*f3

      VELIN(3) = w(i,j,k)*(1-(f1 + f2 + f3))
     &           + w(iii,j,k)*f1 + w(i,jjj,k)*f2
     &           + w(i,j,kkk)*f3

      xs = xa + VELIN(1)*DDT
      ys = ya + VELIN(2)*DDT
      zs = za + VELIN(3)*DDT
      if( xs .GE. 1) xs = 1
      if( ys .GE. 2) ys = 2
      if( zs .GE. 1) zs = 1
      if( xs .LE. 0) xs = 0
      if( ys .LE. 0) ys = 0
      if( zs .LE. 0) zs = 0 
 5010 FORMAT(3(XF10.8))

      xa = xs
      ya = ys
      za = zs
       
      end
