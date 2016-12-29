program sde

    use mtmod
    use sde_par
    
    implicit none

    integer :: seed, i, j, k, i1, i2, j1, k1, k2
    
  9 FORMAT (1X,A13,1X,A20,1X,A12,1X,A16,1X,A23,1X,A9,1X,A16,1X,A21,1X,A11,1X,A16)
 10 FORMAT (1X,A13,1X,A24,1X,A8,1X,A16,1X,A22,1X,A10,1X,A15)  
 11 FORMAT (1X,A18,1X,A16,1X,A16,1X,A16,1X,A16,1X,A16,1X,A16,1X,A16,1X,A16,1X,A16)
100 FORMAT(1X,8000E17.6E3) !Printing format
	CHARACTER(20) filename1, filename2, filename3, filename4, filename5, filename6, filename7 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	r_1 = at_1 - e_1
	r_2 = at_2 - e_2

	alpha_1 = tpi_h*K_1/H_0
	tsigma_2 = (tpi_v)*(A_2/(A_2+(1-eta_1)*w_1*A_1))
	alpha_2 = pi_h*K_2/((1-eta_1)*w_1*A_1+A_2)
	tsigma_1 = (tpi_v)*((1-eta_1)*w_1*A_1/(A_2+(1-eta_1)*w_1*A_1))
	sigma_11 = (tpi_v*eta_1*w_1)/((1-w_1)+eta_1*w_1)
	sigma_12 = (tpi_v*(1-w_1)*w_2)/((1-w_1)+eta_1*w_1)
	talpha_11 = ((pi_h*eta_1*w_1)/(eta_1*w_1*A_1+(1-w_1)*A_1))*(K_1)
	talpha_12 = ((pi_h*(1-w_1)*w_2)/(eta_1*w_1*A_1+(1-w_1)*A_1))*(K_1)
	talpha_2 = (pi_h*(1-eta_1)*w_1/(A_2+(1-eta_1)*w_1*A_1))*(K_2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Dt_reso = T/N
    Dt = R1*Dt_reso
    t_h(1,1) = 0.0_dp
    t_h_apr(1,1) = 0.0_dp

    seed = getseed()                                           !
    call sgrnd(seed)
    
   	filename1 = 'Ih.dat'
	filename2 = 'Ia1.dat'
    filename3 = 'Iv1.dat'
	filename4 = 'Ia2.dat'
    filename5 = 'Iv2.dat'
    filename6 = 'estimadores_1.dat'
    filename7 = 'estimadores_2.dat' 
 
	OPEN(UNIT=45,FILE=filename1,position="APPEND")                                              !
	OPEN(UNIT=46,FILE=filename2,position="APPEND")
	OPEN(UNIT=47,FILE=filename3,position="APPEND")
	OPEN(UNIT=48,FILE=filename4,position="APPEND")
	OPEN(UNIT=49,FILE=filename5,position="APPEND")
	OPEN(UNIT=50,FILE=filename6,position="APPEND")
	OPEN(UNIT=51,FILE=filename7,position="APPEND")

    DO k1=1,L+1
       esp_x(k1,:) = 0.0_dp
       var_x(k1,:) = 0.0_dp
    END DO

    DO k2=1,L
       t_h(k2+1,1) = k2*Dt      
    END DO
    
	DO i=1,N1

    	j1 = 1

		X_temp = reshape((/ 0.1, 0.25, 0.1, 0.25, 0.1 /), (/ 5,1 /))
		X_exact_temp = reshape((/ 0.1, 0.25, 0.1, 0.25, 0.1 /), (/ 5,1 /))

        DO k=1,5
           solu_exact(1,k) = X_exact_temp(k,1)

           esp_x(1,k) = esp_x(1,k) + X_temp(k,1)
           var_x(1,k) = var_x(1,k) + X_exact_temp(k,1)**2                   
        END DO
        
        X1(1,i) = X_temp(1,1)
        X2(1,i) = X_temp(2,1)
        X3(1,i) = X_temp(3,1)
        X4(1,i) = X_temp(4,1)
        X5(1,i) = X_temp(5,1)

		DO j=1,L
        
        	DO i1=1,3
            	IF (j==1) THEN
                	dW(:,1) = 0.0d0
                	DO i2=2,R1
                    	dW(i1,i2) = sqrt(Dt_reso)*gaussrnd()
                    END DO
                ELSE
                    DO i2=1,R1
                    	dW(i1,i2) = sqrt(Dt_reso)*gaussrnd()
                    END DO
          		END IF
            END DO

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Wdisc = reshape((/ sum(dW(1,1:R1)), sum(dW(2,1:R1)), sum(dW(3,1:R1)) /),&
                    &(/ 3,1 /))

            F = reshape((/ (z_1*alpha_1*X_temp(3,1)*(1-X_temp(1,1))-mu_h*X_temp(1,1)),&
            	&(((talpha_11*tz_1+talpha_12*tz_2)*X_temp(3,1)+talpha_2*tz_1*X_temp(5,1))*&
                &(1-X_temp(2,1))-mu_1*X_temp(2,1)), ((z_1*pi_v*X_temp(1,1)+(tz_1*sigma_11+&
                &tz_2*sigma_12)*X_temp(2,1))*(1-X_temp(3,1))-(e_1+r_1)*X_temp(3,1)), &
            	&(tz_2*alpha_2*(1-X_temp(4,1))*X_temp(5,1)-mu_2*X_temp(4,1)), ((tz_1*tsigma_1*&
                &X_temp(2,1)+tz_2*tsigma_2*X_temp(4,1))*(1-X_temp(5,1))-(e_2+r_2)*X_temp(5,1)) /), (/5,1/)) 

            G = reshape((/ beta1*alpha_1*X_temp(3,1)*(1-X_temp(1,1)), 0.0d0, beta1*pi_v*X_temp(1,1)*&
            	&(1-X_temp(3,1)), 0.0d0, 0.0d0,&
                &0.0d0, beta2*(talpha_11*X_temp(3,1)+talpha_2*X_temp(5,1))*(1-X_temp(2,1)), beta2*sigma_11*&
                &X_temp(2,1)*(1-X_temp(3,1)), 0.0d0, beta2*tsigma_1*X_temp(2,1)*(1-X_temp(5,1)),&
                &0.0d0, talpha_12*X_temp(3,1)*beta3*(1-X_temp(2,1)), beta3*sigma_12*X_temp(2,1)*(1-X_temp(3,1)),&
                & beta3*alpha_2*X_temp(5,1)*(1-X_temp(4,1)), beta3*tsigma_2*X_temp(4,1)*(1-X_temp(5,1)) /),&
                & (/5,3/))

            X = X_temp + F*Dt + matmul(G,Wdisc)

            X_temp = X

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            F1 = reshape((/ (z_1*alpha_1*X_exact_temp(3,1)*(1-X_exact_temp(1,1))-mu_h*X_exact_temp(1,1)),&
            	&(((talpha_11*tz_1+talpha_12*tz_2)*X_exact_temp(3,1)+talpha_2*tz_1*X_exact_temp(5,1))*&
                &(1-X_exact_temp(2,1))-mu_1*X_exact_temp(2,1)), ((z_1*pi_v*X_exact_temp(1,1)+(tz_1*sigma_11+&
                &tz_2*sigma_12)*X_exact_temp(2,1))*(1-X_exact_temp(3,1))-(e_1+r_1)*X_exact_temp(3,1)), &
            	&(tz_2*alpha_2*(1-X_exact_temp(4,1))*X_exact_temp(5,1)-mu_2*X_exact_temp(4,1)), ((tz_1*tsigma_1*&
                &X_exact_temp(2,1)+tz_2*tsigma_2*X_exact_temp(4,1))*(1-X_exact_temp(5,1))-(e_2+r_2)*&
                &X_exact_temp(5,1)) /), (/5,1/)) 

            X_exact = X_exact_temp + F1*Dt

            X_exact_temp = X_exact

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
           
            DO k=1,5
              solu_exact(j+1,k) = X_exact(k,1)

           	  esp_x(j+1,k) = esp_x(j+1,k) + X_temp(k,1)
              var_x(j+1,k) = var_x(j+1,k) + X_temp(k,1)**2                 
            END DO

            IF (j==(100*j1)+1) THEN
              
                X1(j1+1,i) = X(1,1)
                X2(j1+1,i) = X(2,1)
	          	X3(j1+1,i) = X(3,1)
        	  	X4(j1+1,i) = X(4,1)
        	  	X5(j1+1,i) = X(5,1)

                IF (i==1) THEN
                    t_h_apr(j1+1,1) = t_h(j,1)
                END IF

                j1 = j1+1

            END IF

            IF (j==L) THEN
              
                X1(j1+1,i) = X(1,1)
                X2(j1+1,i) = X(2,1)
	          	X3(j1+1,i) = X(3,1)
        	  	X4(j1+1,i) = X(4,1)
        	  	X5(j1+1,i) = X(5,1)

                IF (i==1) THEN
                    t_h_apr(j1+1,1) = t_h(j+1,1)
                END IF

            END IF
            
                      
        END DO

    END DO

    esp_x = esp_x/N1

    DO i=1,(L/100)+1      
      write(45,100) t_h_apr(i,1), (X1(i,j),j=1,N1)
      write(46,100) t_h_apr(i,1), (X2(i,j),j=1,N1)
      write(47,100) t_h_apr(i,1), (X3(i,j),j=1,N1)
      write(48,100) t_h_apr(i,1), (X4(i,j),j=1,N1)
      write(49,100) t_h_apr(i,1), (X5(i,j),j=1,N1)
    END DO    

	WRITE(50,9) 'Tiempo','Humanos inf.','Esp.','Var.','Animales dom. inf.','Esp.','Var.','Vect. dom inf.',&
    		   &'Esp.','Var.'

	WRITE(50,11) '===============','===============','===============','===============','===============',&
    		    &'===============','===============','===============','===============','==============='

	WRITE(51,10) 'Tiempo','Animales silv. inf.','Esp.','Var.','Vect. silv. inf.','Esp.','Var.'

	WRITE(51,11) '===============','===============','===============','===============','===============',&
    			&'===============','==============='

    DO i=1,L+1
      DO j=1,5
         var_x(i,j) = (var_x(i,j)/N1) - esp_x(i,j)**2
      END DO

      write(50,100) t_h(i,1), solu_exact(i,1), esp_x(i,1), var_x(i,1), solu_exact(i,2), esp_x(i,2), var_x(i,2),&
       			   &solu_exact(i,3), esp_x(i,3), var_x(i,3)

      write(51,100) t_h(i,1), solu_exact(i,4), esp_x(i,4), var_x(i,4), solu_exact(i,5), esp_x(i,5), var_x(i,5)                   
      
    END DO

    close(45)
    close(46)
	close(47)
    close(48)
    close(49)
    close(50)
    close(51)
    
end program sde
