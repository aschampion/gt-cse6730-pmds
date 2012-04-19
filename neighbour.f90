	SUBROUTINE Neighbour
	
		USE Globals
		IMPLICIT NONE
		
		INTEGER :: numcell,celli,cellj,cellID,i,j,oi,oj,ocelli,ocellj,s,k
		DOUBLE PRECISION:: Lc,xr,yr,rr
		
		!Bin the atoms in the cells
		Lc = MAXVAL(Rcut)
		Lc = MAX(Lc,Rcut_bond,Rcut_soft)+Rskin
		Lc = Box/INT(Box/Lc)
		numcell = Box/Lc
		DO i=1,numcell
			DO j=1,numcell
				hoc(i,j) = 0
			END DO
		END DO
	
		DO i=1,Natom
			celli = INT(Xx(i)/Lc)+1
			cellj = INT(Yy(i)/Lc)+1
			LL(i) = hoc(celli,cellj)
			hoc(celli,cellj) = i
			nlist(i) = 0
		END DO
		
		!Construction of Verlet List
		DO i=1,Natom
			celli = INT(Xx(i)/Lc)+1
			cellj = INT(Yy(i)/Lc)+1
			DO oi=-1,1						!offset to neighbouring cell in x direction
				DO oj=-1,1					!offset to neighbouring cell in y direction
					ocelli = celli+oi
					ocellj = cellj+oj
					IF(ocelli .GT. numcell) THEN
						ocelli = ocelli-numcell
					ELSE IF(ocelli .LE. 0) THEN
						ocelli = ocelli+numcell
					ENDIF
					IF(ocellj .GT. numcell) THEN
						ocellj = ocellj-numcell
					ELSE IF(ocellj .LE. 0) THEN
						ocellj = ocellj+numcell
					ENDIF
					j = hoc(ocelli,ocellj)
					DO WHILE(j .NE. 0)
						IF(j .GT. i) THEN
							xr = Xx(i) - Xx(j)
							yr = Yy(i) - Yy(j)
							IF(xr .GT. Box/2.0D0) THEN
								xr = xr-Box
							ELSE IF(xr .LE. -Box/2.0D0) THEN
								xr = xr+Box
							ENDIF
							IF(yr .GT. Box/2.0D0) THEN
								yr = yr-Box
							ELSE IF(yr .LE. -Box/2.0D0) THEN
								yr = yr+Box
							ENDIF
							rr = SQRT(xr**2+yr**2)
							IF(rr .LT. Lc) THEN
								nlist(i) = nlist(i)+1
								nlist(j) = nlist(j)+1
								list(i,nlist(i)) = j
								list(j,nlist(j)) = i
							END IF
						END IF
						j = LL(j)
					END DO
				END DO
			END DO
		END DO
		
	END SUBROUTINE Neighbour
