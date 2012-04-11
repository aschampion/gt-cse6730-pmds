	SUBROUTINE Neighbour
	
		USE Globals
		IMPLICIT NONE
		
		INTEGER :: numcell,celli,cellj,cellID,i,j,oi,oj,ocelli,ocellj
		DOUBLE PRECISION:: Lc,xr,yr,rr
		
		!Bin the atoms in the cells
		Lc = MAXVAL(Rcut)+Rskin
		Lc = Box/INT(Box/Lc)
		numcell = Box/Lc
		DO i=1,numcell**2
			hoc(i) = 0
		END DO
		DO i=1,Natom
			celli = INT(Xx(i)/Lc)
			cellj = INT(Yy(i)/Lc)
			cellID= cellj*numcell+celli
			LL(i) = hoc(cellID)
			hoc(cellID) = i
		END DO
		
		!Construction of Verlet List
		DO i=1,Natom
			nlist(i) = 0
			celli = INT(Xx(i)/Lc)
			cellj = INT(Yy(i)/Lc)
			DO oi=-1,1
				DO oj=-1,1
					ocelli = celli+oi
					ocellj = cellj+oj
					IF(ocelli .GT. numcell) THEN
						ocelli = ocelli-numcell
					ELSE IF(ocelli .LE. 0) THEN
						ocelli = ocelli+numcell
					ENDIF
					IF(ocellj .GT. numcell) THEN
						ocellj = ocellj-numcell
					ELSE IF(ocelli .LE. 0) THEN
						ocellj = ocellj+numcell
					ENDIF
					cellID = ocellj*numcell+ocelli
					j = hoc(cellID)
					DO WHILE(LL(j) .NE. 0)
						IF(LL(j) .GT. i) THEN
							xr = Xx(i) - Xx(LL(j))
							yr = Yy(i) - Yy(LL(j))
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
							IF(rr < (MAXVAL(Rcut) + Rskin)) THEN
								nlist(i) = nlist(i)+1
								nlist(LL(j)) = nlist(LL(j))+1
								list(i,nlist(i)) = LL(j)
								list(LL(j),nlist(LL(j))) = i
							END IF
						END IF
						j = LL(j)
					END DO
				END DO
			END DO
		END DO
		
	
	END SUBROUTINE Neighbour
