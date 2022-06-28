;PROGRAM TO OBTAIN POWER SPECTRA OF IPS DATA AND CALCULATE m (Mejía- Ambriz, Gonzales-Esparza; june 2016)
;(all lines indicated with asterisks (*) can be for a particular observation)


;================================================================================================
;(I) SELECT DATA TO BE ANALYZED AND SET UP SAMPLING CHARACTERISTICS
;================================================================================================
;read data in the indicated path, normalize the off-set flux to 0 and plot at the screen

file = DIALOG_PICKFILE(/READ, PATH = '~/Desktop/mexart-test/transits/', GET_PATH=dir)
file_name = FILE_BASENAME(file)

READCOL, file, time, flux, FORMAT='D,D'
IF (MIN(flux) LT 0) THEN flux=flux+ABS(MIN(flux)) ELSE flux=flux-MIN(flux)

n_int=12	;number of intervals ;**********************************
secs=10.24	;seconds in each subinterval(2^9 points) *************
T = 0.02	;MEXART and STEL's sampling (0.2 sec or 50 Hz) ************
n=LONG(secs/T)	;n is the number of poits for each subinterval ***********
interval=n_int*n	; nuber of points in allintervals;

PLOT,time,flux, TITLE=file_name
;By using the cursor over the choose the initial point of an interval,
PRINT, "Click cusor at ON-SOURCE start time"
	CURSOR, X, y, /DOWN
i0=0L
WHILE (time[i0] LT x) DO i0=i0+1
t1=time[i0]
t2=time[i0+interval-1]
OPLOT,[t1,t1],[y, y], LINESTYLE = 5, THICK = 3	;Plot an horzontal line on the interval

time1=time[i0:i0+interval-1];
flux1=flux[i0:i0+interval-1]; th chosen interval

;=========================================================================================
; (IIa) FOURIER TRANSFORMS USING THE SUBINTERVALS
;smooth + interpolate low freq + averge + interpolate low + smooth high freq
;=========================================================================================

FT = FINDGEN((n/2) + 1)*(1.0/(n*t)) ;associates frequencies for n points

parts_pws = DBLARR(n_int, n/2 + 1)
s_part_pws= DBLARR(n_int, n/2 + 1)	;(n_int columns) X (n/2 + 1 rows)

FOR i=0, n_int-1 DO BEGIN
    part_flux = flux1[i*n  :  i*n +  -1]
    F_part=FFT(part_flux)
    part_pws[i,*]=SMOOTH(part_pws[i,*],5)	;smooth of x/- 0.2 Hz
ENDFOR

;interpolate at low freq subintervals (first three points)
FOR i=0, n_int-1 DO BEGIN
    k=0
    WHILE (FT[k] LT 0.25) DO k=k+1
    XP=FT[k:k+2]
    YP=s_part_pws[i,k:k+2]
    fit_l=POLY_FIT(XP,YP,1)
    FOR j=k-1,0,-1 DO s_part_pws[i,j]=fit_l[1]*FT[j]
ENDFOR

avg_s_part_pws=( TOTAL(s_part_pws,1) )/n_int;

;interpolate at low freq the sum
k=0
WHILE (FT[k] LT 0.25) DO k=k+1
XP=FT[k:k+2]
YP=avg_s_part_pws[k:k+2]
Fit_l=POLY_FIT(XP,YP,1)
FOR j=k-1,0,-1 DO avg_s_part_pws[j]=fit_l[0]+fit_l[1]*FT[j]

;=======================================================================
;(IIb) QUIT BAD SPECTRA
;=======================================================================


noise_level_h=DBLARR(n_int); array to estimate level of noise at high frequencies
noise_level_l=DBLARR(n_int); array to estimate level of noise at high frec
























