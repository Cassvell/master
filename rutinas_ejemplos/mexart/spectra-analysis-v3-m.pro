;PROGRAM TO OBTAIN POWER SPECTRA OF IPS DATA AND CALCULATE m (Mejia-Ambriz, Gonzalez-Esparza; june 2016)
;(all lines indicated with asterisks (*) can be changed for a particular observation)

;==============================================================================
; (I) SELECT DATA TO BE ANALYZED AND SET UP SAMPLING CHARACTERISTICS
;==============================================================================
;read data in the indicated path, normalize the off-set flux to 0 and plot at the screen
                 
file = DIALOG_PICKFILE(/READ, GET_PATH=dir)
file_name = FILE_BASENAME(file) 

READCOL, file, time, flux,FORMAT='D,D' 
IF (MIN(flux) LT 0) THEN flux=flux+ABS(MIN(flux)) ELSE flux=flux-MIN(flux)

n_int=12  ;number of subintervals ;*******************************
secs=10.24; seconds in each subinterval (2^9 points) *********************
T = 0.02 ;MEXART and STEL's sampling (.02 sec or 50 Hz) ***********
n=LONG(secs/T) ; n is the number of points for each subinterval 
interval=n_int*n  ;number of points in all the interval;

PLOT,time,flux, TITLE=file_name
;By using the cursor over the plot choose the initial point of an interval, 
PRINT, "Click cursor at ON-SOURCE start time"    
      CURSOR, x, y, /DOWN
i0=0L
WHILE (time[i0] LT x) DO i0=i0+1
t1=time[i0]
t2=time[i0+interval-1]                              
OPLOT,[t1,t2],[y, y], LINESTYLE = 5, THICK = 3 ;Plot an horizontal line on the interval
 
time1=time[i0:i0+interval-1] ; 
flux1=flux[i0:i0+interval-1] ;the chosen interval                           

;=================================================================================
; (IIa) FOURIER TRANSFORMS USING THE SUBINTERVALS 
;smooth + interpolate low freq + average + interpolate low + smooth high freq
;=================================================================================

FT =FINDGEN((n/2) + 1)*(1.0/(n*T)) ;associates frequencies for n points

part_pws = DBLARR(n_int , n/2 + 1)
s_part_pws=DBLARR(n_int , n/2 + 1) ; (n_int columns) X (n/2 + 1 rows) 

FOR i=0, n_int-1 DO BEGIN
    part_flux = flux1[i*n : i*n + n -1]
    F_part=FFT(part_flux)
    part_pws[i,*]=(ABS(F_part[0:n/2]))^2
    s_part_pws[i,*]=SMOOTH(part_pws[i,*],5) ;smooth of +/- 0.2 Hz
ENDFOR

;interpolate at low freq subintervals (first three points) 
FOR i=0,n_int-1 DO BEGIN
    k=0
    WHILE (FT[k] LT 0.25) DO k=k+1
    XP=FT[k:k+2]
    YP=s_part_pws[i,k:k+2]
    fit_l=POLY_FIT(XP,YP,1)
    FOR j=k-1,0,-1 DO s_part_pws[i,j]=fit_l[0]+fit_l[1]*FT[j]
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
noise_level_l=DBLARR(n_int); array to estimate level of noise at low frequencies
part_pws_ln=DBLARR(n_int , n/2 + 1) ; array to allocate spectra without noise

PRINT,'Noise level at high freq for each spectrum '
kh=WHERE(FT GT 2.49 AND FT LT 10.01)
FOR i=0,(n_int)-1 DO BEGIN
   noise_level_h[i]=TOTAL(s_part_pws[i,kh])
   PRINT,STRTRIM(i), noise_level_h[i]
ENDFOR
min_noise_h=MIN(noise_level_h)

PRINT,'Noise level at low freq for each spectrum '
kl=WHERE(FT GT 0.2 AND FT LT 0.91) 
FOR i=0,(n_int)-1 DO BEGIN
   noise_level_l[i]=TOTAL(s_part_pws[i,kl])
   PRINT, STRTRIM(i), noise_level_l[i]
ENDFOR
avg_noise_l=AVG(noise_level_l)

nul0=0
FOR i=0, (n_int)-1 DO BEGIN
    IF (noise_level_h[i] lt 4.0*min_noise_h) AND (ABS( (avg_noise_l - noise_level_l[i])/avg_noise_l) LT 0.65 ) $  
   THEN part_pws_ln[i,*]=s_part_pws[i,*] ELSE nul0=nul0+1 ;noise level at low freq is at least .35 times the average and lower than 4 times the minimum at high freq
ENDFOR 
PRINT,'number of discarded spectra: ', nul0

avg_s_part_pws_ln=TOTAL(part_pws_ln,1)/(n_int - nul0)

;fit low freq
k=0
WHILE (FT[k] LT 0.25) DO k=k+1
XPL0=FT[k:k+2]
YPL0=avg_s_part_pws_ln[k:k+2]
FIT_LOW00=POLY_FIT(XPL0,YPL0,1)
FOR j=k-1,0,-1 DO avg_s_part_pws_ln[j]=FIT_LOW00[0]+FIT_LOW00[1]*FT[j];

;Smooth at high freq
knee=WHERE(ceil(10.0*FT) EQ 11)
sh_avg_s_part_pws_ln=[avg_s_part_pws_ln[0:knee],SMOOTH( avg_s_part_pws_ln[knee+1:n/2],3 )]
s_avg_s_part_pws_ln=SMOOTH( avg_s_part_pws_ln,3 )
;==========================================================================
;(IIc) CORRECT THE SPECTRA AT HIGH FREQUENCIES 
;==========================================================================
;first correction (time constant receiver)
rctc=0.047
tptau=(2.0*!PI*rctc)^2
tau= 1.0/(1+FT^2*tptau)

pws_tau_corr_ln=avg_s_part_pws_ln/tau

;second correction
high_i=WHERE(FT GT 2.0 AND FT LT 8)
poff=AVG(pws_tau_corr_ln[high_i])

pws_corr_ln = (pws_tau_corr_ln - poff)/Poff

s_pws_corr_ln=SMOOTH(pws_corr_ln,3)
;================================================================================
;(IId) Create a file with the name of the time-series
;file_mkdir, dir+'plots/'file_name
;================================================================================
;(IIe) Calculate scintillation index^2 with area under the curve of power spectra
;================================================================================
ips_i=WHERE(FT GT 0.25 AND FT LE 2.0)
m_index=SQRT(INT_TABULATED(FT[ips_i],pws_corr_ln[ips_i]))
PRINT,'m_index is',m_index

OPENW, oindex, dir+'plots/index-'+file_name+'.dat', /GET_LUN
PRINTF, oindex, 'm_index: ', m_index
FREE_LUN, oindex

;==============================================================================
;(III) PLOT TIMESERIES, CO-ADDED SPECTRA (-BAD SPECTRA), AND SPECTRA OF SUBINTERVALS.
;==============================================================================
 old_plot=!D.NAME; store current graphycal device in a variable
   !P.MULTI=[0,2,3] 
 SET_PLOT,'PS' 
   
   DEVICE, FILENAME=dir+'/plots/'+file_name+'-coad-bad' ;*********************
   
   PLOT,time,flux,XRANGE=[time[0],time[N_ELEMENTS(time)-1]],XSTYLE=1,TITLE='time series',$
   POS=[0.1,0.70,0.9,0.95]
   
   PLOT,time,flux,XRANGE=[time[0],time[N_ELEMENTS(time)-1]],XSTYLE=1.0,TITLE='time series',$
   POS=[0.1,0.70,0.9,0.95]
   
   OPLOT,[t1,t2],[y,y]   
   FOR l=0L,interval-1,n DO OPLOT, [time1[0+l],time1[0+l]], [y,y],PSYM=1
   
   PLOT, FT, s_part_pws[0,*], /XLOG,/YLOG, YTICKFORMAT='logticks_exp',$
   TITLE=STRTRIM(n_int,2)+' spectra taken each '+STRTRIM(secs,2)+' secs',XSTYLE=1, XRANGE=[0.1,10.0],$
   YRANGE=[MIN(s_part_pws[*,25:100]), MAX(s_part_pws[*,1:10])],YSTYLE=1
   FOR i=1, (n_int)-1 DO BEGIN
     OPLOT,FT,s_part_pws[i,*];  
   ENDFOR

   PLOT,FT,avg_s_part_pws,/XLOG,/YLOG, XRANGE=[0.1,10.0],XSTYLE=1, YTICKFORMAT='logticks_exp',$
   TITLE='sum'


   PLOT, FT, part_pws_ln[0,*],XSTYLE=1,XRANGE=[0.1,10.0],/XLOG,/YLOG, YTICKFORMAT='logticks_exp',$
   YRANGE=[MIN(s_part_pws[*,25:100]), MAX(s_part_pws[*,1:10])],YSTYLE=1, TITLE='drop bad spectra '
   FOR i=1, (n_int)-1 DO BEGIN
     OPLOT,FT,part_pws_ln[i,*];  
   ENDFOR

   PLOT, FT, avg_s_part_pws_ln,/XLOG, /YLOG,ytickformat='logticks_exp',xrange=[0.1,10.0],XSTYLE=1,$
   TITLE='sum'
   
DEVICE,/CLOSE_FILE
   
!P.MULTI=0   

;PLOT NO NOISY INDIVIDUAL POWER SPECTRA OF n POINTS EACH
FOR i=0, (n_int)-1 DO BEGIN
    SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-p'+STRTRIM(i,2)
       PLOT,FT,part_pws_ln[i,*],XRANGE=[0.1,10],/XLOG,/YLOG,TITLE='pws'+STRTRIM(i,2); 
       DEVICE, /CLOSE_FILE
ENDFOR

;PLOT FIRST CORRECTION TO THE SPECTRA
;SET_PLOT,'PS'
;       DEVICE,FILENAME=dir+'/plots/'+file_name+'-pws-tau-corr'
;       PLOT,FT,pws_tau_corr/pws_tau_corr[3],XRANGE=[0.1,10],/XLOG,/YLOG
;       ;YRANGE=[MIN(pws_tau_corr[3:103]), MAX(pws_tau_corr[3:103])], YSTYLE=1
;       DEVICE, /CLOSE_FILE
;PLOT SECOND CORRECTION
;SET_PLOT,'PS'
;       DEVICE,FILENAME=dir+'/plots/'+file_name+'-pws-corr'
;       PLOT,FT,pws_corr/pws_corr[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
;       YRANGE=[1e-4,1.5], YSTYLE=1
;       DEVICE, /CLOSE_FILE

;PLOT FIRST CORRECTION TO THE SPECTRA
SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-pws-ln-tau-corr'
       PLOT,FT,SMOOTH(pws_tau_corr_ln/pws_tau_corr_ln[3], 3),XRANGE=[0.1,10],/XLOG,/YLOG,$
       YSTYLE=1
       DEVICE, /CLOSE_FILE

;PLOT SECOND CORRECTION
;SET_PLOT,'PS'
;       DEVICE,FILENAME=dir+'/plots/'+file_name+'-pws-ln-corr'
;       PLOT,FT,pws_corr_ln/pws_corr_ln[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
;       YRANGE=[1.E-4,1.5], YSTYLE=1
;       DEVICE, /CLOSE_FILE

;overplot
SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-test'
       PLOT,FT, pws_tau_corr_ln/pws_tau_corr_ln[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
       YRANGE=[1.E-3,1.5], YSTYLE=1
       OPLOT,FT, pws_corr_ln/pws_corr_ln[3]
       OPLOT, FT, avg_s_part_pws_ln/avg_s_part_pws_ln[3]
       OPLOT, FT, s_pws_corr_ln/s_pws_corr_ln[3],LINESTYLE=2
;      DEVICE, /CLOSE_FILE

SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-smooth-high-freq'
       PLOT,FT,sh_avg_s_part_pws_ln/sh_avg_s_part_pws_ln[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
       YRANGE=[1.E-3,1.5], YSTYLE=1
       DEVICE, /CLOSE_FILE

SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-smooth'
       PLOT,FT,s_avg_s_part_pws_ln/s_avg_s_part_pws_ln[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
       YRANGE=[1.E-3,1.5], YSTYLE=1
       DEVICE, /CLOSE_FILE

;overplot
SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-overplot'
       PLOT,FT,sh_avg_s_part_pws_ln/sh_avg_s_part_pws_ln[3],XRANGE=[0.1,10],/XLOG,/YLOG,$
       YRANGE=[1.E-3,1.5], YSTYLE=1
       OPLOT,FT, s_avg_s_part_pws_ln/s_avg_s_part_pws_ln[3],LINESTYLE=2
       DEVICE, /CLOSE_FILE

SET_PLOT,'PS'
       DEVICE,FILENAME=dir+'/plots/'+file_name+'-pws-no-log'
       PLOT,FT,pws_corr_ln, XRANGE=[0.1,5],XSTYLE=1
DEVICE,/CLOSE
;====================================================================
;(IV) WRITE SPECTRA IN ASCII FILES
;====================================================================

;WRITE POINTS OF THE SMOOTHED SUM OF PWS
OPENW, o1, dir+'plots/pws-sum-'+file_name+'.dat', /GET_LUN
FOR i=0D, N_ELEMENTS(avg_s_part_pws)-1 DO PRINTF,$
    o1, FT[i], avg_s_part_pws[i];
FREE_LUN, o1

;WRITE POINTS OF THE SUM OF PWS -BAD SPECTRA
OPENW, o2, dir+'plots/pws-sum-bad-'+file_name+'.dat', /GET_LUN
FOR i=0d, N_ELEMENTS(avg_s_part_pws_ln)-1 DO PRINTF,$
    o2, FT[i], avg_s_part_pws_ln[i];
FREE_LUN, o2

;WRITE POINTS OF THE CORRECTED
;OPENW, o3, path+'plots/pws-sum-corr'+file_name+'.dat', /GET_LUN
;FOR i=0D, N_ELEMENTS(pws_corr)-1 DO PRINTF,$
;    o3, FT[i], pws_corr[i];
;FREE_LUN, o3

;WRITE POINTS OF THE CORRECTED (SMOOTHED) LESS NOISY SPECTRA
OPENW, o4, dir+'plots/s-pws-corr'+file_name+'.dat', /GET_LUN
FOR i=0D, N_ELEMENTS(s_pws_corr_ln)-1 DO PRINTF,$
    o4, FT[i], s_pws_corr_ln[i];
FREE_LUN, o4


;WRITE FIRST CORRECTION TO THE SPECTRA
p_test=FINDGEN(N_ELEMENTS(pws_tau_corr_ln))
p_test=SMOOTH(pws_tau_corr_ln,3)
OPENW, o5, dir+'plots/pws-sum-ln-tau'+file_name+'.dat', /GET_LUN
FOR i=0D, N_ELEMENTS(pws_tau_corr_ln)-1 DO PRINTF,$
    o5, FT[i], p_test[i];
FREE_LUN, o5


;WRITE SMOOTH AT HIGH FREQ
OPENW, o6, dir+'plots/s-h-'+file_name+'.dat', /GET_LUN
FOR i=0D,n/2-1 DO PRINTF,$
    o6, FT[i], sh_avg_s_part_pws_ln[i];
FREE_LUN, o6

;WRITE SMOOTH (3 POINTS) TO ALL THE SPECTRA
OPENW, o7, dir+'plots/s-'+file_name+'.dat', /GET_LUN
FOR i=0D, n/2 - 1 DO PRINTF,$
    o7, FT[i], s_avg_s_part_pws_ln[i];
FREE_LUN, o7



SET_PLOT, old_plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;,
END





