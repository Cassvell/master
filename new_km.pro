;Name: new_km.pro
;	
;purpose:
;	this routine will call read a certain number of files containing 
;   magnetic data index to derive a new kmex index considering diono effects on Kp
;author:
;	Carlos Isaac Castellanos Velazco
;	Estudiante de Maestría en Ciencias de la Tierra
;	Instituto de Geofísica, Unidad Michoacan
;	UNAM
;	ccastellanos@igeofisica.unam.mx
;
;category:
;   data analysis
;
;calling sequence:
;   .r new_km
;   new_km, r_dst, Bsq, DOY, idate[yyyy,mm,dd], fdate[yyyy,mm,dd]
;parameters:
;
;
;dependencies:
;   Instituto de geofísica, UNAM
;
;input files
;   DST index, DH index, SQ_Baseline
;
;output files:
;   n_kmex index en archivos de un día con 3 horas de muestreo.
FUNCTION dst_data, initial
	ON_ERROR, 2
	COMPILE_OPT idl2, HIDDEN

	year = STRING(initial, format = '(I4)')
	;type_data = string(tp, format = '(A)')
		file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_D.dat'
		
        header = 25             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		
	    IF opened_files NE N_ELEMENTS(file) THEN BEGIN
	        file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_P.dat'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN BEGIN
    	        file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_Q.dat'
	            file = FILE_SEARCH(file_name, COUNT=opened_files)    	        
    	        IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
    	    ENDIF    	    	    
	    ENDIF

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		OPENR, LUN, file, /GET_LUN, ERROR=err
		READF, LUN, data, FORMAT = '(A)'
		CLOSE, LUN
		FREE_LUN, LUN
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
        DOY : 0, Dst: 0}

		r_dst = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], r_dst, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,8X,I3,X,I6)'
		
		RETURN, r_dst
END

FUNCTION DH_teo, date
	ON_ERROR, 2
	COMPILE_OPT idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = STRING(year, month, day, format = '(I4, I02, I02)')		
		file_name = '../rutidl/dH_teo/'+'teo_'+date+'.dst.early'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		OPENR, LUN, file, /GET_LUN, ERROR=err
		READF, LUN, data, FORMAT = '(A)'
		CLOSE, LUN
		FREE_LUN, LUN
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        DStruct = {hora : 0, D_stdesv : 0., D : 0., H_stdesv : 0., H : 0., $
        Z_stdesv : 0., Z : 0., N_stdesv : 0., N : 0., F_stdesv : 0., F : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(I2, F10, F8, F10, F10, F10, F10, F10, F10, F10, F10)'		
		RETURN, teo_mag		
END

FUNCTION kp_data, in
	ON_ERROR, 2
	COMPILE_OPT idl2, HIDDEN

        header = 36             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        year = string(in, format = '(I4)')
		file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_D.dat'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)

	    IF opened_files NE N_ELEMENTS(file) THEN BEGIN
	        file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_P.dat'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN BEGIN
    	        file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_Q.dat'
	            file = FILE_SEARCH(file_name, COUNT=opened_files)    	        
    	        IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
    	    ENDIF    	    	    
	    ENDIF     
        
		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		OPENR, LUN, file, /GET_LUN, ERROR=err
		READF, LUN, data, FORMAT = '(A)'
		CLOSE, LUN
		FREE_LUN, LUN
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
                      DOY : 0, Kp: 0, Kp_str: '', Ap: 0}

		resulting_data0 = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], resulting_data0, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,8X,I3,X,I1,A1,X,I3)'
		
		indexes_07 = WHERE(resulting_data0[*].Kp_str EQ '-')
		indexes_03 = WHERE(resulting_data0[*].Kp_str EQ '+')
		indexes_00 = WHERE(resulting_data0[*].Kp_str EQ 'o')
		
		Kp_tmp = resulting_data0[*].Kp*10
		
		Kp_tmp[indexes_07] = Kp_tmp[indexes_07]-3
		Kp_tmp[indexes_03] = Kp_tmp[indexes_03]+3


        DataStruct2 = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
                      DOY : 0, Kp: 0, Ap: 0}

		r_kp = REPLICATE(DataStruct2, number_of_lines-header)	

		r_kp[*].year   = resulting_data0[*].year
		r_kp[*].month  = resulting_data0[*].month
		r_kp[*].day    = resulting_data0[*].day
		r_kp[*].hour   = resulting_data0[*].hour
		r_kp[*].minute = resulting_data0[*].minute
		r_kp[*].DOY    = resulting_data0[*].DOY
		r_kp[*].Kp     = Kp_tmp[*]
		r_kp[*].Ap     = resulting_data0[*].AP				
		RETURN, r_kp
END    

FUNCTION baseline_sq, date
	ON_ERROR, 2
	COMPILE_OPT idl2, HIDDEN
	
	yr	= date[0]
	mh	= date[1]

        date = string(yr, mh, format = '(I4, "-", I02)')
        header=0

		file_name = '../rutidl/output/Bsq_baselines/'+'Bsq_'+date+'.txt'
	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		OPENR, LUN, file, /GET_LUN, ERROR=err
		READF, LUN, data, FORMAT = '(A)'
		CLOSE, LUN
		FREE_LUN, LUN
		
        DStruct = {year : 0, month : 0, day : 0, hour : 0, doy : 0, Bsq : 0.}                                   

		B_sq = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], B_sq, $
	 format='(I4,X, I02,X, I02, 2X, I02, 2X, I03, 2X, F08.4)' 		
		RETURN, B_sq
END


FUNCTION Date2DOY, idate
;	Check data type of input set ascII flag and convert to yy,mm,dd:
	info = SIZE(idate)
	IF (info(0) EQ 0) THEN BEGIN
	  scalar = 1				;scalar flag set
	ENDIF ELSE BEGIN
	  scalar = 0				;vector input
	ENDELSE

	IF (info(info(0) + 1) EQ 7) THEN BEGIN
	  ascII = 1				;ascII input flag set
	  yy = FIX(STRMID(idate,0,2))		;extract year
	  mm = FIX(STRMID(idate,2,2))		;extract month
	  dd = FIX(STRMID(idate,4,2))		;extract day
	ENDIF ELSE BEGIN			;should be a longWord
	  ascII = 0				;non-ascII input
	  sdate = STRTRIM(STRING(idate),2)	;convert to string 
	  yy = FIX(STRMID(sdate,0,2))		;extract year
	  mm = FIX(STRMID(sdate,2,2))		;extract month
	  dd = FIX(STRMID(sdate,4,2))		;extract day
	ENDELSE

;	Check for leap year and compute DOY:
;       	      J   F   M   A   M   J   J   A   S   O   N   D
	imonth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	IF (scalar) THEN BEGIN			;scalar input
	  IF ((yy MOD 4) EQ 0) THEN BEGIN	;leap year
	    imonth(2) = 29			;set feb
	  ENDIF
	  DOY = FIX( TOTAL(imonth(0:mm-1)) ) + dd
	ENDIF ELSE BEGIN
	  DOY = dd				;set correct len on vector
	  leapYrs = WHERE( (yy MOD 4) EQ 0)	;index of leap years
	  nonLeap = WHERE( (yy MOD 4) NE 0)	;index of non-leap years
	  IF (nonLeap(0) NE -1) THEN BEGIN
	    FOR i=0, N_elements(nonLeap)-1 DO BEGIN
	      DOY(nonLeap(i)) = FIX( TOTAL(imonth(0:mm(nonLeap(i))-1)) ) + $
				dd(nonLeap(i))
	    ENDFOR
	  ENDIF
	  IF (leapYrs(0) NE -1) THEN BEGIN
	    imonth(2) = 29			;set feb
	    FOR i =0, N_elements(leapYrs)-1 DO BEGIN
	      DOY(leapYrs(i)) = FIX( TOTAL(imonth(0:mm(leapYrs(i))-1)) ) + $
				dd(leapYrs(i))
	    ENDFOR
	  ENDIF
	ENDELSE

	IF (N_PARAMS() EQ 3) THEN BEGIN         ;pass year back to caller
          IF (ascII) THEN BEGIN
	    DOY = STRTRIM( STRING(DOY), 2)	;convert to string	    
	    yr = STRTRIM( STRING(yy), 2)	;convert to string	  
	  ENDIF ELSE BEGIN
	    yr = yy	
	  ENDELSE			
	ENDIF ELSE BEGIN			;pass DOY only
	  IF (ascII) THEN BEGIN
	    DOY = STRTRIM( STRING(DOY), 2)	;convert to string
	  ENDIF
	ENDELSE  
	RETURN, DOY
	END
	
	
PRO new_km, r_dst, r_kp, B_sq, DOY, date_i, date_f
	ON_ERROR, 2
	COMPILE_OPT idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
;###############################################################################
idate0 = STRING(yr_i, mh_i, format='(I4,I02)')
TGM_n = idate0
CASE TGM_n of
    '200311' : TGM_n = 1
    '200411' : TGM_n = 2
    '200505' : TGM_n = 3
    '201503' : TGM_n = 4
    '201705' : TGM_n = 5
    '201709' : TGM_n = 6
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
    d_dst = dst_data(yr_i)     
    i_dst   = d_dst.Dst
    hour    = d_dst.hour
    dst_doy = d_dst.DOY          
;###############################################################################   
    d_dst = dst_data(yr_i)
    t = N_ELEMENTS(d_dst.year)    
    i_dst = d_dst.Dst       
;###############################################################################    
    d_kp = kp_data(yr_i)     
    i_kp = d_kp.Kp       
;###############################################################################      
    year = d_dst.year
    tiempo = TIMEGEN(t, START=julday(d_dst.month[0], d_dst.day[0],  $
                     d_dst.year[0], d_dst.hour[0]), UNITS='Hours')  
                                        
        iyear = STRMID(STRING(yr_i, format='(I4)'),2,2)
        fyear = STRMID(STRING(yr_f, format='(I4)'),2,2)
                
        idoy      = Date2DOY(STRING(iyear, mh_i, dy_i,format = '(I02,I02,I02)'))
        fdoy      = Date2DOY(STRING(fyear, mh_f, dy_f,format = '(I02,I02,I02)'))         
            
    time_w  = tiempo[idoy:fdoy]
    tw      = N_ELEMENTS(time_w)
    tot_days= FINDGEN(tw*24)/24.0
    k_days  = FINDGEN(tw*8)/8.0    
    dst     = i_dst[(idoy*24)-24:fdoy*24-1]    
    kp      =  i_kp[(idoy*8)-8:fdoy*8-1]
    kp      = kp/10.0    
    Date    = STRING(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')    
;###############################################################################
; define DH variables
;###############################################################################
       file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
        string_date        = STRARR(file_number)                       
        data_file_name_dh  = STRARR(file_number)        
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                
                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'		                                    
        ENDFOR

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))
        
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF
                        
        H    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_dh)-1 DO BEGIN
                IF exist_data_file_dh[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_dh = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_dh.H[*]
                                                                                                                       
                ENDIF ELSE BEGIN
                        H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
        
        i_nan1 = WHERE(H EQ 999999.0, ncount)
        i_nan2 = WHERE(H GT 100.0, n2count)
        
        prcent_nan = FLOAT(ncount+n2count)*100.0
        PRINT,'porcentaje de valores NaN:', prcent_nan/N_ELEMENTS(H),'%'
        
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 999999.0 THEN BEGIN
                H[WHERE(H[*] EQ 999999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR
        
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] GE 100.0 THEN BEGIN
                H[WHERE(H[*] GE 100.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR                
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo        
    H_tmp   = H
    H_exist = WHERE(finite(H_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
    ; interpolate at the locations of the bad data using the good data    
    IF nbaddata GT 0 THEN H_tmp[baddata] = interpol(H_tmp[H_exist], H_exist, baddata)
    H = H_tmp        
;############################################################################### 
; define diurnal baseline
;###############################################################################  
    sqline = baseline_sq([yr_i, mh_i])
    
    sq_doy  = sqline.doy
    Bsq_ln   = sqline.Bsq

    tmp_doy = dst_doy[(idoy*24)-24:fdoy*24-1]
    Bsq     = FLTARR(n_elements(Bsq_ln))
    
    tmp_doy = tmp_doy[uniq(tmp_doy, sort(tmp_doy))]
    
    FOR i=0, N_ELEMENTS(tmp_doy)-1 DO BEGIN
            FOR j=0, N_ELEMENTS(sq_doy)-1 DO BEGIN
                IF tmp_doy[i] EQ sq_doy[j] THEN BEGIN
               Bsq[j]    = Bsq_ln[j]                            
                ENDIF 
            ENDFOR
    ENDFOR
b_sq = WHERE(Bsq EQ 0, zcount, complement=val, ncomplement=valcount)
Bsq      = Bsq[val]
;###############################################################################
; define Diono
;###############################################################################  
    mlat        = 28.06*!PI
    ld          = COS(mlat/180)
    p_a         = dst*ld
    baseline    = Bsq + p_a             
    diono       = H-baseline
    time = 3600.0
    fn      = FLOAT(1.0/(2.0*time)) ; frecuencia de Nyquist    
;###############################################################################
; define pass band frequencies
;###############################################################################  
passband_l = idate0
CASE passband_l of
    '200311' : passband_l = 1.2e-5
    '200411' : passband_l = 1e-5
    '200505' : passband_l = 1.2e-5
    '201503' : passband_l = 1.14e-5
    '201705' : passband_l = 0.9e-5
    '201709' : passband_l = 1.2e-5
    ELSE: PRINT, 'fuera de rango'
ENDCASE  

passband_u = idate0
CASE passband_u of
    '200311' : passband_u = 3.2e-5
    '200411' : passband_u = 2.5e-5
    '200505' : passband_u = 2.8e-5
    '201503' : passband_u = 2.8e-5
    '201705' : passband_u = 3.2e-5
    '201709' : passband_u = 3.4e-5
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
; define high band frequencies
;###############################################################################
highpass_l = idate0
CASE highpass_l of
    '200311' : highpass_l = 7e-5
    '200411' : highpass_l = 7.5e-5
    '200505' : highpass_l = 7.3e-5
    '201503' : highpass_l = 7.5e-5
    '201705' : highpass_l = 7.5e-5
    '201709' : highpass_l = 7.e-5
    ELSE: PRINT, 'fuera de rango'
ENDCASE
;###############################################################################
; define filtering
;############################################################################### 
coeff_ddyn    = DIGITAL_FILTER(passband_l/fn, passband_u/fn, 50, 18)
coeff_dp2   = DIGITAL_FILTER(highpass_l/fn, 1.0, 50, 4)
;###############################################################################
; define disturbing effects
;############################################################################### 
ddyn        = CONVOL(diono, coeff_ddyn, /edge_wrap)
dp2         = CONVOL(diono, coeff_dp2, /edge_wrap)        
;###############################################################################
;Computing new DH
;###############################################################################
new_H       = baseline+ddyn+dp2
;Pasar de horario a tri-horario
print, '########################################################################'
;###############################################################################
;a partir de Kp, obtenemos dos valores DH superior e inferior
dH_tab      = [4.8, 5.5, 6.5, 7.3, 8.3, 9.8, 11.1, 12.6, 14.9, 16.9, 19.2,  $
               22.7, 25.7, 29.2, 34.5, 39.1, 44.4, 52.4, 59.5, 67.4, 79.7, 90.4, $
               102.5, 121.2, 137.4, 155.8, 184.2, 208.9, 500.0]

dh1         = FLTARR(N_ELEMENTS(kp))   
dh2         = FLTARR(N_ELEMENTS(kp))              
FOR i=0, N_ELEMENTS(kp)-1 DO BEGIN
    IF kp[i] EQ 0.0 THEN dh1[i] = 0.0     
    IF kp[i] EQ 0.3 THEN dh1[i] = dH_tab[1] 
    IF kp[i] EQ 0.7 THEN dh1[i] = dH_tab[2]
    IF kp[i] EQ 1.0 THEN dh1[i] = dH_tab[3]   
    IF kp[i] EQ 1.3 THEN dh1[i] = dH_tab[4] 
    IF kp[i] EQ 1.7 THEN dh1[i] = dH_tab[5]  
    IF kp[i] EQ 2.0 THEN dh1[i] = dH_tab[6]  
    IF kp[i] EQ 2.3 THEN dh1[i] = dH_tab[7]   
    IF kp[i] EQ 2.7 THEN dh1[i] = dH_tab[8]  
    IF kp[i] EQ 3.0 THEN dh1[i] = dH_tab[9]  
    IF kp[i] EQ 3.3 THEN dh1[i] = dH_tab[10] 
    IF kp[i] EQ 3.7 THEN dh1[i] = dH_tab[11] 
    IF kp[i] EQ 4.0 THEN dh1[i] = dH_tab[12]  
    IF kp[i] EQ 4.3 THEN dh1[i] = dH_tab[13] 
    IF kp[i] EQ 4.7 THEN dh1[i] = dH_tab[14] 
    IF kp[i] EQ 5.0 THEN dh1[i] = dH_tab[15] 
    IF kp[i] EQ 5.3 THEN dh1[i] = dH_tab[16] 
    IF kp[i] EQ 5.7 THEN dh1[i] = dH_tab[17] 
    IF kp[i] EQ 6.0 THEN dh1[i] = dH_tab[18] 
    IF kp[i] EQ 6.3 THEN dh1[i] = dH_tab[19]
    IF kp[i] EQ 6.7 THEN dh1[i] = dH_tab[20] 
    IF kp[i] EQ 7.0 THEN dh1[i] = dH_tab[21] 
    IF kp[i] EQ 7.3 THEN dh1[i] = dH_tab[22] 
    IF kp[i] EQ 7.7 THEN dh1[i] = dH_tab[23] 
    IF kp[i] EQ 8.0 THEN dh1[i] = dH_tab[24] 
    IF kp[i] EQ 8.3 THEN dh1[i] = dH_tab[25] 
    IF kp[i] EQ 8.7 THEN dh1[i] = dH_tab[26] 
    IF kp[i] EQ 9.0 THEN dh1[i] = dH_tab[27]                                                                  
;###############################################################################
    IF kp[i] EQ 0.0 THEN dh2[i] = dH_tab[1]    
    IF kp[i] EQ 0.3 THEN dh2[i] = dH_tab[2] 
    IF kp[i] EQ 0.7 THEN dh2[i] = dH_tab[3]
    IF kp[i] EQ 1.0 THEN dh2[i] = dH_tab[4]    
    IF kp[i] EQ 1.3 THEN dh2[i] = dH_tab[5] 
    IF kp[i] EQ 1.7 THEN dh2[i] = dH_tab[6]   
    IF kp[i] EQ 2.0 THEN dh2[i] = dH_tab[7] 
    IF kp[i] EQ 2.3 THEN dh2[i] = dH_tab[8]    
    IF kp[i] EQ 2.7 THEN dh2[i] = dH_tab[9] 
    IF kp[i] EQ 3.0 THEN dh2[i] = dH_tab[10]  
    IF kp[i] EQ 3.3 THEN dh2[i] = dH_tab[11]
    IF kp[i] EQ 3.7 THEN dh2[i] = dH_tab[12]   
    IF kp[i] EQ 4.0 THEN dh2[i] = dH_tab[13]
    IF kp[i] EQ 4.3 THEN dh2[i] = dH_tab[14] 
    IF kp[i] EQ 4.7 THEN dh2[i] = dH_tab[15] 
    IF kp[i] EQ 5.0 THEN dh2[i] = dH_tab[16] 
    IF kp[i] EQ 5.3 THEN dh2[i] = dH_tab[17] 
    IF kp[i] EQ 5.7 THEN dh2[i] = dH_tab[18]  
    IF kp[i] EQ 6.0 THEN dh2[i] = dH_tab[19]
    IF kp[i] EQ 6.3 THEN dh2[i] = dH_tab[20] 
    IF kp[i] EQ 6.7 THEN dh2[i] = dH_tab[21] 
    IF kp[i] EQ 7.0 THEN dh2[i] = dH_tab[22] 
    IF kp[i] EQ 7.3 THEN dh2[i] = dH_tab[23] 
    IF kp[i] EQ 7.7 THEN dh2[i] = dH_tab[24] 
    IF kp[i] EQ 8.0 THEN dh2[i] = dH_tab[25] 
    IF kp[i] EQ 8.3 THEN dh2[i] = dH_tab[26] 
    IF kp[i] EQ 8.7 THEN dh2[i] = dH_tab[27]  
    IF kp[i] EQ 9.0 THEN dh2[i] = dH_tab[28]                                                                            
ENDFOR
;###############################################################################
;Computing new kmex
;###############################################################################
H_3h    = FLTARR(N_ELEMENTS(kp))
n_km1   = FLTARR(N_ELEMENTS(kp))
n_km2   = FLTARR(N_ELEMENTS(kp))
err_sup = FLTARR(N_ELEMENTS(kp))
err_inf = FLTARR(N_ELEMENTS(kp))
FOR i=0, N_ELEMENTS(kp)-1 DO BEGIN
    H_3h[i]    = ABS(MAX(new_H[(i*3):(i+1)*3-1])-MIN(new_H[(i*3):(i+1)*3-1]))
    err_inf[i] = dh1[i]-H_3h[i]
    err_sup[i] = dh2[i]+H_3h[i] 
    
    IF err_sup[i] LE dH_tab[1]  THEN n_km1[i] = 0 
    IF err_sup[i] GT dH_tab[1]  AND err_sup[i] LE dh_tab[2]  THEN n_km2[i] = 3 
    IF err_sup[i] GT dH_tab[2]  AND err_sup[i] LE dh_tab[3]  THEN n_km2[i] = 7
    IF err_sup[i] GT dH_tab[3]  AND err_sup[i] LE dH_tab[4]  THEN n_km2[i] = 10 
    IF err_sup[i] GT dH_tab[4]  AND err_sup[i] LE dH_tab[5]  THEN n_km2[i] = 13
    IF err_sup[i] GT dH_tab[5]  AND err_sup[i] LE dH_tab[6]  THEN n_km2[i] = 17
    IF err_sup[i] GT dH_tab[6]  AND err_sup[i] LE dH_tab[7]  THEN n_km2[i] = 20
    IF err_sup[i] GT dH_tab[7]  AND err_sup[i] LE dH_tab[8]  THEN n_km2[i] = 23
    IF err_sup[i] GT dH_tab[8]  AND err_sup[i] LE dH_tab[9]  THEN n_km2[i] = 27
    IF err_sup[i] GT dH_tab[9]  AND err_sup[i] LE dH_tab[10] THEN n_km2[i] = 30
    IF err_sup[i] GT dH_tab[10] AND err_sup[i] LE dH_tab[11] THEN n_km2[i] = 33
    IF err_sup[i] GT dH_tab[11] AND err_sup[i] LE dH_tab[12] THEN n_km2[i] = 37
    IF err_sup[i] GT dH_tab[12] AND err_sup[i] LE dH_tab[13] THEN n_km2[i] = 40
    IF err_sup[i] GT dH_tab[13] AND err_sup[i] LE dH_tab[14] THEN n_km2[i] = 43
    IF err_sup[i] GT dH_tab[14] AND err_sup[i] LE dH_tab[15] THEN n_km2[i] = 47
    IF err_sup[i] GT dH_tab[15] AND err_sup[i] LE dH_tab[16] THEN n_km2[i] = 50 
    IF err_sup[i] GT dH_tab[16] AND err_sup[i] LE dH_tab[17] THEN n_km2[i] = 53
    IF err_sup[i] GT dH_tab[17] AND err_sup[i] LE dH_tab[18] THEN n_km2[i] = 57
    IF err_sup[i] GT dH_tab[18] AND err_sup[i] LE dH_tab[19] THEN n_km2[i] = 60
    IF err_sup[i] GT dH_tab[19] AND err_sup[i] LE dH_tab[20] THEN n_km2[i] = 63
    IF err_sup[i] GT dH_tab[20] AND err_sup[i] LE dH_tab[21] THEN n_km2[i] = 67
    IF err_sup[i] GT dH_tab[21] AND err_sup[i] LE dH_tab[22] THEN n_km2[i] = 70
    IF err_sup[i] GT dH_tab[22] AND err_sup[i] LE dH_tab[23] THEN n_km2[i] = 73
    IF err_sup[i] GT dH_tab[23] AND err_sup[i] LE dH_tab[24] THEN n_km2[i] = 77
    IF err_sup[i] GT dH_tab[24] AND err_sup[i] LE dH_tab[25] THEN n_km2[i] = 80
    IF err_sup[i] GT dH_tab[25] AND err_sup[i] LE dH_tab[26] THEN n_km2[i] = 83
    IF err_sup[i] GT dH_tab[26] AND err_sup[i] LE dH_tab[27] THEN n_km2[i] = 87
    IF err_sup[i] GT dH_tab[27] THEN n_km2[i] = 90 
;###############################################################################
    IF err_inf[i] LE dH_tab[1]  THEN n_km1[i] = 0 
    IF err_inf[i] GT dH_tab[1]  AND err_inf[i] LE dh_tab[2]  THEN n_km1[i] = 3 
    IF err_inf[i] GT dH_tab[2]  AND err_inf[i] LE dh_tab[3]  THEN n_km1[i] = 7
    IF err_inf[i] GT dH_tab[3]  AND err_inf[i] LE dH_tab[4]  THEN n_km1[i] = 10 
    IF err_inf[i] GT dH_tab[4]  AND err_inf[i] LE dH_tab[5]  THEN n_km1[i] = 13
    IF err_inf[i] GT dH_tab[5]  AND err_inf[i] LE dH_tab[6]  THEN n_km1[i] = 17
    IF err_inf[i] GT dH_tab[6]  AND err_inf[i] LE dH_tab[7]  THEN n_km1[i] = 20
    IF err_inf[i] GT dH_tab[7]  AND err_inf[i] LE dH_tab[8]  THEN n_km1[i] = 23
    IF err_inf[i] GT dH_tab[8]  AND err_inf[i] LE dH_tab[9]  THEN n_km1[i] = 27
    IF err_inf[i] GT dH_tab[9]  AND err_inf[i] LE dH_tab[10] THEN n_km1[i] = 30
    IF err_inf[i] GT dH_tab[10] AND err_inf[i] LE dH_tab[11] THEN n_km1[i] = 33
    IF err_inf[i] GT dH_tab[11] AND err_inf[i] LE dH_tab[12] THEN n_km1[i] = 37
    IF err_inf[i] GT dH_tab[12] AND err_inf[i] LE dH_tab[13] THEN n_km1[i] = 40
    IF err_inf[i] GT dH_tab[13] AND err_inf[i] LE dH_tab[14] THEN n_km1[i] = 43
    IF err_inf[i] GT dH_tab[14] AND err_inf[i] LE dH_tab[15] THEN n_km1[i] = 47
    IF err_inf[i] GT dH_tab[15] AND err_inf[i] LE dH_tab[16] THEN n_km1[i] = 50 
    IF err_inf[i] GT dH_tab[16] AND err_inf[i] LE dH_tab[17] THEN n_km1[i] = 53
    IF err_inf[i] GT dH_tab[17] AND err_inf[i] LE dH_tab[18] THEN n_km1[i] = 57
    IF err_inf[i] GT dH_tab[18] AND err_inf[i] LE dH_tab[19] THEN n_km1[i] = 60
    IF err_inf[i] GT dH_tab[19] AND err_inf[i] LE dH_tab[20] THEN n_km1[i] = 63
    IF err_inf[i] GT dH_tab[20] AND err_inf[i] LE dH_tab[21] THEN n_km1[i] = 67
    IF err_inf[i] GT dH_tab[21] AND err_inf[i] LE dH_tab[22] THEN n_km1[i] = 70
    IF err_inf[i] GT dH_tab[22] AND err_inf[i] LE dH_tab[23] THEN n_km1[i] = 73
    IF err_inf[i] GT dH_tab[23] AND err_inf[i] LE dH_tab[24] THEN n_km1[i] = 77
    IF err_inf[i] GT dH_tab[24] AND err_inf[i] LE dH_tab[25] THEN n_km1[i] = 80
    IF err_inf[i] GT dH_tab[25] AND err_inf[i] LE dH_tab[26] THEN n_km1[i] = 83 
    IF err_inf[i] GT dH_tab[26] AND err_inf[i] LE dH_tab[27] THEN n_km1[i] = 87
    IF err_inf[i] GT dH_tab[27] THEN n_km1[i] = 90    
ENDFOR
outfile = STRARR(file_number)
        Device_bak2 = !D.Name         
        SET_PLOT, 'X'
       ; print, n_km1    
        ;print, '################################################################'
       ; print, n_km2    
       ; print, '################################################################'        
    ;plot, tot_days, dst
    ;plot, tot_days, new_H, linestyle=0 
    ;oplot, tot_days, H, linestyle=2 
    ;oplot, k_days, n_km2/10.0, psym=2;linestyle=2 
   ; plot, k_days, kp, yrange=[0,9], ystyle=1
   ; oplot, k_days, kp, psym=6  
    ;oplot, k_days, n_km1/10.,psym=1 
    ;oplot, k_days, n_km2/10.,psym=2 
FOR i=0, file_number-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                
    outfile[i] = '../rutidl/Kmex/new_kmex/new_km'+string_date[i]+'.txt'
    OPENW, LUN, outfile[i], /GET_LUN
    PRINTF, LUN, n_km1[i*8:(i+1)*8-1], FORMAT='(8(I02,X))'
    PRINTF, LUN, n_km2[i*8:(i+1)*8-1], FORMAT='(8(I02,X))' 
    CLOSE, LUN
    FREE_LUN, LUN 
ENDFOR
END		
