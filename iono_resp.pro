;
;Name:
;	get_data_date
;purpose:
;	this routine will call read a certain number of files containing 
;   interplanetary data measurements
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
;
;
;parameters:
;
;
;dependencies:
;
;
;input files
;
;
;output files:
;
function tec_data, idate
	On_error, 2
	compile_opt idl2, HIDDEN

	iyear	= idate[0]
	imonth	= idate[1]
	iday 	= idate[2]		

        header = 1      ; Defining number of lines of the header 

        idate = string(iyear, imonth, iday, format = '(I4, "-", I02, "-", I02)')
		file_name = '../rutidl/tec/'+'tec_'+idate+'.txt'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DStruct = {doy : 0, tec : 0., med : 0.}                                   
		r_tec = REPLICATE(DStruct, number_of_lines-header)	      
		READS, data[header:number_of_lines-1], r_tec, $
	format='(F3, X, F5, X, F5)'		
		return, r_tec
end

function dst_data, initial

	On_error, 2
	compile_opt idl2, HIDDEN

	year = string(initial, format = '(I4)')
	;type_data = string(tp, format = '(A)')
		file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_D.dat'
		
        header = 25             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_P.dat'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN begin
    	        file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_Q.dat'
	            file = FILE_SEARCH(file_name, COUNT=opened_files)    	        
    	        IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
    	    ENDIF    	    	    
	    ENDIF

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
        DOY : 0, Dst: 0}

		r_dst = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], r_dst, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,8X,I3,X,I6)'
		
		return, r_dst
end

function DH_teo, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4, I02, I02)')
       ; sts  = string(stats, format = '(A5)')
		
		file_name = '../rutidl/dH_teo/'+'teo_'+date+'.dst.early'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        DStruct = {hora : 0, D_stdesv : 0., D : 0., H_stdesv : 0., H : 0., $
        Z_stdesv : 0., Z : 0., N_stdesv : 0., N : 0., F_stdesv : 0., F : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(I2, F10, F8, F10, F10, F10, F10, F10, F10, F10, F10)'		
		return, teo_mag		
end



function baseline_sq, date
	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr	= date[0]
	mh	= date[1]

        date = string(yr, mh, format = '(I4, "-", I02)')
        header=0

		file_name = '../rutidl/output/Bsq_baselines/'+'Bsq_'+date+'.txt'
	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
		
        DStruct = {year : 0, month : 0, day : 0, hour : 0, doy : 0, Bsq : 0.}                                   

		B_sq = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], B_sq, $
	 format='(I4,X, I02,X, I02, 2X, I02, 2X, I03, 2X, F08.4)' 		
		return, B_sq
end


function Date2DOY, idate
;	Check data type of input set ascII flag and convert to yy,mm,dd:
	info = SIZE(idate)
	IF (info(0) eq 0) THEN BEGIN
	  scalar = 1				;scalar flag set
	ENDIF ELSE BEGIN
	  scalar = 0				;vector input
	ENDELSE

	IF (info(info(0) + 1) eq 7) THEN BEGIN
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
	  IF ((yy MOD 4) eq 0) THEN BEGIN	;leap year
	    imonth(2) = 29			;set feb
	  ENDIF
	  DOY = FIX( TOTAL(imonth(0:mm-1)) ) + dd
	ENDIF ELSE BEGIN
	  DOY = dd				;set correct len on vector
	  leapYrs = WHERE( (yy MOD 4) eq 0)	;index of leap years
	  nonLeap = WHERE( (yy MOD 4) ne 0)	;index of non-leap years
	  IF (nonLeap(0) ne -1) THEN BEGIN
	    FOR i=0, N_elements(nonLeap)-1 DO BEGIN
	      DOY(nonLeap(i)) = FIX( TOTAL(imonth(0:mm(nonLeap(i))-1)) ) + $
				dd(nonLeap(i))
	    ENDFOR
	  ENDIF
	  IF (leapYrs(0) ne -1) THEN BEGIN
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
	;print, uint(DOY)  
	return, DOY
	END
	
	
pro iono_resp, r_dst, B_sq, DOY, date_i, date_f 

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
;###############################################################################
;###############################################################################
idate0 = string(yr_i, mh_i, format='(I4,I02)')
TGM_n = idate0
case TGM_n of
    '200311' : TGM_n = 1
    '200411' : TGM_n = 2
    '200505' : TGM_n = 3
    '201503' : TGM_n = 4
    '201705' : TGM_n = 5
    '201709' : TGM_n = 6
    else: print, 'fuera de rango'
endcase  
;###############################################################################
;###############################################################################
    d_dst = dst_data(yr_i)     
    i_dst   = d_dst.Dst
    hour    = d_dst.hour
    dst_doy = d_dst.DOY          
;###############################################################################   
    d_dst = dst_data(yr_i)
    t = n_elements(d_dst.year)    
    i_dst = d_dst.Dst
   
    year = d_dst.year
    tiempo = TIMEGEN(t, START=julday(d_dst.month[0], d_dst.day[0],  $
                     d_dst.year[0], d_dst.hour[0]), UNITS='Hours')  
                                        
        iyear = strmid(string(yr_i, format='(I4)'),2,2)
        fyear = strmid(string(yr_f, format='(I4)'),2,2)
                
        idoy      = Date2DOY(string(iyear, mh_i, dy_i,format = '(I02,I02,I02)'))
        fdoy      = Date2DOY(string(fyear, mh_f, dy_f,format = '(I02,I02,I02)'))         
            
    time_w  = tiempo[idoy:fdoy]
    tw      = n_elements(time_w)
    tot_days= findgen(tw*24)/24.0
    
    dst     = i_dst[(idoy*24)-24:fdoy*24-1]    
    Date    = string(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
;###############################################################################
       file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
               
        data_file_name_dh  = strarr(file_number)        
        string_date_dh        = strarr(file_number)

        data_file_name_tec  = strarr(file_number)        
        string_date_tec        = strarr(file_number)        
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date_dh[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date_tec[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4, "-", I02, "-", I02)')
        
                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date_dh[i]+'.dst.early'
                data_file_name_tec[i] = '../rutidl/tec/'+'tec_'+string_date_tec[i]+'.txt'
		        file_dh = FILE_SEARCH(data_file_name_dh[i], COUNT=opened_files)
		        file_tec = FILE_SEARCH(data_file_name_tec[i], COUNT=opened_files)
		        
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date_dh[i]+'.dst.early'    
	            ENDIF
	            
	            IF opened_files NE N_ELEMENTS(file_tec) THEN begin
	                data_file_name_tec[i] = '../rutidl/tec/'+'tec_'+string_date_tec[i]+'.txt'
	            ENDIF 	                            
        ENDFOR

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        exist_data_file_tec   = FILE_TEST(data_file_name_tec)
        capable_to_plot_tec   = N_ELEMENTS(where(exist_data_file_tec EQ 1))
        
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_tec NE N_ELEMENTS(data_file_name_tec) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF
                        
        H    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_dh)-1 DO BEGIN
                IF exist_data_file_dh[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date_dh[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_dh = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_dh.H[*]
                                                                                                                       
                ENDIF ELSE BEGIN
                        H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
        
        tec  = fltarr(file_number*12)
        med  = fltarr(file_number*12)
        FOR i = 0, N_ELEMENTS(exist_data_file_tec)-1 DO BEGIN
                IF exist_data_file_tec[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date_tec[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,X,I02,X,I02)'                 
                        d_tec = tec_data([tmp_year, tmp_month, tmp_day])
                        
                        tec[i*12:(i+1)*12-1] = d_tec.tec[*]
                        med[i*12:(i+1)*12-1] = d_tec.med[*] 
                ENDIF ELSE BEGIN
                        tec[i*12:(i+1)*12-1] = 999.0
                        med[i*12:(i+1)*12-1] = 999.0
                ENDELSE                
        ENDFOR

        for i=0, n_elements(tec)-1 do begin
            if tec[i] eq 999.0 then begin
                tec[where(tec[*] eq 999.0)] = !Values.F_NAN          
            endif
        endfor    

        for i=0, n_elements(med)-1 do begin
            if med[i] eq 999.0 then begin
                med[where(med[*] eq 999.0)] = !Values.F_NAN          
            endif
        endfor
        
        i_nan1 = where(H eq 999999.0, ncount)
        i_nan2 = where(H gt 100.0, n2count)
        
        prcent_nan = FLOAT(ncount+n2count)*100.0
        print,'porcentaje de valores NaN:', prcent_nan/n_elements(H),'%'
        
        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(H)-1 do begin
            if H[i] ge 100.0 then begin
                H[where(H[*] ge 100.0)] = !Values.F_NAN          
            endif
        endfor        
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo
       
    H_tmp   = H
    H_exist = where(finite(H_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
    ; interpolate at the locations of the bad data using the good data
    
    if nbaddata gt 0 then H_tmp[baddata] = interpol(H_tmp[H_exist], H_exist, baddata)
    H = H_tmp
;###############################################################################      
    tec_days= findgen(tw*12)/12.0                         
    tec_diff = tec-med    
    new_tecdays = findgen(tw*1440)/1440.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_tecdiff = FLTARR(N_ELEMENTS(new_tecdays))     	    
    tmp_tecdif  = INTERPOL(tec_diff, N_ELEMENTS(new_tecdays))
    new_tecdiff = tmp_tecdif
;###############################################################################
    new_dstdays = findgen(tw*1440)/1440.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_dst = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_dst  = INTERPOL(dst, N_ELEMENTS(new_dstdays))
    new_dst = tmp_dst 
    
    new_H = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_H  = INTERPOL(H, N_ELEMENTS(new_dstdays))
    new_H = tmp_H 
;###############################################################################    
;############################################################################### 
; define diurnal baseline
;###############################################################################  
    sqline = baseline_sq([yr_i, mh_i])
    
    sq_doy  = sqline.doy
    Bsq_ln   = sqline.Bsq


    tmp_doy = dst_doy[(idoy*24)-24:fdoy*24-1]
    Bsq     = fltarr(n_elements(Bsq_ln))
    
    tmp_doy = tmp_doy[uniq(tmp_doy, sort(tmp_doy))]
   ; print, tmp_doy
    
   ; tmp_doy2= intarr(n)
    for i=0, n_elements(tmp_doy)-1 do begin
        ;print, tmp_doy[i]
            for j=0, n_elements(sq_doy)-1 do begin
           ; print, ip_doy[j]
                if tmp_doy[i] eq sq_doy[j] then begin
               ; ip_doy[j]   = tmp_doy[i]   
               Bsq[j]    = Bsq_ln[j]                            
                endif 
            endfor
    endfor
b_sq = where(Bsq eq 0, zcount, complement=val, ncomplement=valcount)

Bsq      = Bsq[val]
;###############################################################################
; define frequencies
;###############################################################################  
   mlat         = 28.06*!pi
   ld           = cos(mlat/180)
   p_a          = dst*ld
   baseline     = Bsq + p_a             
        diono   = H-baseline
    n           = n_elements(diono) 
time = 3600.0

fn      = float(1.0/(2.0*time)) ; frecuencia de Nyquist
y       = FFT(diono)

pws     = abs(y[0:n/2])^2
pws_s   = smooth(pws, 1)

f_k     = (1+findgen(n))/(n*time)
print, 'Nyquist freq: ', fn, 'Hz'
;###############################################################################
    i_diff = diono
    new_idiff = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_idiff  = INTERPOL(i_diff, N_ELEMENTS(new_dstdays))
    new_idiff = tmp_idiff  
;###############################################################################
; define pass band frequencies
;###############################################################################  
passband_l = idate0
case passband_l of
    '200311' : passband_l = 1.22e-5
    '200411' : passband_l = 1.32e-5
    '200505' : passband_l = 1.2e-5
    '201503' : passband_l = 1.28e-5
    '201705' : passband_l = 1.25e-5
    '201709' : passband_l = 1.22e-5
    else: print, 'fuera de rango'
endcase  

passband_u = idate0
case passband_u of
    '200311' : passband_u = 1.5e-5
    '200411' : passband_u = 1.55e-5
    '200505' : passband_u = 1.42e-5
    '201503' : passband_u = 1.6e-5
    '201705' : passband_u = 1.55e-5
    '201709' : passband_u = 1.58e-5
    else: print, 'fuera de rango'
endcase  
;###############################################################################
; define high band frequencies
;###############################################################################
highpass_l = idate0
case highpass_l of
    '200311' : highpass_l = 7e-5
    '200411' : highpass_l = 7.5e-5
    '200505' : highpass_l = 7.3e-5
    '201503' : highpass_l = 7.5e-5
    '201705' : highpass_l = 7.5e-5
    '201709' : highpass_l = 7.e-5
    else: print, 'fuera de rango'
endcase
;###############################################################################
; define filtering
;############################################################################### 
coeff_ddyn    = digital_filter(passband_l/fn, passband_u/fn, 50, 18)
coeff_dp2   = digital_filter(highpass_l/fn, 1.0, 50, 4)
;###############################################################################
; define disturbing effects
;############################################################################### 
ddyn        = convol(diono, coeff_ddyn, /edge_wrap)
dp2         = convol(diono, coeff_dp2, /edge_wrap)  
;###############################################################################      
    new_ddyn = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_ddyn  = INTERPOL(ddyn, N_ELEMENTS(new_dstdays))
    new_ddyn = tmp_ddyn           
;###############################################################################
;############################################################################### 
    new_dp2 = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_dp2  = INTERPOL(dp2, N_ELEMENTS(new_dstdays))
    new_dp2 = tmp_dp2         
;###############################################################################
;###############################################################################
; define device and color parameters 
;###############################################################################      
        Device_bak2 = !D.Name 
        
        SET_PLOT, 'Z'      
        
        Xsize=fix(1600)
        Ysize=1000
        DEVICE, SET_RESOLUTION = [Xsize,Ysize],Set_Pixel_Depth=24, DECOMPOSED=0  
        DEVICE, z_buffer=1
        DEVICE, set_character_size = [10, 12] 
        
        ;SET_PLOT, 'X'            
        ;DEVICE, RETAIN=2
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 190
        verde     = 180
        negro     = 0
        azul      = 70
        blanco    = 255
        gris      = 220
        morado    = 16
        
  ;  TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 0, /SILENT
     X_label   = STRARR(tw+1)+' '
        months    = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic']
        old_month = mh_i
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR 
        
case old_month of
    1: old_month = 'Enero'
    2: old_month ='Febrero'
    3: old_month ='Marzo'
    4: old_month ='Abril'
    5: old_month ='Mayo'
    6: old_month ='Junio'
    7: old_month ='Julio'
    8: old_month ='Agosto'
    9: old_month ='Septiembre'
    10:old_month ='Octubre'
    11:old_month ='Noviembre'
    12:old_month ='Diciembre'
    else: print, 'fuera de rango'
endcase             
;###############################################################################
;###############################################################################
    spam_i = idate0
case spam_i of
    '200311' : spam_i = 1000
    '200411' : spam_i = 50
    '200505' : spam_i = 100
    '201503' : spam_i = 0
    '201705' : spam_i = 0
    '201709' : spam_i = 100
    else: print, 'fuera de rango'
endcase 

    spam_f = idate0
case spam_f of
    '200311' : spam_f = 3100
    '200411' : spam_f = 4400
    '200505' : spam_f = 800
    '201503' : spam_f = 3300
    '201705' : spam_f = 1440
    '201709' : spam_f = 1880
    else: print, 'fuera de rango'
endcase    
;############################################################################### 
;###############################################################################
       days = intarr(tw+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       days = days*24/24. 
       day_time = findgen(24)   
;############################################################################### 
    time_title = ' Tiempo universal (dias).'
    window_title = 'TGM'+ string(TGM_n, format='(I01)')+', '+ $
                    string(old_month, yr_i, format='(A, X, I4)')

    plot_title1 = 'Espectro de potencias de la perturbacion ionosferica (Diono)'
    plot_title2 = 'Perturbaciones Ionosfericas DP2 y DDyn'
    
    periodo = 'Periodo (Hr)'        
;###############################################################################
;###############################################################################               
    plot, f_k, pws_s, /xlog, /ylog, POSITION=[0.07,0.1,0.95,0.9],$
    BACKGROUND = blanco, color=negro, $
    CHARSIZE = chr_size1, xstyle=5, ystyle=5, subtitle='', thick=4, /NODATA    

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.95   
   XYOuts, X, y, window_title, /Normal, $
   color=negro, Alignment=0.5, Charsize=1.65               
;###############################################################################
;###############################################################################           
    ysup = max(pws)+10
    yinf = min(pws)-0.0001
    
    plot, f_k, pws_s, /xlog, /ylog, xrange = [min(f_k), fn], POSITION=[0.07,0.1,0.45,0.9],$
    yrange=[yinf, ysup], BACKGROUND = blanco, color=negro, $
    CHARSIZE = chr_size1, xstyle=5, ystyle=5, subtitle='', thick=4, /NODATA,$
    /NOERASE

    ;LOADCT, 0, /SILENT

    POLYFILL, [passband_l, passband_u ,passband_u, passband_l], $
              [!Y.CRANGE[0], !Y.CRANGE[0], ysup, ysup], color=amarillo

    POLYFILL, [highpass_l, fn ,fn, highpass_l], $
              [!Y.CRANGE[0], !Y.CRANGE[0], ysup, ysup], color=amarillo
    oplot, f_k, pws_s, color=negro, thick=5

    ;LOADCT, 0, /SILENT
        AXIS, XAXIS = 0, XRANGE=[min(f_k), fn], $
                         /xlog,$
                         xstyle=1,$
                         xTITLE = 'frecuencias [Hz]',$
                         COLOR=negro, $
                         CHARSIZE = 1.0, $
                         TICKLEN=0.04
      
                         
        AXIS, XAXIS = 1, XRANGE=[min(f_k), fn], $
                         /xlog,$
                         XTICKN=['27:46','02:46'],$
                         xstyle=1,$
                         CHARSIZE = 1.0,$
                         COLOR=negro, $
                         TICKLEN=0.04                     

        AXIS, YAXIS = 0, yrange=[yinf, ysup], $
                         YTITLE = 'Componente espectral [nT]', $
                         ystyle=1,$                          
                         COLOR=negro, $
                         /ylog,$
                         CHARSIZE = 1.0;, $
                        
        AXIS, YAXIS = 1, yrange=[yinf, ysup], $
                         COLOR=negro, $
                         /ylog,$
                         ystyle=1, $
                         CHARSIZE = 1.0;, $


   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.035, plot_title1, /Normal, $
   color=negro, Alignment=0.5, Charsize=1.45  
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.924, periodo, /Normal, $
   color=negro, Alignment=0.5, Charsize=1.0   
;###############################################################################
;###############################################################################
    med_idx = MEDIAN(new_idiff)
    std_idx = stddev(new_idiff, /NAN)
    
    i_out = WHERE(new_idiff GE med_idx+std_idx OR new_idiff LE med_idx-std_idx)
    i_in  = WHERE(new_idiff LE med_idx+std_idx AND new_idiff GE med_idx-std_idx)
    id_diff_out = new_idiff
    id_diff_out[i_in]=!Values.F_NAN
    
    id_diff_in  = new_idiff
    id_diff_in[i_out]=!Values.F_NAN

    sup0 = med_idx+std_idx
    inf0 = med_idx-std_idx
    
    lim_sup = fltarr(n_elements(new_idiff))
    lim_sup[*] = sup0

    lim_inf = fltarr(n_elements(new_idiff))
    lim_inf[*] = inf0   
                    
      
    up_diono = max(diono)
    down_diono = min(diono)

    up_tecdiff = max(tec_diff) 
    down_tecdiff = min(tec_diff)     
;###############################################################################         
     up = max(H)
     down=min(H)
     plot, tot_days, H, XTICKS=file_number, xminor=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = 0.9, CHARTHICK=chr_thick1, $
     POSITION=[0.55,0.73,0.95,0.9], XSTYLE = 5, XRANGE=[0, tw], ySTYLE = 6,$
     XTICKNAME=REPLICATE(' ', tw+1), yrange=[down,up], /NOERASE, THICK=3, /NODATA
       
    POLYFILL, [new_dstdays[i_out[0]+spam_i], new_dstdays[i_out[0]+spam_f] ,$
              new_dstdays[i_out[0]+spam_f], new_dstdays[i_out[0]+spam_i]], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], color=amarillo     
     
     OPLOT, new_dstdays, new_dst, COLOR=azul, THICK=3
     OPLOT, new_dstdays, new_H, COLOR=negro, THICK=3  
     
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         XTITLE=time_title, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=tw, $
                         XTICKV=days,$
                         XTICKN=fix(days),$                         
                         CHARSIZE = 0.8,$
                         XMINOR=8, $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down,up], $
                         YTITLE = 'DH [nT]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $
                         YTITLE = 'Dst [nT]', $         
                         COLOR=azul, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $

;###############################################################################      
     up_diono=max(diono)
     down_diono=min(diono)          
     plot, tot_days, diono, XTICKS=file_number, xminor=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = 0.6, CHARTHICK=chr_thick1, $
     POSITION=[0.55,0.49,0.95,0.66], XSTYLE = 5, XRANGE=[0, tw], ySTYLE = 6,$
     XTICKNAME=REPLICATE(' ', tw+1), yrange=[down_diono,up_diono], /NOERASE,$
     /NODATA

       ; oplot, new_dstdays[i_out[0]+spam_i:i_out[0]+spam_f], $
       ; id_diff_in[i_out[0]+spam_i:i_out[0]+spam_f], $
        ;color=negro, linestyle=0, thick=4
        
        ;oplot, new_dstdays[i_out[0]+spam_i:i_out[0]+spam_f], $
        ;id_diff_out[i_out[0]+spam_i:i_out[0]+spam_f], $
        ;color=negro, linestyle=0, thick=4
        
        oplot, new_dstdays, id_diff_in, color=negro, linestyle=3
        oplot, new_dstdays, id_diff_out, color=negro, linestyle=0, thick=4       
;###############################################################################
;###############################################################################
    med_tec = MEDIAN(new_tecdiff)
    std_tec = stddev(new_tecdiff)
    
    index_out = WHERE(new_tecdiff GE med_tec+std_tec OR new_tecdiff LE med_tec-std_tec)
    index_in  = WHERE(new_tecdiff LE med_tec+std_tec AND new_tecdiff GE med_tec-std_tec)
    tec_diff_out = new_tecdiff
    tec_diff_out[index_in]=!Values.F_NAN
    
    tec_diff_in  = new_tecdiff
    tec_diff_in[index_out]=!Values.F_NAN

    sup = med_tec+std_tec
    inf = med_tec-std_tec
    
    l_sup = fltarr(n_elements(tec_diff))
    l_sup[*] = sup

    l_inf = fltarr(n_elements(tec_diff))
    l_inf[*] = inf   
       
    plot, tec_days, tec_diff, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=rojo,$
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.55,0.49,0.95,0.66], $
     XSTYLE = 5, XRANGE=[0, tw], XTICKNAME=REPLICATE(' ', tw+1), ySTYLE = 6,$
     /NOERASE, /NODATA, YRANGE=[down_tecdiff, up_tecdiff]
    
        tecdiff_inicio = index_out[0]        
        oplot, new_tecdays, tec_diff_in, color=rojo, linestyle=0
        ;oplot, new_tecdays[tecdiff_inicio:tecdiff_inicio+2880], tec_diff_out[tecdiff_inicio:tecdiff_inicio+2880], color=rojo, linestyle=0, thick=4
        oplot, new_tecdays, tec_diff_out, color=rojo, linestyle=0, thick=4                
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=tw, $
                         XTICKV=days,$
                         XTICKN=fix(days),$                         
                         XMINOR=8, $ 
                         CHARSIZE = 0.8, $                       
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_diono,up_diono], $
                         YTITLE = 'Diono [nT]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down_tecdiff, up_tecdiff], $
                         YTITLE = 'TEC-<TEC> [TECu]', $          
                         COLOR=rojo, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $                                               
;###############################################################################
;###############################################################################                
    if max(ddyn) gt max(dp2) then up = max(ddyn) else up = max(dp2)
    if min(ddyn) lt min(dp2) then down = min(ddyn) else down = min(dp2)
;###############################################################################
;###############################################################################
    med_ddyn = MEDIAN(new_ddyn)
    std_ddyn = stddev(new_ddyn, /NAN)
    
    ddyn_out = WHERE(new_ddyn GE med_ddyn+std_ddyn OR new_ddyn LE med_ddyn-std_ddyn)
    ddyn_in  = WHERE(new_ddyn LE med_ddyn+std_ddyn AND new_ddyn GE med_ddyn-std_ddyn)
    
    ddyn_diff_out = new_ddyn
    ddyn_diff_out[ddyn_in]=!Values.F_NAN
    
    ddyn_diff_in  = new_ddyn
    ddyn_diff_in[ddyn_out]=!Values.F_NAN

    supddyn = med_ddyn+std_ddyn
    infddyn = med_ddyn-std_ddyn
    
    lim_supddyn = fltarr(n_elements(new_ddyn))
    lim_supddyn[*] = supddyn

    lim_infddyn = fltarr(n_elements(new_ddyn))
    lim_infddyn[*] = infddyn   
                    
     upddyn     = max(ddyn)
     downddyn   = min(ddyn)
     
     updp2     = max(dp2)
     downdp2   = min(dp2)     

    if upddyn gt updp2 then up = upddyn else up=updp2 
    if downddyn gt downdp2 then down = down else down=downdp2 
                               
     plot, tot_days, ddyn, XTICKS=file_number, xminor=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
     POSITION=[0.55,0.1,0.95,0.42], XSTYLE = 5, XRANGE=[0, tw], ySTYLE = 6,$
     XTICKNAME=REPLICATE(' ', tw+1), yrange=[down,up], /NOERASE, /NODATA
    
     ddyn_inicio = ddyn_out[150]
     
       ; oplot, new_dstdays[ddyn_inicio:ddyn_inicio+spam_f], $
       ; ddyn_diff_out[ddyn_inicio:ddyn_inicio+spam_f], color=negro, $
       ; linestyle=0, thick=5   
                    
       ; oplot, new_dstdays[ddyn_inicio:ddyn_inicio+spam_f], $
       ; ddyn_diff_in[ddyn_inicio:ddyn_inicio+spam_f], color=negro, linestyle=0, $
       ; thick=5          

        oplot, new_dstdays, ddyn_diff_in, color=negro, linestyle=3
        oplot, new_dstdays, ddyn_diff_out, color=negro, linestyle=0, thick=5
;###############################################################################
;###############################################################################     
;###############################################################################
    med_dp2 = MEDIAN(new_dp2)
    std_dp2 = stddev(new_dp2, /NAN)
    
    dp2_out = WHERE(new_dp2 GE med_dp2+std_dp2 OR new_dp2 LE med_dp2-std_dp2)
    dp2_in  = WHERE(new_dp2 LE med_dp2+std_dp2 AND new_dp2 GE med_dp2-std_dp2)
    
    dp2_diff_out = new_dp2
    dp2_diff_out[dp2_in]=!Values.F_NAN
    
    dp2_diff_in  = new_dp2
    dp2_diff_in[dp2_out]=!Values.F_NAN

    supdp2 = med_dp2+std_dp2
    infdp2 = med_dp2-std_dp2
    
    lim_supdp2 = fltarr(n_elements(new_dp2))
    lim_supdp2[*] = supdp2

    lim_infdp2 = fltarr(n_elements(new_dp2))
    lim_infdp2[*] = infdp2       
;###############################################################################
    dp2_inicio = dp2_out[300]   
;###############################################################################     
     ;oplot, tot_days, dp2, color=rojo
      ;  oplot, new_dstdays[dp2_inicio+spam_i:dp2_inicio+spam_f], $
       ; new_dp2[dp2_inicio+spam_i:dp2_inicio+spam_f], $
;        color=rojo, linestyle=0, thick=4      
        
 ;       oplot, new_dstdays[dp2_inicio+spam_i:dp2_inicio+spam_f], $
  ;      new_dp2[dp2_inicio+spam_i:dp2_inicio+spam_f], color=rojo, $
   ;     linestyle=0, thick=4    
        
        oplot, new_dstdays, dp2_diff_in, color=rojo, linestyle=1         
        oplot, new_dstdays, dp2_diff_out, color=rojo, linestyle=0, thick=4       

        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKV=days,$
                         XTICKN=fix(days),$  
                         CHARSIZE = 0.8, $                                                
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, yrange=[down,up], $ 
                         ystyle=2, $  
                         YTITLE = 'Ddyn [nT]', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, yrange=[down,up], $ 
                         ystyle=2, $ 
                         YTITLE = 'DP2 [nT]', $                           
                         COLOR=rojo, $
                         CHARSIZE = 0.9;, $      
;###############################################################################    
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.921, 'Tiempo Local (dias)', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9  

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.679, 'Tiempo Local (dias)', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9  
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.441, 'Tiempo Local (dias)', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9  
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.035, plot_title2, /Normal, $
   color=negro, Alignment=0.5, Charsize=1.45  
    LOADCT, 39, /SILENT     
;###############################################################################
; saving png
;###############################################################################     
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!    
    TVLCT, R_bak, G_bak, B_bak, /get  
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak2  
    path = '../rutidl/output/eventos_tgm/'
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = R_bak[image]
                true_image[1,*,*] = G_bak[image]
                true_image[2,*,*] = B_bak[image]
                write_jpeg, path+'iono_resp_V5_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'iono_resp_V5_'+Date+'.png', Image, R_bak, G_bak, B_bak
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'iono_resp_V5_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	
end	
