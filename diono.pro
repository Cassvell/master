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
		file_name = '../rutidl/tec/tec_newformat/'+'tec_'+idate+'.txt'
		
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


function ip_data, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	

        header = 60      ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

        date = string(year, month, day, format = '(I4, "-", I02, "-", I02)')
		
		file_name = '../rutidl/ip/'+date+'.dat'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DStruct = {YEAR : 0, DOY : 0, hour : 0, bartels : 0, id_ifm : 0, $
                      id_sw : 0, n_points_ifm: 0, n_points_plasna : 0, B_esc : 0., $
                      B_vec : 0., B_lat : 0., B_long : 0., Bx : 0., By : 0., Bz : 0., $
                      By0 : 0., Bz0 : 0., B_rms_esc : 0., B_rms_vec : 0., Bx_rms : 0., $
                      By_rms : 0., Bz_rms : 0., Temp : 0., density_p : 0., v_p : 0., $
                      flow_long : 0., flow_lat : 0., alfa_prot : 0., sigma_T : 0., $
                      sigma_n : 0., sigma_v : 0., sigma_phi : 0., sigma_theta : 0., $
                      sigma_ratio : 0., flow_P : 0., E : 0., beta_p : 0., alfven_M : 0., $
                      mag_M : 0., q_invariant : 0., kp : 0, R : 0, dst : 0, ap : 0, $
                      f10_idx : 0., AE : 0, AL : 0, AU : 0, pc_idx : 0., lyman_alfa : 0.,$
                      p_flux_1MeV : 0., p_flux_2MeV : 0., p_flux_4MeV : 0., p_flux_10MeV : 0.,$
                      p_flux_30MeV : 0., p_flux_60MeV : 0., flux_FLAG : 0}
                                      

		r_ip = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], r_ip, $
FORMAT='(I4,I4,I3,I5,I3,I3,I4,I4,F6,F6,F6,F6,F6,F6,F6,F6,F6,F6,F6,F6,F6,F6,'+$
'F9,F6,F6,F6,F6,F6,F9,F6,F6,F6,F6,F6,F6,F7,F7,F6,F5,F7,I3,I4,I6,I4,F6,I5,'+$
'I6,I6,F6,F9,F10,F9,F9,F9,F9,F9,I3)'		
		return, r_ip
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


pro diono, r_dst, r_ip, B_sq, DOY, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]

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
                data_file_name_tec[i] = '../rutidl/tec/tec_newformat/'+'tec_'+string_date_tec[i]+'.txt'
		        file_dh = FILE_SEARCH(data_file_name_dh[i], COUNT=opened_files)
		        file_tec = FILE_SEARCH(data_file_name_tec[i], COUNT=opened_files)
		        
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date_dh[i]+'.dst.early'    
	            ENDIF
	            
	            IF opened_files NE N_ELEMENTS(file_tec) THEN begin
	                data_file_name_tec[i] = '../rutidl/tec/tec_newformat/'+'tec_'+string_date_tec[i]+'.txt'
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
    dif_tec = tec-med    
;###############################################################################
; define diurnal baseline
;###############################################################################  
    sqline = baseline_sq([yr_i, mh_i])
    
    sq_doy  = sqline.doy
    Bsq_ln   = sqline.Bsq


    tmp_doy = dst_doy[(idoy*24)-23:fdoy*24-1]
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
; define time window
;###############################################################################  

   mlat         = 28.10*!pi
   ld           = cos(mlat/180)
   p_a          = dst*ld
   baseline     = Bsq + p_a             
        diono  = H-baseline         
;###############################################################################
; define ip parameters
;###############################################################################  
    ip = ip_data([yr_i, mh_i, dy_i])
    
    ip_year = ip.year
    ip_doy  = ip.DOY
    ip_hour = ip.hour
    ip_Ey   = ip.E
    ip_flow = ip.flow_p
    ip_AU   = ip.AU
    ip_AL   = ip.AL
    ip_AE   = ip.AE

    tmp_doy = dst_doy[(idoy*24)-24:idoy*24-1]
    Ey      = fltarr(n_elements(ip_doy))
    p_dyn   = fltarr(n_elements(ip_doy))
    AE      = fltarr(n_elements(ip_doy))
    AL      = fltarr(n_elements(ip_doy))
    AU      = fltarr(n_elements(ip_doy))
        
    tmp_doy = tmp_doy[uniq(tmp_doy, sort(tmp_doy))]
   ; print, tmp_doy
    
    ;tmp_doy2= intarr(n)
    for i=0, n_elements(tmp_doy)-1 do begin
        ;print, tmp_doy[i]
            for j=0, n_elements(ip_doy)-1 do begin
           ; print, ip_doy[j]
                if tmp_doy[i] eq ip_doy[j] then begin
               ; ip_doy[j]   = tmp_doy[i]   
               p_dyn[j]    = ip_flow[j]
               Ey[j]       = ip_Ey[j]
               AE[j]       = ip_AE[j]
               AL[j]       = ip_AL[j]
               AU[j]       = ip_AU[j]                              
                endif 
            endfor
    endfor
    
e_z = where(Ey eq 0, zcount, complement=val, ncomplement=valcount)
p_z = where(p_dyn eq 0, zcount2, complement=val_p, ncomplement=valcount_p)
aez = where(AE eq 0, zcount3, complement=val_AE, ncomplement=valcount)
alz = where(AL eq 0, zcount4, complement=val_AL, ncomplement=valcount)
auz = where(AU eq 0, zcount5, complement=val_AU, ncomplement=valcount)

Ey      = Ey[val]
p_dyn   = p_dyn[val]
AE      = AE[val_AE]
AL      = AL[val_AL]
AU      = AU[val_AU]

        for i=0, n_elements(Ey)-1 do begin
            if Ey[i] eq 999.990 then begin
                Ey[where(Ey[*] eq 999.990)] = !Values.F_NAN          
            endif
        endfor

        for i=0, n_elements(p_dyn)-1 do begin
            if p_dyn[i] eq 99.99 then begin
                p_dyn[where(p_dyn[*] eq 99.99)] = !Values.F_NAN          
            endif
        endfor  
 
;###############################################################################
; define device and color parameters 
;###############################################################################      
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        
        Xsize=fix(800*tmp_spam)
        Ysize=800
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
             
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 220
        verde     = 180
        negro     = 0
        azul      = 80
        blanco    = 255
        gris      = 130
        morado    = 16
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT

    path = '../rutidl/output/eventos_tgm/'
;###############################################################################
; Time label
;###############################################################################    
     X_label   = STRARR(tw+1)+' '
    ; print, n_elements(x_label)
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = mh_i
       ; print, old_month
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
      
      

;###############################################################################
; Plot data
;###############################################################################
    plot_title = 'Perturbación Ionosférica asociada a la TGM'   

    if max(dst) gt max(H) then up = max(dst) else up = max(H)
    if min(dst) lt min(H) then down = min(dst) else down = min(H)
    
window_title = 'from '+string(yr_i, mh_i, dy_i, $
                FORMAT='(I4, "/", I02, "/", I02)')+' to '$
                +string(yr_f, mh_f, dy_f, $
                FORMAT='(I4, "/", I02, "/", I02)')

time_title = ' Time  [UT]'                
                
MAG_source = 'Source: International Service of Geomagnetic Indices'  
        
    plot, tot_days, H, XTICKS=tw, xminor = 8, POSITION=[0.07,0.83,0.95,0.95],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down, up],$
    title = window_title, BACKGROUND = blanco, COLOR=negro, XRANGE=[0, tw],$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', tw+1)

    oplot, tot_days, dst, COLOR=azul, linestyle=0
    ;oplot, tot_days, sym_H, COLOR=azul, linestyle=0    
    
    if up-down gt 300 then begin
        for i = -600., 100., 100. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -600., 100., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse
        
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down,up], $
                         YTITLE = 'Dst and DH [nT]', $
                         ystyle=2, $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]

up_ae       = max(AE)
down_ae     = min(AL)
    plot, tot_days, AE, XTICKS=tw, xminor = 8, POSITION=[0.07,0.70,0.95,0.82],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down_ae, up_ae],$
    BACKGROUND = blanco, COLOR=negro, XRANGE=[0, tw],$
	ytitle = '',  XTICKNAME=REPLICATE(' ', tw+1), /noerase

    oplot, tot_days, AL, COLOR=azul, linestyle=0
    oplot, tot_days, AU, COLOR=rojo, linestyle=0    

        for i = -4000., 4000., 1000. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro

        
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_ae,up_ae], $
                         YTITLE = 'AE, AU and AL [nT]', $
                         ystyle=2, $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down_ae,up_ae], $
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]
  
 up_idx = max(diono)
 down_idx = min(diono)
      
    plot, tot_days, diono, XTICKS=tw, xminor = 8, POSITION=[0.07,0.57,0.95,0.69],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, BACKGROUND = blanco, $
    COLOR=negro, XRANGE=[0, tw], ytitle = 'Indice DST [nT]',$
	XTICKNAME=REPLICATE(' ', tw+1), /noerase, YRANGE=[down_idx, up_idx]
    
    ;oplot, tot_days, diono, COLOR=azul, linestyle=0
    
    if up_idx-down_idx ge 100 then begin
        for i = -600., 100., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -600., 100., 25. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse
    
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0,  YRANGE=[down_idx, up_idx] , $
                         ystyle=2, $ 
                         YTITLE = 'diono [nT]', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down_idx, up_idx] , $
                         ystyle=2, $ 
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]   	

    if max(med) gt max(tec) then up0 = max(med) else up0 = max(tec)
    if min(med) lt min(tec) then down0 = min(med) else down0 = min(tec)
    
    plot, tec_days, tec, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.07,0.44,0.95,0.56], $
     XSTYLE = 5, XRANGE=[0, tw], XTICKNAME=REPLICATE(' ', tw+1), ySTYLE = 6, $
     /noerase, YRANGE=[down0,up0]
     
    oplot, tec_days, med, COLOR=rojo
    
    if up0-down0 gt 120 then begin
        for i = -20., 200., 40. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -20., 200., 20. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 
          
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0], $
                         YTITLE = 'TECU', $
                         ystyle=2, $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        
                        
        AXIS, YAXIS = 1, YRANGE=[down0,up0], $
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
                         


    up_diff = max(dif_tec) 
    down_diff = min(dif_tec)
    
    plot, tec_days, dif_tec, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.07,0.31,0.95,0.43], $
     XSTYLE = 5, XRANGE=[0, tw], XTICKNAME=REPLICATE(' ', tw+1), ySTYLE = 6,$
     /noerase, YRANGE=[down_diff, up_diff]

    if up_diff-down_diff gt 80 then begin
        for i = -20., 200., 20. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -20., 200., 10. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 
    
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_diff, up_diff], $
                         ystyle=2, $  
                         YTITLE = 'TECU diference', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down_diff, up_diff], $
                         ystyle=2, $  
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]     

    up_E     = max(Ey)
    down_E   = min(Ey)
    plot, tot_days, Ey, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.07,0.18,0.95,0.30], $
     XSTYLE = 5, XRANGE=[0, tw],  XTICKNAME=REPLICATE(' ', tw+1), ySTYLE = 6,$
     /noerase, YRANGE=[down_E,up_E]  

    if up_E-down_E gt 30 then begin
        for i = -50., 50., 10. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -50., 50., 5. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 
    
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_E,up_E], $
                         ystyle=2, $   
                         YTITLE = 'Ey [mV/m]', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down_E,up_E], $ 
                         ystyle=2, $  
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]         	
    up_p    = max(p_dyn)
    down_p  = min(p_dyn) 
    plot, tot_days, p_dyn, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.07,0.05,0.95,0.17], $
     XSTYLE = 5, XRANGE=[0, tw], YRANGE=[down_p,up_p], XTICKNAME=REPLICATE(' ', tw+1), $
     ySTYLE = 6, /noerase;, SUBTITLE = time_title 

    if up_p-down_p gt 30 then begin
        for i = 0., 50., 10. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = 0., 50., 5. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 

        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         ;XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         xtitle=time_title,$
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                        ; xtitle=time_title,$
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_p,up_p], $
                         ystyle=2, $  
                         YTITLE = 'P [nPa]', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        ;YRANGE=[down,up]
                        
        AXIS, YAXIS = 1, YRANGE=[down_p,up_p], $
                         ystyle=2, $  
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    
;                         YRANGE=[down0,up0]   
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak   
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; open the post stript device
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'mag_tec_V1_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'mag_tec_V1_'+Date+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'mag_tec_V1_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	
                    
end


