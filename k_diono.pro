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

function kmex, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4, I02, I02)')
       ; sts = string(stat, format='(A5)')
		
		name = 'teo_'+date+'.index.'
		
		file_name = '../rutidl/Kmex/'+name+'final'
       ; file_name = '../rutidl/ip/'+yr+'-'+mt+'-'+dy+'.csv'		
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/Kmex/'+name+'early'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
	    
	    endif

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        idx_kmex      = {k_mex      : intarr(8), $
                         k_mex_sum  : 0, $
                         a_mex      : intarr(8), $
                         a_med      : 0}
        
        struct = {x : [0, 0, 0, 0, 0, 0, 0, 0], y : 0}
        
        tmp_var = replicate(struct, 2)

	;	idx_kmex = REPLICATE(DStruct, number_of_lines-header)	
  
		READS, data, tmp_var, FORMAT='(I3, I4, I4, I4, I4, I4, I4, I4, I4)'
		
		idx_kmex.k_mex[*]   = tmp_var[0].x
		idx_kmex.a_mex[*]   = tmp_var[1].x
        idx_kmex.k_mex_sum  = tmp_var[0].y
        idx_kmex.a_med      = tmp_var[1].y		
		
		return, idx_kmex	
end

function kp_data, in

	On_error, 2
	compile_opt idl2, HIDDEN

        header = 36             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        year = string(in, format = '(I4)')
		file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_D.dat'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)

	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_P.dat'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN begin
    	        file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_Q.dat'
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
		return, r_kp
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
	
	
pro k_diono, r_dst, r_kp, B_sq, DOY, date_i, date_f, JPEG = jpeg 

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
    d_kp = kp_data(yr_i)     
    i_kp = d_kp.Kp        
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
    
    Kp      =  i_dst[(idoy*8)-8:fdoy*8-1]
    Kp      = Kp/10.0
    Date    = string(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
;###############################################################################
       file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
        string_date        = STRARR(file_number)    
                   
        data_file_name_dh  = strarr(file_number)
        data_file_name_km  = strarr(file_number)
        
        string_date        = strarr(file_number)
        
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                
                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
                data_file_name_km[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        file = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name_km[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF 	                            
        ENDFOR

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))
        
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
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
        
        k_mex    = fltarr(file_number*8)
        a_mex    = INTARR(file_number*8)
        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_km[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'                 
                        d_km = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF             
        ENDFOR
    k_mex = k_mex/10.0
    k_days = findgen(tw*8)/8.0 
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
; define diurnal baseline
;###############################################################################  
    sqline = baseline_sq([yr_i, mh_i])
    
    sq_doy  = sqline.doy
    Bsq_ln   = sqline.Bsq


    tmp_doy = dst_doy[(idoy*24)-24:fdoy*24-1]
    Bsq     = fltarr(n_elements(Bsq_ln))
    
    tmp_doy = tmp_doy[uniq(tmp_doy, sort(tmp_doy))]
    
    for i=0, n_elements(tmp_doy)-1 do begin
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
    '200311' : passband_l = 1.2e-5
    '200411' : passband_l = 1e-5
    '200505' : passband_l = 1.2e-5
    '201503' : passband_l = 1.14e-5
    '201705' : passband_l = 0.9e-5
    '201709' : passband_l = 1.2e-5
    else: print, 'fuera de rango'
endcase  

passband_u = idate0
case passband_u of
    '200311' : passband_u = 3.2e-5
    '200411' : passband_u = 2.5e-5
    '200505' : passband_u = 2.8e-5
    '201503' : passband_u = 2.8e-5
    '201705' : passband_u = 3.2e-5
    '201709' : passband_u = 3.4e-5
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
        DEVICE, SET_RESOLUTION = [Xsize,Ysize],Set_Pixel_Depth=24, DECOMPOSED=1  
        DEVICE, z_buffer=4
        DEVICE, set_character_size = [10, 12] 
        
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 190
        verde     = 150
        negro     = 0
        azul      = 70
        blanco    = 255
        gris      = 110
        morado    = 16
        naranja  = 220
                
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT
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
    '200311' : spam_f = 2800
    '200411' : spam_f = 4400
    '200505' : spam_f = 1800
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
    time_title = ' Tiempo Universal ['+textoidl("dias")+' de '+old_month+'].'
    window_title = 'TGM'+ STRING(TGM_n, FORMAT='(I01)')+', '+ $
                    STRING(old_month, yr_i, FORMAT='(A, X, I4)')
;###############################################################################                                                           
;###############################################################################[0.1,0.73,0.95,0.9]
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
 
    ; upddyn     = max(ddyn)
    ; downddyn   = min(ddyn)
     
     updp2     = max(dp2)
     downdp2   = min(dp2)     

   ; IF upddyn GT updp2 THEN up = upddyn ELSE up=updp2 
   ; IF downddyn LT downdp2 THEN down = downddyn ELSE down=downdp2 
                               
     PLOT, tot_days, ddyn, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
     POSITION=[0.1,0.58,0.95,0.9], XSTYLE = 5, XRANGE=[0, tw], YSTYLE = 6,$
     XTICKNAME=REPLICATE(' ', tw+1), YRANGE=[downdp2,updp2], /NODATA
;###############################################################################
    ddyn_i = idate0
CASE ddyn_i of
    '200311' : ddyn_i = ddyn_out[0]
    '200411' : ddyn_i = ddyn_out[0]
    '200505' : ddyn_i = ddyn_out[100]
    '201503' : ddyn_i = ddyn_out[0]
    '201705' : ddyn_i = ddyn_out[0]
    '201709' : ddyn_i = ddyn_out[0]
    ELSE: PRINT, 'fuera de rango'
ENDCASE 

    ddyn_si = idate0
CASE ddyn_si of
    '200311' : ddyn_si = -200
    '200411' : ddyn_si = -100
    '200505' : ddyn_si = -160
    '201503' : ddyn_si = -100
    '201705' : ddyn_si = -130
    '201709' : ddyn_si = -180
    ELSE: PRINT, 'fuera de rango'
ENDCASE 

    ddyn_sf = idate0
CASE ddyn_sf of
    '200311' : ddyn_sf = -800
    '200411' : ddyn_sf = 300
    '200505' : ddyn_sf = 470
    '201503' : ddyn_sf = 170
    '201705' : ddyn_sf = 280
    '201709' : ddyn_sf = 350
    ELSE: PRINT, 'fuera de rango'
ENDCASE       
;###############################################################################      
       ; OPLOT, new_dstdays[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
;        ddyn_diff_out[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
 ;       COLOR=negro, LINESTYLE=0, THICK=5   
                    
  ;      OPLOT, new_dstdays[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
   ;     ddyn_diff_in[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
    ;    COLOR=negro, LINESTYLE=0, THICK=5            

     ;   OPLOT, new_dstdays, ddyn_diff_in, COLOR=negro, LINESTYLE=3
      ;  OPLOT, new_dstdays, ddyn_diff_out, COLOR=negro, LINESTYLE=0, THICK=5
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
;###############################################################################
;###############################################################################
    dp2_i = idate0
case dp2_i of
    '200311' : dp2_i = dp2_out[6]
    '200411' : dp2_i = dp2_out[50]
    '200505' : dp2_i = dp2_out[0]
    '201503' : dp2_i = dp2_out[0]
    '201705' : dp2_i = dp2_out[0]
    '201709' : dp2_i = dp2_out[350]
    else: print, 'fuera de rango'
endcase  

    dp2_si = idate0
case dp2_si of
    '200311' : dp2_si = 0
    '200411' : dp2_si = -100
    '200505' : dp2_si = -200
    '201503' : dp2_si = -100
    '201705' : dp2_si = -100
    '201709' : dp2_si = -100
    else: print, 'fuera de rango'
endcase 

    dp2_sf = idate0
case dp2_sf of
    '200311' : dp2_sf = 20
    '200411' : dp2_sf = 100
    '200505' : dp2_sf = 50
    '201503' : dp2_sf = 100
    '201705' : dp2_sf = 130
    '201709' : dp2_sf = -400
    else: print, 'fuera de rango'
endcase      
;############################################################################### 
;###############################################################################     
        OPLOT, tot_days, dp2, COLOR=rojo
        OPLOT, new_dstdays[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], $
        new_dp2[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], COLOR=rojo, $
        LINESTYLE=0, THICK=4      
        
        OPLOT, new_dstdays[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], $
        new_dp2[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], COLOR=rojo, $
        LINESTYLE=0, THICK=4            
;############################################################################### 
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
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $  
                         CHARSIZE = 0.8, $                                                
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[downdp2,updp2], $ 
                         YSTYLE=2, $  
                         YTITLE = 'DP2 [nT]', $                          
                         COLOR=negro, $
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[downdp2,updp2], $ 
                         YSTYLE=2, $                           
                         COLOR=negro, $
                         CHARSIZE = 0.9;, $      

;###############################################################################        
    PLOT, k_days, k_mex, PSYM=6, /NODATA, MAX_VALUE=9., XTICKS=file_number, XMINOR=8, $
                    TITLE = '', SUBTITLE = '', XTITLE = 's', YTITLE = 's', $
                    BACKGROUND = blanco, COLOR=negro, YRANGE=[0,9], YTICKS=9, $
                    YMINOR=0, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                    POSITION=[0.1,0.34,0.95,0.51], XSTYLE = 5, YSTYLE = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], $
                    /NOERASE
        j = N_ELEMENTS(k_mex)
        FOR i = 0, j-1 DO BEGIN
                IF k_mex[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (k_mex[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex[i] EQ 4 : COLOR = amarillo
                                                k_mex[i] GE 4 : COLOR = rojo
                                                ELSE       : COLOR = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], COLOR=color
                                      ENDIF $
                                      ELSE BEGIN                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., linestyle = 0
                                                                                  
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.],COLOR=morado, /LINE_FILL, $
                                         ORIENTATION=-45., linestyle = 0
                                      ENDELSE
        ENDFOR
 
 FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,9.], linestyle=1, COLOR=negro
        ENDFOR

        FOR i=5, 8, 1 DO OPLOT, [0,file_number], [i,i], linestyle=1, COLOR=negro
        ;OPLOT, [0,file_number], [6.,6.], linestyle=1, COLOR=negro

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTITLE=time_title, $ 
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $ 
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTITLE = 'Kmex', $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $

        XYOUTS, 0.01, 0.031 , /NORMAL, $
                'Codigo de Color:        quieto,           perturbado,         tormenta,            datos no disponibles.', COLOR=negro, $
                CHARSIZE = 0.9, $   
                CHARTHICK=chr_thick1
        POLYFILL, [0.09,0.12,0.12,0.09], [0.025,0.025,0.045,0.045], color = verde, /NORMAL
        POLYFILL, [0.17,0.20,0.20,0.17], [0.025,0.025,0.045,0.045], color = amarillo, /NORMAL
        POLYFILL, [0.27,0.3,0.3,0.27], [0.025,0.025,0.045,0.045], color = rojo, /NORMAL
        POLYFILL, [0.52,0.55,0.55,0.52], [0.025,0.025,0.045,0.045], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., linestyle = 0
        POLYFILL, [0.52,0.55,0.55,0.52], [0.025,0.025,0.045,0.045], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., linestyle = 0  
;###############################################################################                    
;###############################################################################       
    PLOT, k_days, k_mex, PSYM=6, /NODATA, MAX_VALUE=9., XTICKS=file_number, XMINOR=8, $
                    TITLE = '', SUBTITLE = '', XTITLE = 's', YTITLE = 's', $
                    BACKGROUND = blanco, COLOR=negro, YRANGE=[0,9], YTICKS=9, $
                    YMINOR=0, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                    POSITION=[0.1,0.1,0.95,0.27]  , XSTYLE = 5, YSTYLE = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number],$
                    /NOERASE
        j = N_ELEMENTS(k_mex)
        FOR i = 0, j-1 DO BEGIN
                IF k_mex[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (k_mex[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex[i] EQ 4 : COLOR = amarillo
                                                k_mex[i] GE 4 : COLOR = rojo
                                                ELSE       : COLOR = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], COLOR=color
                                      ENDIF $
                                      ELSE BEGIN                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., linestyle = 0
                                                                                  
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.],COLOR=morado, /LINE_FILL, $
                                         ORIENTATION=-45., linestyle = 0
                                      ENDELSE
        ENDFOR
 
 FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,9.], linestyle=1, COLOR=negro
        ENDFOR

        FOR i=5, 8, 1 DO OPLOT, [0,file_number], [i,i], linestyle=1, COLOR=negro
        ;OPLOT, [0,file_number], [6.,6.], linestyle=1, COLOR=negro

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTITLE=time_title, $ 
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=0.9, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $ 
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTITLE = 'Kp', $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $

;###############################################################################
;###############################################################################          
;###############################################################################    
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.95   
   XYOUTS, X, y, window_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.65  

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.921, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9  

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.53, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9  
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.29, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9   
;###############################################################################     
    ;    POLYFILL, [0.88,0.91,0.91,0.88], [0.405,0.405,0.407,0.407], color = rojo, /NORMAL
    ;    POLYFILL, [0.88,0.91,0.91,0.88], [0.375,0.375,0.377,0.377], color = negro, /NORMAL        

    ;    XYOUTS, 0.91, 0.4 , /NORMAL, $
    ;            ' DP2', COLOR=negro, $
    ;            CHARSIZE = chr_size1, $
    ;            CHARTHICK=chr_thick1 

    ;    XYOUTS, 0.91, 0.37 , /NORMAL, $
    ;            ' Ddyn', COLOR=negro, $
    ;            CHARSIZE = chr_size1, $
    ;            CHARTHICK=chr_thick1                                    
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
                write_jpeg, path+'k_vs_diono_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'k_vs_diono_'+Date+'.png', Image, R_bak, G_bak, B_bak
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'k_vs_diono_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	
end	
