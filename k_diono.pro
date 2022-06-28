;Name: k_diono.pro
;purpose:
;	this routine will generate a figure of geomagnetic disturbances and K index 
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
;   .r k_diono
;   k_diono, r_dst, k_dst, Bsq, DOY, idate[yyyy,mm,dd], fdate[yyyy,mm,dd]
;parameters:
;   magnetic indexes
;
;dependencies:
;   Instituto de Geofísica, Unidad Michoacan
;
;input files
;   Dst Index, DH index, Kp Index, Kmex index, ksigma[K1,K2]
;
;output files:
;   DP2, and K indexes
function dst_data, initial

	On_error, 2
	compile_opt idl2, HIDDEN

	year = string(initial, format = '(I4)')
		file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_D.dat'
		
        header = 25             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
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
        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
        DOY : 0, Dst: 0}

		r_dst = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], r_dst, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,8X,I3,X,I6)'
		
		return, r_dst
end

function teo, date
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
		
		file_name = '../rutidl/teoloyucan/hourly/'+'teo_'+date+'h.dat'
		
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
        DStruct = {H : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(F10)'		
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
        date = string(year, month, day, format = '(I4, I02, I02)')		
		name = 'teo_'+date+'.index.'		
		file_name = '../rutidl/Kmex/'+name+'final'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/Kmex/'+name+'early'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
	    
	    endif

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
        idx_kmex      = {k_mex      : intarr(8), $
                         k_mex_sum  : 0, $
                         a_mex      : intarr(8), $
                         a_med      : 0}
        
        struct = {x : [0, 0, 0, 0, 0, 0, 0, 0], y : 0}
        
        tmp_var = replicate(struct, 2)
  
		READS, data, tmp_var, FORMAT='(I3, I4, I4, I4, I4, I4, I4, I4, I4)'
		
		idx_kmex.k_mex[*]   = tmp_var[0].x
		idx_kmex.a_mex[*]   = tmp_var[1].x
        idx_kmex.k_mex_sum  = tmp_var[0].y
        idx_kmex.a_med      = tmp_var[1].y		
		
		return, idx_kmex	
end

function new_kmex, date
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
		
		name = 'new_km'+date+'.txt'
		
		file_name = '../rutidl/Kmex/new_kmex/'+name		
		file = FILE_SEARCH(file_name, COUNT=opened_files)	

	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/Kmex/new_kmex/'+name
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
	    
	    endif

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        idx_newkmex      = {n_kmex1      : intarr(8), $
                            n_kmex2      : intarr(8)}
        
        struct = {x : [0, 0, 0, 0, 0, 0, 0, 0]}
        
        tmp_var = replicate(struct, 2)  
		READS, data, tmp_var, FORMAT='(I3, I3, I3, I3, I3, I3, I3, I3)'
		
		idx_newkmex.n_kmex1[*]   = tmp_var[0].x
		idx_newkmex.n_kmex2[*]   = tmp_var[1].x		
		
		return, idx_newkmex	
end
function kp_data, in

	On_error, 2
	compile_opt idl2, HIDDEN

        header = 36             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
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
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
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
		resulting_data = REPLICATE(DataStruct2, number_of_lines-header)	

		resulting_data[*].year   = resulting_data0[*].year
		resulting_data[*].month  = resulting_data0[*].month
		resulting_data[*].day    = resulting_data0[*].day
		resulting_data[*].hour   = resulting_data0[*].hour
		resulting_data[*].minute = resulting_data0[*].minute
		resulting_data[*].DOY    = resulting_data0[*].DOY
		resulting_data[*].Kp     = Kp_tmp[*]
		resulting_data[*].Ap     = resulting_data0[*].AP				
		return, resulting_data
end    

function baseline_sq, date
	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr	= date[0]
	mh	= date[1]
	dy  = date[2]

        date = string(yr, mh, dy, format = '(I4,I02,I02)')
        header=0

		file_name = '../rutidl/output/Bsq_baselines/'+'Bsq_'+date+'h.dat
	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
		
        ;DStruct = {year : 0, month : 0, day : 0, hour : 0, doy : 0, Bsq : 0.}                                   
        DStruct = {Bsq : 0.}
		B_sq = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], B_sq, $
	 format='(F08.4)' 		
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
	return, DOY
	END
	
	
pro k_diono, r_dst, r_kp, DOY, date_i, date_f, JPEG = jpeg 
	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
;###############################################################################
idate0 = STRING(yr_i, mh_i, FORMAT='(I4,I02)')
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
    d_kp = kp_data(yr_i)     
    i_kp = d_kp.Kp       
;###############################################################################  
    d_dst = dst_data(yr_i)
    t = N_ELEMENTS(d_dst.year)    
    i_dst = d_dst.Dst
   
    year = d_dst.year
    tiempo = TIMEGEN(t, START=julday(d_dst.month[0], d_dst.day[0],  $
                     d_dst.year[0], d_dst.hour[0]), UNITS='Hours')  
                                        
        iyear = STRMID(STRING(yr_i, FORMAT='(I4)'),2,2)
        fyear = STRMID(STRING(yr_f, FORMAT='(I4)'),2,2)
                
        idoy      = Date2DOY(STRING(iyear, mh_i, dy_i,FORMAT = '(I02,I02,I02)'))
        fdoy      = Date2DOY(STRING(fyear, mh_f, dy_f,FORMAT = '(I02,I02,I02)'))         
            
    time_w  = tiempo[idoy:fdoy]
    tw      = N_ELEMENTS(time_w)
    tot_days= findgen(tw*24)/24.0
    
    dst     = i_dst[(idoy*24)-24:fdoy*24-1]   
    
    Kp      =  i_kp[(idoy*8)-8:fdoy*8-1]
    Kp      = Kp/10.0
    Date    = STRING(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
       file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
        string_date        = STRARR(file_number)    
                   
        data_file_name_h  = STRARR(file_number)
        data_file_name_km  = STRARR(file_number)
        data_file_name_kmn = STRARR(file_number)                
        data_file_name_bsq  = strarr(file_number)                               
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                
                data_file_name_h[i] = '../rutidl/teoloyucan/hourly/'+'teo_'+string_date[i]+'h.dat'
                data_file_name_km[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        data_file_name_km[i] = '../rutidl/Kmex/new_kmex/new_km'+string_date[i]+'.txt'
                data_file_name_bsq[i] = '../rutidl/output/Bsq_baselines/BsqV2_'+string_date[i]+'.txt'
                		        
		        file = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN BEGIN
	                data_file_name_km[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF 	                            
        ENDFOR

        exist_data_file_h   = FILE_TEST(data_file_name_h)
        capable_to_plot_h   = N_ELEMENTS(where(exist_data_file_h EQ 1))

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        exist_data_file_bsq   = FILE_TEST(data_file_name_bsq)
        capable_to_plot_bsq   = N_ELEMENTS(where(exist_data_file_bsq EQ 1))   
                
        IF capable_to_plot_h NE N_ELEMENTS(data_file_name_h) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF
                        
        H    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_h)-1 DO BEGIN
                IF exist_data_file_h[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_h = teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_h.H[*]
                                                                                                                       
                ENDIF ELSE BEGIN
                        H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
        
        k_mex    = FLTARR(file_number*8)
        a_mex    = INTARR(file_number*8)
        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_km[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'                 
                        d_km = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]/10.0
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF             
        ENDFOR

        new_kmex1    = FLTARR(file_number*8)
        new_kmex2    = FLTARR(file_number*8)
        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_km[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'                 
                        d_km = new_kmex([tmp_year, tmp_month, tmp_day])
                        
                        new_kmex1[i*8:(i+1)*8-1] = d_km.n_kmex1[*]/10.0
                        new_kmex2[i*8:(i+1)*8-1] = d_km.n_kmex2[*]/10.0
                ENDIF             
        ENDFOR
        
    k_days = findgen(tw*8)/8.0 
        i_nan1 = WHERE(H EQ 999999.0, ncount)
        i_nan2 = WHERE(H EQ 99999.0, n2count)
        
        prcent_nan = FLOAT(ncount+n2count)*100.0
        print,'porcentaje de valores NaN:', prcent_nan/n_elements(H),'%'
        
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 999999.0 THEN BEGIN
                H[WHERE(H[*] EQ 999999.0)] = !Values.F_NAN          
            endif
        ENDFOR
        
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 99999.0 THEN BEGIN
                H[WHERE(H[*] EQ 99999.0)] = !Values.F_NAN          
            endif
        ENDFOR               
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo       
    H_tmp   = H
    H_exist = WHERE(finite(H_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
   
    ; interpolate at the locations of the bad data using the good data    
    IF nbaddata GT 0 THEN H_tmp[baddata] = interpol(H_tmp[H_exist], H_exist, baddata)
    H = H_tmp      
;###############################################################################
    new_dstdays = findgen(tw*1440)/1440.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_dst = FLTARR(N_ELEMENTS(new_dstdays))     	       
;############################################################################### 
; define diurnal baseline  
        Bsq    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_bsq)-1 DO BEGIN
                IF exist_data_file_bsq[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        sqline = baseline_sq([tmp_year, tmp_month, tmp_day])
                        
                        Bsq[i*24:(i+1)*24-1] = sqline.Bsq[*]
                                                                                                                       
                ENDIF                
        ENDFOR
;###############################################################################
; define frequencies  
   mlat         = 28.06*!pi
   ld           = cos(mlat/180)
   p_a          = dst*ld
   baseline     = Bsq + p_a             
   diono   = H-baseline
  ; print, Bsq
   time = 3600.0
   fn      = FLOAT(1.0/(2.0*time)) ; frecuencia de Nyquist       
;###############################################################################
    i_diff = diono
    new_idiff = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_idiff  = INTERPOL(i_diff, N_ELEMENTS(new_dstdays))
    new_idiff = tmp_idiff  
;###############################################################################
; define pass band frequencies  
passband_l = idate0
CASE passband_l of
    '200311' : passband_l = 1.2e-5  ;23:08 hr
    '200411' : passband_l = 1e-5    ;27:46
    '200505' : passband_l = 1.2e-5  ;23:08
    '201503' : passband_l = 1.14e-5 ;24:21
    '201705' : passband_l = 0.9e-5  ;30:51
    '201709' : passband_l = 0.95e-5  ;23:08
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
passband_u = idate0
CASE passband_u of
    '200311' : passband_u = 3.2e-5  ;08:40
    '200411' : passband_u = 2.5e-5  ;11:06
    '200505' : passband_u = 2.8e-5  ;09:55
    '201503' : passband_u = 2.8e-5  ;09:55
    '201705' : passband_u = 3.2e-5  ;08:40
    '201709' : passband_u = 3.4e-5  ;08:10
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
; define high band frequencies
highpass_l = idate0
CASE highpass_l of
    '200311' : highpass_l = 7e-5    ;03:58
    '200411' : highpass_l = 7.5e-5  ;03:42
    '200505' : highpass_l = 7.3e-5  ;03:48
    '201503' : highpass_l = 7.5e-5  ;03:42
    '201705' : highpass_l = 7.5e-5  ;03:42
    '201709' : highpass_l = 7e-5    ;03:58
    ELSE: PRINT, 'fuera de rango'
ENDCASE
;###############################################################################
; define filtering DIGITAL_FILTER
    coeff_ddyn    = DIGITAL_FILTER(passband_l/fn, passband_u/fn, 50, 18)
    coeff_dp2   = DIGITAL_FILTER(highpass_l/fn, 1.0, 50, 4)
;###############################################################################
; define disturbing effects 
    ddyn        = convol(diono, coeff_ddyn, /edge_wrap)
    dp2         = convol(diono, coeff_dp2, /edge_wrap)  
;###############################################################################      
    new_ddyn = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_ddyn  = INTERPOL(ddyn, N_ELEMENTS(new_dstdays))
    new_ddyn = tmp_ddyn           
;############################################################################### 
    new_dp2 = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_dp2  = INTERPOL(dp2, N_ELEMENTS(new_dstdays))
    new_dp2 = tmp_dp2         
;###############################################################################
; define device and color parameters       
        Device_bak2 = !D.Name         
        SET_PLOT, 'Z'      
        
        Xsize=FIX(1600)
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
        azul      = 90
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
                        X_label[i]  = (tmp_month EQ old_month) ? STRING(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+STRING(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR 
        
CASE old_month of
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
    ELSE: PRINT, 'fuera de rango'
ENDCASE             
;###############################################################################
    spam_i = idate0
CASE spam_i of
    '200311' : spam_i = 1000
    '200411' : spam_i = 50
    '200505' : spam_i = 100
    '201503' : spam_i = 0
    '201705' : spam_i = 0
    '201709' : spam_i = 0
    ELSE: PRINT, 'fuera de rango'
ENDCASE 
;###############################################################################
    spam_f = idate0
CASE spam_f of
    '200311' : spam_f = 2800
    '200411' : spam_f = 4400
    '200505' : spam_f = 0
    '201503' : spam_f = 0
    '201705' : spam_f = 0
    '201709' : spam_f = 0
    ELSE: PRINT, 'fuera de rango'
ENDCASE    
;###############################################################################
       days = intarr(tw+1)
       FOR i=0, n_elements(days)-1 DO BEGIN
            days[i] = dy_i+i
       ENDFOR
       days = days*24/24. 
       day_time = findgen(24)   
;############################################################################### 
    time_title = ' Tiempo Universal ['+textoidl("dias")+' de '+old_month+'].'
    window_title = 'TGM'+ STRING(TGM_n, FORMAT='(I01)')+', '+ $
                    STRING(old_month, yr_i, FORMAT='(A, X, I4)')
;###############################################################################
    med_ddyn = MEDIAN(new_ddyn)
    std_ddyn = stddev(new_ddyn, /NAN)
    
    ddyn_out = WHERE(new_ddyn GE med_ddyn+std_ddyn OR new_ddyn LE med_ddyn-std_ddyn)
    ddyn_in  = WHERE(new_ddyn LE med_ddyn+std_ddyn AND new_ddyn GE med_ddyn-std_ddyn)
    
    ddyn_diff_out = new_ddyn
    ddyn_diff_out[ddyn_in]=!Values.F_NAN
    
    ddyn_diff_in  = new_ddyn
    ddyn_diff_in[ddyn_out]=!Values.F_NAN
     
     updp2     = max(diono-ddyn)
     downdp2   = min(diono-ddyn)     
                               
     PLOT, tot_days, dp2, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
     POSITION=[0.1,0.58,0.95,0.9], XSTYLE = 5, XRANGE=[0, tw], YSTYLE = 6,$
     XTICKNAME=REPLICATE(' ', tw+1), YRANGE=[downdp2,updp2], /NODATA          
;###############################################################################
    med_dp2 = MEDIAN(new_dp2)
    std_dp2 = STDDEV(new_dp2, /NAN)
    
    dp2_out = WHERE(new_dp2 GE med_dp2+std_dp2 OR new_dp2 LE med_dp2-std_dp2)
    dp2_in  = WHERE(new_dp2 LE med_dp2+std_dp2 AND new_dp2 GE med_dp2-std_dp2)
    
    dp2_diff_out = new_dp2
    dp2_diff_out[dp2_in]=!Values.F_NAN
    
    dp2_diff_in  = new_dp2
    dp2_diff_in[dp2_out]=!Values.F_NAN     
;###############################################################################
    dp2_i = idate0
case dp2_i of
    '200311' : dp2_i = dp2_out[6]
    '200411' : dp2_i = dp2_out[50]
    '200505' : dp2_i = dp2_out[100]
    '201503' : dp2_i = dp2_out[500]
    '201705' : dp2_i = dp2_out[100]
    '201709' : dp2_i = dp2_out[500]
    else: print, 'fuera de rango'
endcase  

    dp2_si = idate0
case dp2_si of
    '200311' : dp2_si = 0
    '200411' : dp2_si = -100
    '200505' : dp2_si = -100
    '201503' : dp2_si = -730
    '201705' : dp2_si = -50
    '201709' : dp2_si = 0
    else: print, 'fuera de rango'
endcase 

    dp2_sf = idate0
case dp2_sf of
    '200311' : dp2_sf = 20
    '200411' : dp2_sf = 100
    '200505' : dp2_sf = 350
    '201503' : dp2_sf = 550
    '201705' : dp2_sf = 900
    '201709' : dp2_sf = 990
    else: print, 'fuera de rango'
endcase 
;############################################################################### 
        OPLOT, new_dstdays ,new_idiff-new_ddyn, COLOR=negro, LINESTYLE=0, THICK=4          
        OPLOT, new_dstdays, new_dp2, COLOR=rojo

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
                         CHARSIZE = 0.9, $                                                
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[downdp2,updp2], $ 
                         YSTYLE=2, $  
                         YTITLE = '[nT]', $                          
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
                    POSITION=[0.1,0.1,0.95,0.51], XSTYLE = 5, YSTYLE = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], $
                    /NOERASE
        j = N_ELEMENTS(k_mex)
        FOR i = 0, j-1 DO BEGIN
                IF k_mex[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (k_mex[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex[i] EQ 4 : COLOR = azul
                                                k_mex[i] GE 4 : COLOR = azul
                                                ELSE       : COLOR = azul
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], COLOR=color
                                      ENDIF $
                                      ELSE BEGIN                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., LINESTYLE = 0
                                                                                  
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.],COLOR=morado, /LINE_FILL, $
                                         ORIENTATION=-45., LINESTYLE = 0
                                      ENDELSE
        ENDFOR

        XYOUTS, 0.01, 0.031 , /NORMAL, $
                'Codigo de Color:        Kp,             Kp y Kmex,           Kmex,                datos no disponibles.', COLOR=negro, $
                CHARSIZE = 0.9, $   
                CHARTHICK=chr_thick1
        POLYFILL, [0.09,0.12,0.12,0.09], [0.025,0.025,0.045,0.045], COLOR = verde, /NORMAL
        POLYFILL, [0.17,0.20,0.20,0.17], [0.025,0.025,0.045,0.045], COLOR = amarillo, /NORMAL
        POLYFILL, [0.27,0.3,0.3,0.27], [0.025,0.025,0.045,0.045], COLOR = azul, /NORMAL
        POLYFILL, [0.52,0.55,0.55,0.52], [0.025,0.025,0.045,0.045], COLOR = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., LINESTYLE = 0
        POLYFILL, [0.52,0.55,0.55,0.52], [0.025,0.025,0.045,0.045], COLOR = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., LINESTYLE = 0  
;###############################################################################                           
    PLOT, k_days, Kp, PSYM=6, /NODATA, MAX_VALUE=9., XTICKS=file_number, XMINOR=8, $
                    TITLE = '', SUBTITLE = '', XTITLE = 's', YTITLE = 's', $
                    BACKGROUND = blanco, COLOR=negro, YRANGE=[0,9], YTICKS=9, $
                    YMINOR=0, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                    POSITION=[0.1,0.1,0.95,0.51]  , XSTYLE = 5, YSTYLE = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number],$
                    /NOERASE
        j = N_ELEMENTS(Kp)
        FOR i = 0, j-1 DO BEGIN
                IF Kp[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (Kp[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                Kp[i] EQ 4 : COLOR = verde
                                                Kp[i] GE 4 : COLOR = verde
                                                ELSE       : COLOR = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,Kp[i]$
                                        +step,Kp[i]+step], COLOR=color
                                      ENDIF $
                                      ELSE BEGIN                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., LINESTYLE = 0
                                                                                  
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.],COLOR=morado, /LINE_FILL, $
                                         ORIENTATION=-45., LINESTYLE = 0
                                      ENDELSE
        ENDFOR
        
        FOR i = 0, j-1 DO BEGIN
                IF k_mex[i] LT Kp[i] THEN BEGIN
                                        color = 0
                                        step  = (k_mex[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex[i] EQ 4 : COLOR = azul
                                                k_mex[i] GE 4 : COLOR = azul
                                                ELSE       : COLOR = azul
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], COLOR=color
                                      ENDIF 
                                      IF  k_mex[i] EQ Kp[i] THEN BEGIN
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], COLOR=amarillo
                                        ENDIF                                      
        ENDFOR        
 
 FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,9.], LINESTYLE=1, COLOR=negro
 ENDFOR

        FOR i=5, 8, 1 DO OPLOT, [0,file_number], [i,i], LINESTYLE=1, COLOR=negro
;###############################################################################
FOR i=0, N_ELEMENTS(Kp)-1 DO BEGIN         
     IF Kp[i] GT k_mex[i] THEN BEGIN
        ERRPLOT,k_days[i]+1.5/24, new_kmex1[i], Kp[i], COLOR=negro, THICK=3
        ERRPLOT, k_days[i]+1.5/24, Kp[i], new_kmex2[i], COLOR=negro, LINESTYLE=0
    ENDIF  
        IF Kp[i] LT k_mex[i] THEN BEGIN
        ERRPLOT,k_days[i]+1.5/24, new_kmex1[i], Kp[i], COLOR=negro, LINESTYLE=0
        ERRPLOT, k_days[i]+1.5/24, Kp[i], new_kmex2[i], COLOR=negro, THICK=3
    ENDIF
ENDFOR     
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
                         YTITLE = 'indices K', $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.01

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.01;, $
;###############################################################################    
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.95   
   XYOUTS, X, y, window_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.65  

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.921, 'Tiempo Local', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=0.9  

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.53, 'Tiempo Local', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=0.9   
;###############################################################################                     
;first panel legend                   
        POLYFILL, [0.79,0.82,0.82,0.79], [0.654,0.654,0.657,0.657], color = rojo, /NORMAL
        POLYFILL, [0.79,0.82,0.82,0.79], [0.624,0.624,0.627,0.627], color = negro, /NORMAL        

        XYOUTS, 0.825, 0.65 , /NORMAL, $
                'DP2', COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1 
                
     D1 = TexToIDL('P_{PI}')
        XYOUTS, 0.825, 0.62 , /NORMAL, $
               D1+'-Ddyn', COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1                                  
;###############################################################################
; saving png     
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!    
    TVLCT, R_bak, G_bak, B_bak, /get  
        
    SET_PLOT, Device_bak2  
    path = '../rutidl/output/eventos_tgm/'
        IF keyword_set(jpeg) THEN BEGIN
                info = SIZE(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = R_bak[image]
                true_image[1,*,*] = G_bak[image]
                true_image[2,*,*] = B_bak[image]
                write_jpeg, path+'k_vs_diono_V4_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN PRINT, '        Setting PNG as default file type.'
                WRITE_PNG, path+'k_vs_diono_V4_'+Date+'.png', Image, R_bak, G_bak, B_bak
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                PRINT, '        Saving: '+path+'k_vs_diono_V4_'+Date+'.png'
                PRINT, ''
        ENDIF
        RETURN 	
END	
