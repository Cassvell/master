;
;Name:
;	magdat_v2
;purpose:
;	this routine will call read a certain number of files containing magnetic 
;   data index.
;author:
;	Carlos Isaac Castellanos Velazco
;	Estudiante de Maestría en Ciencias de la Tierra
;	Instituto de Geofísica, Unidad Michoacan
;	UNAM
;	ccastellanos@igeofisica.unam.mx
;
;category:
;   Getting the data magnetic, date and time within several files. 
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

function Date2DOY, idate
;	------------------------------------------------------------
;+							18-Sep-91
;	NAME: 
;		Date2DOY
;	PURPOSE:
;		Convert yymmdd into DOY (day of the year).  Input can
;		be either string or integer.  The year is an Optional 
;		return parameter.
;	CALLING SEQUENCE:
;		Date2DOY, idate, DOY [, yr]
;	INPUT:
;		idate	input format for the date: yymmdd.
;			Data-type can be either string or integer.
;	OUTPUT:
;		DOY	integer with the day of the year.
;	Output/Optional:
;		yr	year of the returned DOY.
;	Note:	If input data-type is string the returned values are
;		string-type and if input-type is longword the returned
;		parameters (DOY and yr) are integers.
;	HISTORY:
;		written by GAL 18-Sep-91
;-
;	-----------------------------------------------------------------
;	ON_ERROR, 2	;force a return to caller on error

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

function sym_data, date_i, date_f
	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
	
	header = 25
		
        initial_date = string(yr_i, mh_i, dy_i, format = '(I4,"-", I02,"-", I02)')
        final_date = string(yr_f, mh_f, dy_f, format = '(I4,"-", I02,"-", I02)')
        
        
		
		file_name = '../rutidl/sym/'+'ASY_'+initial_date+'_'+final_date+$
		'h_P'+'.dat'
		
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
        DStruct = {year : 0, month : 0, day : 0, hour : 0,  DOY : 0, ASY_D : 0, $
        ASY_H : 0, SYM_D : 0, SYM_H : 0}

		sym_mag = REPLICATE(DStruct, number_of_lines-header)	
  
		READS, data[header:number_of_lines-1], sym_mag, $
		FORMAT='(I4,X,I02,X,I02,X,I02,11X,I03,4X,I03,4X,I03,4X,I03,3X,I04)'
		
		return, sym_mag
end

function tec_data, i_date, f_date

	On_error, 2
	compile_opt idl2, HIDDEN

	iyear	= i_date[0]
	imonth	= i_date[1]
	iday 	= i_date[2]	

	fyear	= f_date[0]
	fmonth	= f_date[1]
	fday 	= f_date[2]	

        header = 1      ; Defining number of lines of the header 
	;	file = DIALOG_PICKFILE(FILTER='*.dat')

;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

        idate = string(iyear, imonth, iday, format = '(I4, "-", I02, "-", I02)')
        fdate = string(fyear, fmonth, fday, format = '(I4, "-", I02, "-", I02)')

		file_name = '../rutidl/tec/tec_newformat/'+'tec_'+idate+'_'+fdate+'.txt'

		
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
        DStruct = {doy : 0., tec : 0., med : 0.}                                   

		r_tec = REPLICATE(DStruct, number_of_lines-header)	
        ;print, number_of_lines-header-1, number_of_lines-header+1

       
		READS, data[header:number_of_lines-1], r_tec, $
	format='(F5, X, F5, X, F6)'	
	
		return, r_tec

end


function list_dst, r_dst, initial
    
    On_error, 2
	compile_opt idl2, HIDDEN

;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;calling the get_data_date function
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
	a = initial
	dat = dst_data(a) 

;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Defining variables to write output files indicating the events where Dst index
; got lower than -150 nT.
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
	tmp = n_elements(dat.year)
	
	indexes = WHERE(dat.Dst LE -150, count)
	IF count LE 0 THEN MESSAGE, 'ERROR'
		
	dst0      = dat.Dst
    doy0     = dat.DOY
    year0    = dat.year
    month0   = dat.month
    day0     = dat.day
    hour0    = dat.hour 
    
    idx_dst  = dst0[indexes]
	doy     = doy0[indexes]
	year    = year0[indexes]
	month   = month0[indexes]
	day     = day0[indexes]
	hour    = hour0[indexes]
    
	idx_dst_tmp   = idx_dst*0
	doy_tmp      = doy*0
	year_tmp     = year*0
	month_tmp    = month*0
	day_tmp      = day*0
	hour_tmp     = hour*0   

    j=0l
    idx_dst_tmp[j] = idx_dst_tmp[0]
	year_tmp[j]   = year[0]
	month_tmp[j]  = month[0]
    day_tmp[j]    = day[0]
	hour_tmp[j]   = hour[0]
    doy_tmp[j]    = doy[0]   
    
    Date    = strmid(string(dat.year[0]), 8, 4)
    outfile = '../rutidl/output/'+'dst_'+Date[0]+'TGM.txt'   ; defining the output file
    ;names and path directories

    FOR i=1l, N_ELEMENTS(idx_dst)-1 DO BEGIN
	IF doy_tmp[j] EQ doy[i] THEN BEGIN
	        IF idx_dst[i] LT idx_dst_tmp[j] THEN BEGIN
        	        idx_dst_tmp[j] = idx_dst[i]
		            year_tmp[j]   = year[i]
	                month_tmp[j]  = month[i]
        	        day_tmp[j]    = day[i]
	                hour_tmp[j]   = hour[i]
	                doy_tmp[j]    = doy[i]
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
	        idx_dst_tmp[j] = idx_dst[i]
	        year_tmp[j]   = year[i]
	        month_tmp[j]  = month[i]
        	day_tmp[j]    = day[i]
	        hour_tmp[j]   = hour[i]
	        doy_tmp[j]    = doy[i]
	ENDELSE
    ENDFOR
   
 ;   OPENW, lun, Outfile, /GET_LUN
   
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice Dst descendió por'
    print, 'debajo de -150 nT, lo que es un indicativo de la ocurrencia de un '
    print, 'un evento de tormenta geomagnética intensa y de interés  '
    print, '                                                                   '
    print, format='("Fecha", 6X, "Hora", X, "DOY", 8X, "Indice Dst")'
   
    for i=0, N_ELEMENTS(indexes)-1 do begin

           if doy_tmp[i] ne 0 then begin

            ;idx_kp[idx] = idx_kp[i]

                print, year_tmp[i], month_tmp[i], day_tmp[i], hour_tmp[i], doy_tmp[i], idx_dst_tmp[i], $
                             FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 2X, I03, X, I)'  
               ; printf, lun, year_tmp[i], month_tmp[i], day_tmp[i], hour_tmp[i], doy_tmp[i], idx_dst_tmp[i], $
                ;            FORMAT = '(I4, X, I02, X, I02, 2X, I02, 2X, I03, X, I4)'             
            endif 
   
        ;printf, lun, year[i], month[i], day[i], hour[i], doy[i], idx_kp[i], $
        ;FORMAT = '(I4, "/", I02, "/", I02, X, I02, 2X, I03, 4X, I)'       
    endfor
    
  ;  close,lun
   ; FREE_LUN, lun
   ;PRINT, DATETIME[i], Dst[i], format = "(A29)"
    print, '                                                                    '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento'     
end 

pro list_dh, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]	
;##############################################################################
; reading data files
;##############################################################################
        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        data_file_name = strarr(file_number)
        string_date     = strarr(file_number)
       
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                data_file_name[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
        ENDFOR

        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        H    = FLTARR(file_number*24)
                      
        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        dat = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = dat.H[*]
                ENDIF ELSE BEGIN
                         H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
        
        for i=0, n_elements(H)-1 do begin
            ;idx = where(H eq 999999.0, inan)
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor

    time_w = timegen(n_elements(file_number), final=julday(mh_f, dy_f, yr_f, 23), $
                start=julday(mh_i, dy_i, yr_i, 0) , units='H')

    caldat, time_w, m, d, y, hr
;   print, m, d, format='(I02, X, I02)'

    y_idx    = y[idx]   
    m_idx    = m[idx]
    d_idx    = d[idx]
    H_idx    = H[idx]
    hr_idx   = hr[idx]
    
    y_idx_tmp    = y_idx*0
    m_idx_tmp    = m_idx*0
    d_idx_tmp    = d_idx*0
    H_idx_tmp    = H_idx*0
    hr_idx_tmp   = hr_idx*0
;print, H[idx]
    j=0l
    y_idx_tmp[j]    = y_idx[0]
    m_idx_tmp[j]    = m_idx[0]
    d_idx_tmp[j]    = d_idx[0]
    H_idx_tmp[j]    = H_idx[0]
    hr_idx_tmp[j]   = hr_idx[0]
    
    FOR i=1l, N_ELEMENTS(H_idx)-1 DO BEGIN
	IF d_idx_tmp[j] EQ d_idx[i] THEN BEGIN
	        IF H_idx[i] GT H_idx_tmp[j] THEN BEGIN
        	        H_idx_tmp[j]   = H_idx[i]
		            y_idx_tmp[j]   = y_idx[i]
	                m_idx_tmp[j]   = m_idx[i]
        	        d_idx_tmp[j]   = d_idx[i]
	                hr_idx_tmp[j]  = hr_idx[i]	               
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
        	    H_idx_tmp[j]   = H_idx[i]
		        y_idx_tmp[j]   = y_idx[i]
	            m_idx_tmp[j]   = m_idx[i]
        	    d_idx_tmp[j]   = d_idx[i]
	            hr_idx_tmp[j]  = hr_idx[i]	
	ENDELSE
    ENDFOR
    
    outfile='../rutidl/output/'+'DH_list_TGM.txt'
        OPENW, lun, outfile, /GET_LUN, /append
        
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice DH descendió por'
    print, 'debajo de -100 nT, significando un evento de tormenta geomagnética '
    print, 'intensa.'
    print, '                                                                    '    
    print, format='(2X, "Fecha", 4X, "Hora", 2X, "Indice DH")'
    print, '---------------------------------'
    print, '                                                                    '
         
    for i=0, N_ELEMENTS(idx)-1 do begin
           if d_idx_tmp[i] ne 0 then begin

                print, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], hr_idx_tmp[i], $
                H_idx_tmp[i], $
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I04)' 
  
                printf, lun, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], $
                hr_idx_tmp[i], H_idx_tmp[i], $ 
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I04)'                         
        endif        
    endfor
    
    close,lun
    FREE_LUN, lun
   ;PRINT, DATETIME[i], Dst[i], format = "(A29)"
    print, '                                                                   '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento'
 
end
  
pro magdat_v2, r_dst, r_kp, r_tec, sym_mag, initial, doy_i, doy_f

	On_error, 2
	compile_opt idl2, HIDDEN
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Generate the time variables to plot time series of Dst Index
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
    yr  = initial
    ;dt  = tp  

    d_dst = dst_data(yr)
    t = n_elements(d_dst.year)
    ;print, t
    
    i_dst = d_dst.Dst
   
    year = d_dst.year
    tiempo = TIMEGEN(t, START=julday(d_dst.month[0], d_dst.day[0],  $
                     d_dst.year[0], d_dst.hour[0]), UNITS='Hours')
        
    ini_time  = JULDAY(1, doy_i, initial)
    fnl_time  = JULDAY(1, doy_f, initial)
    
    time_w  = tiempo[doy_i:doy_f]
    tw = n_elements(time_w)
    days_dst= findgen(tw*24)/24.0
    
    dst = i_dst[(doy_i*24)-24:doy_f*24-1]
                                    
    caldat, ini_time, ini_month, ini_day, ini_year
    caldat, fnl_time, fnl_month, fnl_day, fnl_year
    
    Date = string(year[0], ini_month, ini_day, FORMAT='(I4, "-", I02, "-", I02)')
    
;##############################################################################
; defining TEC and SYM data time & data variables
;############################################################################## 
	yr_i	= ini_year
	mh_i	= ini_month
	dy_i 	= ini_day	

	yr_f	= fnl_year
	mh_f	= fnl_month
	dy_f 	= fnl_day
    
    t_data = tec_data ([yr_i,mh_i,dy_i], [yr_f,mh_f,dy_f])
    
    DOY     = t_data.DOY
    fn = n_elements(DOY)
;###############################################################################   
    ;tec_doy = t_data.doy
    tec     = t_data.tec
    med     = t_data.med
    ;dif_tec = tec-med
    tec_days= findgen(tw*12)/12.0  
;###############################################################################     
    d_sym = sym_data([yr_i,mh_i,dy_i], [yr_f,mh_f,dy_f])
    
    sym_H = d_sym.SYM_H 
;##############################################################################
; defining Kp and time data variables
;##############################################################################   
    d_kp = kp_data(yr)
    t = n_elements(d_kp.year)

    i_kp = d_kp.Kp
    i_ap = d_kp.ap
;###############################################################################          
    k_days= findgen(tw*8)/8.0    
    Kp      = i_kp[(doy_i*8)-7:doy_f*8]
    ap      = i_ap[(doy_i*8)-7:doy_f*8]       
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Generate the time variables to plot time series of DH Index
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        
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

        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        H    = FLTARR(file_number*24)
        D    = FLTARR(file_number*24)
        Z    = FLTARR(file_number*24)
        F    = FLTARR(file_number*24)
        N    = FLTARR(file_number*24)
        
        H_STDESV    = FLTARR(file_number*24)
        D_STDESV    = FLTARR(file_number*24)
        Z_STDESV    = FLTARR(file_number*24) 
        F_STDESV    = FLTARR(file_number*24)
        N_STDESV    = FLTARR(file_number*24)
                       
        FOR i = 0, N_ELEMENTS(exist_data_file_dh)-1 DO BEGIN
                IF exist_data_file_dh[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_dh = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_dh.H[*]
                        D[i*24:(i+1)*24-1] = d_dh.D[*]
                        Z[i*24:(i+1)*24-1] = d_dh.Z[*]
                        F[i*24:(i+1)*24-1] = d_dh.F[*]
                        N[i*24:(i+1)*24-1] = d_dh.N[*]
                                                
                        H_STDESV[i*24:(i+1)*24-1] = d_dh.H_stdesv[*]
                        D_STDESV[i*24:(i+1)*24-1] = d_dh.D_stdesv[*]
                        Z_STDESV[i*24:(i+1)*24-1] = d_dh.Z_stdesv[*]  
                        F_STDESV[i*24:(i+1)*24-1] = d_dh.F_stdesv[*]
                        N_STDESV[i*24:(i+1)*24-1] = d_dh.N_stdesv[*]
                                                                                              
                ENDIF ELSE BEGIN
                         H[i*24:(i+1)*24-1] = 999999.0
                         D[i*24:(i+1)*24-1] = 999999.0
                         Z[i*24:(i+1)*24-1] = 999999.0
                         F[i*24:(i+1)*24-1] = 999999.0
                         N[i*24:(i+1)*24-1] = 999999.0

                         N_STDESV[i*24:(i+1)*24-1] = 999999.0                        
                         N_STDESV[i*24:(i+1)*24-1] = 999999.0
                         N_STDESV[i*24:(i+1)*24-1] = 999999.0
                         N_STDESV[i*24:(i+1)*24-1] = 999999.0
                         N_STDESV[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
        
        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor    
        
        for i=0, n_elements(H)-1 do begin
            ;idx = where(H eq 999999.0, inan)
            if H[i] ge 100.0 then begin
                H[where(H[*] ge 100.0)] = !Values.F_NAN          
            endif
        endfor  
;##############################################################################
; defining Kmex and time data variables
;############################################################################## 
	
        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        k_mex    = fltarr(file_number*8)
        a_mex    = INTARR(file_number*8)

        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_km[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        d_km = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF 
        ENDFOR     
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; initiate the figure Device, defining colors and figure dim
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        IF tw GT 7 THEN tmp_spam = 1.5
        IF tw GT 15 THEN tmp_spam = 2.
        
        Xsize=fix(800*tmp_spam)
        Ysize=1200
        ;DEVICE, decompose=0
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
             
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 198
        verde     = 160
        negro     = 0
        azul      = 80
        blanco    = 255
        gris      = 130
        morado    = 16
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT

    path = '../rutidl/output/globfig_to_reg/'   
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Create a time label based on the DOY initial and DOY final inputs
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
     X_label   = STRARR(tw+1)+' '
    ; print, n_elements(x_label)
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = ini_month
       ; print, old_month
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(1, doy_i, initial)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR    
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; plot the Dst time series for each event
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;file_name = ''
;IF file_name eq '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_D.dat' then plot_title = 'Definitive Dst index'
 
;if file_name eq '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_P.dat' then plot_title = 'Provisional Dst index'

;if file_name eq '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_Q.dat' then plot_title = 'Quick look Dst index'

window_title = 'from '+string(yr_i, mh_i, dy_i, $
                FORMAT='(I4, "/", I02, "/", I02)')+' to '$
                +string(yr_f, mh_f, dy_f, $
                FORMAT='(I4, "/", I02, "/", I02)')

time_title = 'Time [UT]' 
                
MAG_source = 'Source: International Service of Geomagnetic Indices'  

    if max(H) gt max(dst) then up0 = max(H) else up0 = max(dst)
    if min(H) lt min(dst) then down0 = min(H) else down0 = min(dst)
            
    plot, days_dst, H, XTICKS=tw, xminor = 8, POSITION=[0.07,0.66,0.95,0.96],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down0,up0],$
    title = window_title, BACKGROUND = blanco, COLOR=negro, XRANGE=[0, tw],$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', tw+1)

    oplot, days_dst, dst, COLOR=azul
  
    if up0-down0 gt 300 then begin
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
                         XTICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         XTICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0], $
                         YTITLE = 'Dst and DH [nT]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00

                        
        AXIS, YAXIS = 1, YRANGE=[down0,up0], $
                         ystyle=2, $
                         COLOR=negro, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    

     
    up      = max(tec)
    down    = min(tec)
    
    plot, tec_days, tec, XTICKS=tw, xminor = 8, POSITION=[0.07,0.04,0.95,0.34],$[0.07,0.04,0.95,0.34]
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down, up],$
    BACKGROUND = blanco, COLOR=negro, XRANGE=[0, tw], $
	XTICKNAME=REPLICATE(' ', tw+1), /noerase

    oplot, tec_days, med, COLOR=rojo
    
    if up-down gt 120 then begin
        for i = -100., 260., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -100., 260., 10. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 
         
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         ;XTICKFORMAT='(A1)',$
                         xtitle='time [UT]', $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$                         
                         ;XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down,up], $
                         YTITLE = 'TEC [TECu]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.6;, $
                        ; CHARTHICK=chr_thick1;, $    

if max(ap)gt max(a_mex) then up1 = max(ap) else up1=max(a_mex)
if min(ap)gt min(a_mex) then down1 = min(ap) else down1=min(a_mex)
    plot, k_days, ap, XTICKS=tw, xminor=8, XTITLE = 's', YTITLE = 's', /noerase, $
    XRANGE=[0, tw], BACKGROUND = blanco, COLOR=negro, CHARSIZE = chr_size1, $
    CHARTHICK=chr_thick1, POSITION=[0.07,0.35,0.95,0.65], XSTYLE = 5, ystyle = 6,$
                    XTICKNAME=REPLICATE(' ', tw+1), yrange=[down1,up1] 
    oplot, k_days, a_mex, color=verde

    if up1-down1 gt 300 then begin
        for i = -500., 500., 100. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -300., 300., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
    endelse 
    
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKFORMAT="(A1)",$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         ;CHARSIZE = 0.7, $
                         ;CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', tw+1), $
                         COLOR=negro, $
                         TICKLEN=0.04
                         


        AXIS, YAXIS = 0, YRANGE=[down1,up1], $
                         ;YTICKS=9, $
                        ; YMINOR=1, $
                         ystyle=2,$
                         YTITLE = 'Ap and Amex index', $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
                         CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00

        AXIS, YAXIS = 1, YRANGE=[down1,up1], $
                         ;YTICKS=9, $
                         ystyle=2,$
                         ;YMINOR=1, $
                         ;YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
                         CHARTHICK=chr_thick1;, $

;first panel legend
        POLYFILL, [0.77,0.80,0.80,0.77], [0.674,0.674,0.676,0.676], color = negro, /NORMAL
        POLYFILL, [0.86,0.89,0.89,0.86], [0.674,0.674,0.676,0.676], color = azul, /NORMAL        
    if tw gt 7 then begin
        XYOUTS, 0.777, 0.670 , /NORMAL, $
                '    Dst,          DH', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1         
    endif else begin
        XYOUTS, 0.764, 0.671 , /NORMAL, $
                '    Dst,      DH', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1     
    endelse

;second panel legend
        POLYFILL, [0.77,0.80,0.80,0.77], [0.62,0.62,0.622,0.622], color = negro, /NORMAL
        POLYFILL, [0.86,0.89,0.89,0.86], [0.62,0.62,0.622,0.622], color = verde, /NORMAL        
    if tw gt 7 then begin
        XYOUTS, 0.772, 0.617 , /NORMAL, $
                '     Ap,           Amex', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1         
    endif else begin
        XYOUTS, 0.772, 0.617 , /NORMAL, $
                '    Ap,      Amex', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1     
    endelse

;third panel legend

        POLYFILL, [0.70,0.73,0.73,0.70], [0.304,0.304,0.306,0.306], color = negro, /NORMAL
        POLYFILL, [0.83,0.86,0.86,0.83], [0.304,0.304,0.306,0.306], color = rojo, /NORMAL        
    if tw gt 7 then begin
        XYOUTS, 0.732, 0.302 , /NORMAL, $
                'Tec obs,             Tec med', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1         
    endif else begin
        XYOUTS, 0.731,   0.302 , /NORMAL, $
                'Tec obs,     Tec med', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1     
    endelse
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
                write_jpeg, path+'idx_V2'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'idx_V2'+Date+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'mag_idx_'+Date+'_V2.png'
                print, ''
        ENDIF
        RETURN
end

pro magdat_v3, r_dst, r_kp, r_tec, sym_mag, doy, idate, fdate

	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr_i	= idate[0]
	mh_i	= idate[1]
	dy_i 	= idate[2]	

	yr_f	= fdate[0]
	mh_f	= fdate[1]
	dy_f 	= fdate[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Generate the time variables to plot time series of Dst Index
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
    d_dst = dst_data(yr_i)
    t = n_elements(d_dst.year)
    ;print, t
    
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
    tw2     = tw+0.25
    days_dst= findgen(tw*24)/24.0
    
    dst     = i_dst[(idoy*24)-24:fdoy*24-1]
    
    Date    = string(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
    
;##############################################################################
; defining TEC and SYM data time & data variables
;##############################################################################     
    t_data = tec_data ([yr_i,mh_i,dy_i], [yr_f,mh_f,dy_f])
    
    DOY     = t_data.DOY
    fn = n_elements(DOY)
;###############################################################################   
    tec     = t_data.tec
    med     = t_data.med
    dif_tec = tec-med
    tec_days= findgen(tw*12)/12.0  
;###############################################################################     
    d_sym = sym_data([yr_i,mh_i,dy_i], [yr_f,mh_f,dy_f])
    
    sym_H = d_sym.SYM_H 
;##############################################################################
; defining Kp and time data variables
;##############################################################################   
    d_kp = kp_data(yr_i)
    t = n_elements(d_kp.year)

    i_kp = d_kp.Kp
    i_ap = d_kp.ap
;###############################################################################          
    k_days= findgen(tw*8)/8.0    
    Kp      = i_kp[(idoy*8)-8:fdoy*8]
    ap      = i_ap[(idoy*8)-8:fdoy*8]       
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Generate the time variables to plot time series of DH Index
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        
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

        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

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
        
        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor    
        
        for i=0, n_elements(H)-1 do begin
            ;idx = where(H eq 999999.0, inan)
            if H[i] ge 100.0 then begin
                H[where(H[*] ge 100.0)] = !Values.F_NAN          
            endif
        endfor  
;##############################################################################
; defining Kmex and time data variables
;############################################################################## 
	
        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        k_mex    = fltarr(file_number*8)
        a_mex    = INTARR(file_number*8)

        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_km[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        d_km = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF 
        ENDFOR     
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; initiate the figure Device, defining colors and figure dim
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        ;IF tw GT 7 THEN tmp_spam = 1.5
        ;IF tw GT 15 THEN tmp_spam = 2.
        
        Xsize=fix(1000*tmp_spam)
        Ysize=200
        ;DEVICE, decompose=0
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
             
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 198
        verde     = 160
        negro     = 0
        azul      = 80
        blanco    = 255
        gris      = 5
        morado    = 50
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT

    path = '../rutidl/output/eventos_tgm/'   
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Create a time label based on the DOY initial and DOY final inputs
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
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
    
    dy_i2 = dy_i-1    
    X_label2   = STRARR(tw+1)+' '
    ; print, n_elements(x_label)
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = mh_i
       ; print, old_month
        FOR i =0,  N_ELEMENTS(X_label2)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i2, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label2)-1 THEN $
                        X_label2[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR        
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; plot the Dst time series for each event
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
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

idate0 = string(yr_i, mh_i, dy_i, format='(I4,I02,I02)')
TGM_n = idate0
case TGM_n of
    '20031119' : TGM_n = 1
    '20041106' : TGM_n = 2
    '20050514' : TGM_n = 3
    '20150315' : TGM_n = 4
    '20170526' : TGM_n = 5
    '20170906' : TGM_n = 6
    else: print, 'fuera de rango'
endcase

window_title = 'TGM'+ string(TGM_n, format='(I01)')+', '+ $
string(old_month, yr_i, format='(A, X, I4)')

time_title = 'Time [UT]' 
                
    up  = max(dst-H)
    down= min(dst-H)
    up2      = max(tec-med)
    down2    = min(tec-med)
    
    if up2 gt up then up0 = up2 else up0 = up
    if down2 lt down then down0 = down2 else down0 = down
       dH = TeXtoIDL('\DeltaH')   
                
    plot, days_dst, dst-H, XTICKS=tw, xminor = 8, POSITION=[0.1,0.18,0.9,0.85],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down,up],$
    title = '', BACKGROUND = blanco, COLOR=negro, XRANGE=[0, tw],$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', tw+1), CHARSIZE = 0.9,$
	thick=2
	

    ;oplot, days_dst, H, COLOR=rojo
  
 ;   if up-down gt 300 then begin
 ;       for i = -600., 100., 100. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
 ;   endif else begin
 ;       for i = -600., 100., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
  ;  endelse

;#####################################################################################
;TGM1
    ;oplot, [1.3,1.3], [down-100,up+100], linestyle=2, COLOR=negro 
    ;oplot, [2.5,2.5], [down-100,up+100], linestyle=2, COLOR=negro  
;##################################################################################### 
;TGM2 
    ;oplot, [1.5,1.5], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2.55,2.55], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
    
    ;oplot, [4.6,4.6], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [3.55,3.55], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08       
;##################################################################################### 
;TGM3 
    ;oplot, [0.9,0.9], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2,2], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
;#####################################################################################    
;TGM4 
    ;oplot, [3.1,3.1], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
   ; oplot, [2,2], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
;#####################################################################################  
;TGM5 
    ;oplot, [3,3], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [1.7,1.7], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
;#####################################################################################       
;TGM6 
   ; oplot, [3,3], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
   ; oplot, [1.7,1.7], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
;#####################################################################################
    tec_diff = tec-med
    
    med_tec = MEDIAN(tec_diff)
    std_tec = stddev(tec_diff)
    
    sup = med_tec+std_tec
    inf = med_tec-std_tec
    
    l_sup = fltarr(n_elements(tec_diff))
    l_sup[*] = sup

    l_inf = fltarr(n_elements(tec_diff))
    l_inf[*] = inf    
    i = where((tec_diff ge sup or tec_diff le inf), complement=ni)
    tec_diff1 = tec_diff[i]
    tec_diff0 = tec_diff[ni]
    days1      = tec_days[i]
    days0     = tec_days[ni]    
    tec_diff[ni] = !Values.F_NAN

        plot, tec_days, tec_diff, color=rojo, XTICKS=tw, xminor = 8, POSITION=[0.1,0.18,0.9,0.85],$
        XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down2,up2],$
        BACKGROUND = blanco, XRANGE=[0, tw], /noerase,$
        XTICKNAME=REPLICATE(' ', tw+1), CHARSIZE = 0.9, thick=3
        ;oplot, tec_days, tec_diff, color=rojo, psym=2
        
    tec_diff[ni] = tec_diff0
    tec_diff[i]  = !Values.F_NAN      
        oplot, tec_days, tec_diff, color=rojo, linestyle=2
     
    oplot, tec_days, l_sup, color=negro, linestyle=3
    oplot, tec_days, l_inf, color=negro, linestyle=3
 ;   oplot, tec_days, med, COLOR=rojo
    
 ;   if up2-down2 gt 120 then begin
 ;       for i = -100., 260., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
 ;   endif else begin
 ;       for i = -100., 260., 10. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
 ;   endelse 
;#####################################################################################
;TGM1
    ;oplot, [1.3,1.3], [down-100,up+100], linestyle=2, COLOR=negro 
    ;oplot, [2.5,2.5], [down-100,up+100], linestyle=2, COLOR=negro  
;#####################################################################################    
;TGM2
    ;oplot, [1.5,1.5], [down2-200,up2+200], linestyle=2, COLOR=negro 
    ;oplot, [2.55,2.55], [down2-200,up2+200], linestyle=2, COLOR=negro
    
    ;oplot, [4.6,4.6], [down2-100,up2+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [3.55,3.55], [down2-100,up2+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08         
;#####################################################################################  
;TGM3 
    ;oplot, [0.9,0.9], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2,2], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08       
;#####################################################################################   
;TGM4 
    ;oplot, [3.1,3.1], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
   ; oplot, [2,2], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08       
;#####################################################################################         
;TGM5 
   ; oplot, [3,3], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
   ; oplot, [1.7,1.7], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08       
;##################################################################################### 
;TGM6 
   ; oplot, [3,3], [down-100,up+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
   ; oplot, [1.7,1.7], [down-100,up+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08       
;#####################################################################################
        AXIS, XAXIS = 0, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         ;XTICKFORMAT='(A1)',$
                         xtitle='tiempo [UT]', $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,tw], $
                         XTICKS=tw, $
                         XMINOR=8, $
                         ;XTICKFORMAT='(A1)',$                         
                         XTICKNAME=X_label2, $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0], $
                         YTITLE = 'Dst-'+dH+' [nT]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $
                        ; CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00
                        
        AXIS, YAXIS = 1, YRANGE=[down2,up2], $
                         YTITLE = 'TEC-<TEC> [TECu]', $   
                         COLOR=rojo, $
                         ystyle=2, $
                         CHARSIZE = 0.9;, $
                        ; CHARTHICK=chr_thick1;, $    
             
;   a_diff = a_mex-ap
;up3     = max(a_diff)
;down3   = min(a_diff)

;    plot, k_days, a_diff, XTICKS=tw, xminor=8, /noerase, $
 ;    XSTYLE = 5, ystyle = 6, XRANGE=[0, tw], yrange=[down3,up3],$
  ;  BACKGROUND = blanco, COLOR=negro, CHARSIZE = chr_size1, $
   ; CHARTHICK=chr_thick1, POSITION=[0.1,0.35,0.93,0.65]
    ;                XTICKNAME=REPLICATE(' ', tw+1)
   ; oplot, k_days, a_mex, color=verde
;    if up3-down3 gt 300 then begin
 ;       for i = -500., 500., 100. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
 ;   endif else begin
  ;      for i = -300., 300., 50. do oplot, [0,tw], [i,i], linestyle=1, COLOR=negro
   ; endelse 
;#####################################################################################
;TGM1
    ;oplot, [1.3,1.3], [down-100,up+100], linestyle=2, COLOR=negro 
    ;oplot, [2.5,2.5], [down-100,up+100], linestyle=2, COLOR=negro 
;#####################################################################################  
;TGM2
   ; oplot, [1.5,1.5], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2.55,2.55], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08
    
    ;oplot, [4.6,4.6], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [3.55,3.55], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08    
;#####################################################################################      
;TGM3 
    ;oplot, [0.9,0.9], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2,2], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08   
;##################################################################################### 
;TGM4 
   ; oplot, [3.1,3.1], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [2,2], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08   
;#####################################################################################   
;TGM5 
    ;oplot, [3,3], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
    ;oplot, [1.7,1.7], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08   
;#####################################################################################      
;TGM6 
 ;   oplot, [3,3], [down3-100,up3+100], linestyle=2, COLOR=negro ;18:00 hr del nov 07
 ;   oplot, [1.7,1.7], [down3-100,up3+100], linestyle=2, COLOR=negro ;  07:12 hr del nov 08   
;#####################################################################################  

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
                write_jpeg, path+'idx_V4'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'idx_V4'+Date+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'mag_idx_'+Date+'_V4.png'
                print, ''
        ENDIF
        RETURN
end    
;###############################################################################
;############################################################################### 
;###############################################################################
;###############################################################################

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
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
                      DOY : 0, Kp: 0, Kp_str: '', Ap: 0}

		resulting_data0 = REPLICATE(DataStruct, number_of_lines-header)	
        ;print, number_of_lines-header-1, number_of_lines-header+1
                
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
		
		name = 'teo_'+date+'.index.'
		
		file_name = '../rutidl/Kmex/'+name+'final'		
		
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

		READS, data, tmp_var, FORMAT='(I3, I4, I4, I4, I4, I4, I4, I4, I4)'
		
		idx_kmex.k_mex[*]   = tmp_var[0].x
		idx_kmex.a_mex[*]   = tmp_var[1].x
        idx_kmex.k_mex_sum  = tmp_var[0].y
        idx_kmex.a_med      = tmp_var[1].y				
		return, idx_kmex	
end    

pro list_kp, r_kp, in
    
    On_error, 2
	compile_opt idl2, HIDDEN	

	a = in
		
	dat = kp_data(a)    
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Defining variables to write output files indicating the events where Dst index
; got lower than -150 nT.
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
	tmp = n_elements(dat.year)
	indexes = WHERE(dat.Kp GE 70, count)
	IF count LE 0 THEN MESSAGE, 'ERROR'
		
	kp0      = dat.Kp
    doy0     = dat.DOY
    year0    = dat.year
    month0   = dat.month
    day0     = dat.day
    hour0    = dat.hour 
    
    idx_kp  = kp0[indexes]
	doy     = doy0[indexes]
	year    = year0[indexes]
	month   = month0[indexes]
	day     = day0[indexes]
	hour    = hour0[indexes]

	idx_kp_tmp   = idx_kp*0
	doy_tmp      = doy*0
	year_tmp     = year*0
	month_tmp    = month*0
	day_tmp      = day*0
	hour_tmp     = hour*0

    Date    = strmid(string(dat.year[0]), 8, 4)
    outfile = '../rutidl/output/'+'kp_'+Date[0]+'TGM.txt'; defining the 
    ;output file names and path directories
    
    j=0l
    idx_kp_tmp[j] = idx_kp_tmp[0]
	year_tmp[j]   = year[0]
	month_tmp[j]  = month[0]
    day_tmp[j]    = day[0]
	hour_tmp[j]   = hour[0]
    doy_tmp[j]    = doy[0]  

    FOR i=1l, N_ELEMENTS(idx_kp)-1 DO BEGIN
	IF doy_tmp[j] EQ doy[i] THEN BEGIN
	        IF idx_kp[i] GT idx_kp_tmp[j] THEN BEGIN
        	        idx_kp_tmp[j] = idx_kp[i]
		            year_tmp[j]   = year[i]
	                month_tmp[j]  = month[i]
        	        day_tmp[j]    = day[i]
	                hour_tmp[j]   = hour[i]
	                doy_tmp[j]    = doy[i]
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
	        idx_kp_tmp[j] = idx_kp[i]
	        year_tmp[j]   = year[i]
	        month_tmp[j]  = month[i]
        	day_tmp[j]    = day[i]
	        hour_tmp[j]   = hour[i]
	        doy_tmp[j]    = doy[i]
	ENDELSE
    ENDFOR
   
    ;OPENW, lun, outfile, /GET_LUN

    ;PRINTF,lun, tmp[i], i, format = "(A50, 2X, I5)"
    ;idx   = where(idx_kp GT 6, count) 

    ;idx = findgen(n_elements(idx_kp))
 
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice Kp ascendió por'
    print, 'encima de 7, significando un evento de tormenta geomagnética intensa.'
    print, '                                                                    '
    
    print, format='(2X, "Fecha", 4X, "Hora", 3X, "DOY", 2X, "Indice Kp")'
    print, '---------------------------------'
       
    for i=0, N_ELEMENTS(indexes)-1 do begin
        ;idx   = where(idx_kp GT 6, count)
        
            if doy_tmp[i] ne 0 then begin

            ;idx_kp[idx] = idx_kp[i]

                print, year_tmp[i], month_tmp[i], day_tmp[i], hour_tmp[i], doy_tmp[i], idx_kp_tmp[i], $
                             FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I03, 5X, I02)'  
                ;printf, lun, year_tmp[i], month_tmp[i], day_tmp[i], hour_tmp[i], doy_tmp[i], idx_kp_tmp[i], $
                           ; FORMAT = '(I4, X, I02, X, I02, 2X, I02, 4X, I03, X, I02)'                      
        endif        
    endfor    
   ; close,lun
    ;FREE_LUN, lun
    print, '                                                                   '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento'

end     

pro list_kmex, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]	
;##############################################################################
; reading data files
;##############################################################################
        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        data_file_name = strarr(file_number)

        string_date     = strarr(file_number)
       
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        file = FILE_SEARCH(data_file_name[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF		        
	        
        ENDFOR
	
        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        k_mex    = INTARR(file_number*8)
        a_mex    = INTARR(file_number*8)

        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        dat = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = dat.k_mex[*]
                        a_mex[i*8:(i+1)*8-1] = dat.a_mex[*]
                ENDIF 
        ENDFOR
    
    k_mex[where(k_mex[*] eq 99)] = !Values.F_NAN     

	idx = WHERE(k_mex GE 7, count)
	IF count LE 0 THEN MESSAGE, 'ERROR'

    time_w = timegen(n_elements(file_number), final=julday(mh_f, dy_f, yr_f, 23), $
                start=julday(mh_i, dy_i, yr_i, 0) , units='H')

    caldat, time_w, m, d, y, hr

    y_idx        = y[idx]   
    m_idx        = m[idx]
    d_idx        = d[idx]
    k_mex_idx    = k_mex[idx]
    hr_idx       = hr[idx]
    
    y_idx_tmp        = y_idx*0
    m_idx_tmp        = m_idx*0
    d_idx_tmp        = d_idx*0
    k_mex_idx_tmp    = k_mex_idx*0
    hr_idx_tmp       = hr_idx*0
;print, H[idx]
    j=0l
    y_idx_tmp[j]        = y_idx[0]
    m_idx_tmp[j]        = m_idx[0]
    d_idx_tmp[j]        = d_idx[0]
    k_mex_idx_tmp[j]    = k_mex_idx[0]
    hr_idx_tmp[j]       = hr_idx[0]
    
    FOR i=1l, N_ELEMENTS(k_mex_idx)-1 DO BEGIN
	IF d_idx_tmp[j] EQ d_idx[i] THEN BEGIN
	        IF k_mex_idx[i] GT k_mex_idx_tmp[j] THEN BEGIN
        	        k_mex_idx_tmp[j]   = k_mex_idx[i]
		            y_idx_tmp[j]       = y_idx[i]
	                m_idx_tmp[j]       = m_idx[i]
        	        d_idx_tmp[j]       = d_idx[i]
	                hr_idx_tmp[j]      = hr_idx[i]	               
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
        	    k_mex_idx_tmp[j]   = k_mex_idx[i]
		        y_idx_tmp[j]       = y_idx[i]
	            m_idx_tmp[j]       = m_idx[i]
        	    d_idx_tmp[j]       = d_idx[i]
	            hr_idx_tmp[j]      = hr_idx[i]	
	ENDELSE
    ENDFOR
    
    ;outfile='../rutidl/output/'+'Kmex_list_TGM.txt'
     ;   OPENW, lun, outfile, /GET_LUN, /append
        
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice DH ascendió por'
    print, 'encima de 6, significando un evento de tormenta geomagnética '
    print, 'intensa.'
    print, '                                                                    '    
    print, format='(2X, "Fecha", 4X, "Hora", 2X, "Indice Kmex")'
    print, '---------------------------------'
    print, '                                                                    '
         
    for i=0, N_ELEMENTS(idx)-1 do begin
           if d_idx_tmp[i] ne 0 then begin

                print, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], hr_idx_tmp[i], $
                k_mex_idx_tmp[i], $
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I03)' 
  
               ; printf, lun, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], $
               ; hr_idx_tmp[i], k_mex_idx_tmp[i], $ 
               ; FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I03)'                         
        endif        
    endfor
    ;close,lun
    ;FREE_LUN, lun
   ;PRINT, DATETIME[i], Dst[i], format = "(A29)"
   print, '                                                                    '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento' 
end
