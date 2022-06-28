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
FUNCTION dst_data, date

	On_error, 2
	COMPILE_OPT idl2, HIDDEN
    
	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
		file_name = path+'/dst/daily/dst_'+date+'.txt'
		;print, date
        header = 1             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'		
;	    IF opened_files NE N_ELEMENTS(file) THEN BEGIN
;	        file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_P.dat'
;	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    ;IF opened_files NE N_ELEMENTS(file) THEN BEGIN
    	   ;     file_name = '../rutidl/dst/Dst_'+ year+'-01-01_'+year+'-12-31_Q.dat'
	      ;      file = FILE_SEARCH(file_name, COUNT=opened_files)    	        
    	 ;       IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
    	;    ENDIF    	    	    
	  ;  ENDIF

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
        second : 0, DOY : 0, Dst: 0}

		r_dst = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], r_dst, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,X,I2,I04,I5)'
		RETURN, r_dst
END


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

        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'		
		file_name = path+'/dH_teo/'+'teo_'+date+'.dst.early'
		
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


FUNCTION tec_data, idate
	On_error, 2
	compile_opt idl2, HIDDEN

	iyear	= idate[0]
	imonth	= idate[1]
	iday 	= idate[2]		

        header = 1      ; Defining number of lines of the header 
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        idate = string(iyear, imonth, iday, format = '(I4, "-", I02, "-", I02)')
		file_name = path+'/tec/'+'tec_'+idate+'.txt'
		
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
		RETURN, r_tec
END

FUNCTION kp_data, date

	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		file_name = path+'/kp/daily/kp_'+date+'.txt'
		;print, date
        header = 1             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-		
		file = FILE_SEARCH(file_name, COUNT=opened_files)

	  ;  IF opened_files NE N_ELEMENTS(file) THEN begin
	  ;      file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_P.dat'
	  ;      file = FILE_SEARCH(file_name, COUNT=opened_files) 
    ;	    IF opened_files NE N_ELEMENTS(file) THEN begin
    ;	        file_name = '../rutidl/kp/Kp_'+ year+'-01-01_'+year+'-12-31_Q.dat'
	 ;           file = FILE_SEARCH(file_name, COUNT=opened_files)    	        
    	;        IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
    	 ;   ENDIF    	    	    
;	    ENDIF     
        
		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, second : 0,$
                      DOY : 0, Kp: 0, Kp_str: '', Ap: 0}

		resulting_data0 = REPLICATE(DataStruct, number_of_lines-header)	
        ;print, number_of_lines-header-1, number_of_lines-header+1
                
		READS, data[header:number_of_lines-1], resulting_data0, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,X,I2,I04,X,I1,A1,I04)'
		
	;	indexes_07 = WHERE(resulting_data0[*].Kp_str EQ '-')
	;	indexes_03 = WHERE(resulting_data0[*].Kp_str EQ '+')
	;	indexes_00 = WHERE(resulting_data0[*].Kp_str EQ 'o')
		
	;	Kp_tmp = resulting_data0[*].Kp*10
		
	;	Kp_tmp[indexes_07] = Kp_tmp[indexes_07]-3
	;	Kp_tmp[indexes_03] = Kp_tmp[indexes_03]+3

        DataStruct2 = {year : 0, month : 0, day : 0, hour : 0, minute: 0, second : 0,$
                       DOY : 0, Kp: 0, Ap: 0}

		r_kp = REPLICATE(DataStruct2, number_of_lines-header)	

		r_kp[*].year   = resulting_data0[*].year
		r_kp[*].month  = resulting_data0[*].month
		r_kp[*].day    = resulting_data0[*].day
		r_kp[*].hour   = resulting_data0[*].hour
		r_kp[*].minute = resulting_data0[*].minute
		r_kp[*].DOY    = resulting_data0[*].DOY
		;r_kp[*].Kp     = Kp_tmp[*]
		r_kp[*].Ap     = resulting_data0[*].AP				
		return, resulting_data0;r_kp
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
		path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		file_name = path+'/Kmex/'+name+'final'		
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name =  path+'/Kmex/'+name+'early'
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

;###############################################################################
;###############################################################################
;###############################################################################
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

    Date    = string(year[0], mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
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
                data_file_name[i] = '../master_thesis/datos/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        file = FILE_SEARCH(data_file_name[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name[i] = '../master_thesis/datos/Kmex/'+'teo_'+string_date[i]+'.index.early'    
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
                data_file_name[i] = '../master_thesis/datos/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
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
