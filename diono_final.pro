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
;###############################################################################
FUNCTION dst_data, date

	On_error, 2
	COMPILE_OPT idl2, HIDDEN
    
	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
		file_name = path+'/dst/daily/dst_'+date+'.txt'
		;print, date
        header = 1             ; Defining number of lines of the header 
;###############################################################################
;reading data files
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'		

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

FUNCTION DH_teo, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
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
;###############################################################################
;extracting data and denfining an structure data
        DStruct = {hora : 0, D_stdesv : 0., D : 0., H_stdesv : 0., H : 0., $
        Z_stdesv : 0., Z : 0., N_stdesv : 0., N : 0., F_stdesv : 0., F : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(I2, F10, F8, F10, F10, F10, F10, F10, F10, F10, F10)'		
		RETURN, teo_mag		
END


FUNCTION teo, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        date = string(year, month, day, format = '(I4, I02, I02)')
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'		
		file_name = path+'/teoloyucan/hourly/teo_'+date+'h.dat'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;###############################################################################
;extracting data and denfining an structure data
        DStruct = {H : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(F10)'		
		RETURN, teo_mag		
END


FUNCTION kmex, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        date = string(year, month, day, format = '(I4, I02, I02)')
		
		name = 'teo_'+date+'.index.'
		path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		file_name = path+'/Kmex/'+name+'final'		
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name =  path+'/Kmex/'+name+'early'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
	    
	    ENDIF

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;###############################################################################
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
		RETURN, idx_kmex	
END    

FUNCTION new_kmex, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        date = string(year, month, day, format = '(I4, I02, I02)')
		path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		
		name = 'new_km'+date+'.txt'
		
		file_name = path+'/Kmex/new_kmex/'+name		
		file = FILE_SEARCH(file_name, COUNT=opened_files)	

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;###############################################################################
;extracting data and denfining an structure data
        idx_newkmex      = {n_kmex1      : intarr(8), $
                            n_kmex2      : intarr(8)}
        
        struct = {x : [0, 0, 0, 0, 0, 0, 0, 0]}
        
        tmp_var = replicate(struct, 2)  
		READS, data, tmp_var, FORMAT='(I3, I3, I3, I3, I3, I3, I3, I3)'
		
		idx_newkmex.n_kmex1[*]   = tmp_var[0].x
		idx_newkmex.n_kmex2[*]   = tmp_var[1].x		
		
		RETURN, idx_newkmex	
END

FUNCTION kp_data, date

	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		file_name = path+'/kp/daily/kp_'+date+'.txt'
		;print, date
        header = 1             ; Defining number of lines of the header 
;###############################################################################
;reading data files
		file = FILE_SEARCH(file_name, COUNT=opened_files)
        
		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;###############################################################################
;extracting data and denfining an structure data
        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, second : 0,$
                      DOY : 0, Kp: 0, Kp_str: '', Ap: 0}


		r_kp = REPLICATE(DataStruct, number_of_lines-header)	
                
		READS, data[header:number_of_lines-1], r_kp, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,X,I2,I04,I02,A1,I04)'
			
		RETURN, r_kp
END   


FUNCTION baseline_sq, date
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
		
        DStruct = {Bsq : 0.}
		B_sq = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], B_sq, $
	 format='(F08.4)' 		
		RETURN, B_sq
END

	
pro diono_final, date_i, date_f, JPEG = jpeg 

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
    file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 	

    tot_days= findgen(file_number*24)/24.0  
    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')    
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
; define DH variables
;###############################################################################

        
        data_file_name_h  = strarr(file_number)                                   
        data_file_name_dh  = strarr(file_number)
        data_file_name_km  = strarr(file_number)
        data_file_name_kmn = strarr(file_number) 
        data_file_name_bsq  = strarr(file_number) 
        data_file_name_kp  = strarr(file_number)  
        data_file_name_dst = strarr(file_number)                
                                              
        string_date        = strarr(file_number)
        string_date_2        = strarr(file_number)        
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date_2[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')   
                             
                data_file_name_h[i]  = data_path+'/teoloyucan/hourly/'+'teo_'+string_date[i]+'h.dat'                
                data_file_name_dh[i] = data_path+'/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
                data_file_name_dst[i]= data_path+'/dst/daily/dst_'+string_date_2[i]+'.txt'
                data_file_name_km[i] = data_path+'/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        data_file_name_kmn[i] = data_path+'/Kmex/new_kmex/new_km'+string_date[i]+'.txt'
		        data_file_name_kp[i]= data_path+'/kp/daily/kp_'+string_date_2[i]+'.txt' 		        
                data_file_name_bsq[i] = '../rutidl/output/Bsq_baselines/Bsq_'+string_date[i]+'h.dat'
 
		                        		        
		        file = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)
		        file_bsq = FILE_SEARCH(data_file_name_bsq[i], COUNT=opened_files)		        
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name_km[i] = data_path+'/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF 	                            
        ENDFOR

        exist_data_file_h   = FILE_TEST(data_file_name_h)
        capable_to_plot_h   = N_ELEMENTS(where(exist_data_file_h EQ 1))

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        exist_data_file_kmn   = FILE_TEST(data_file_name_kmn)
        capable_to_plot_kmn   = N_ELEMENTS(where(exist_data_file_kmn EQ 1))
 
        exist_data_file_bsq   = FILE_TEST(data_file_name_bsq)
        capable_to_plot_bsq   = N_ELEMENTS(where(exist_data_file_bsq EQ 1))   
        
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.km_index.',A,' impossible to plot all data.')"              
        ENDIF
;###############################################################################
; Generate the time variables to plot dH time series                          
        dH    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_dh)-1 DO BEGIN
                IF exist_data_file_dh[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_dh = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        dH[i*24:(i+1)*24-1] = d_dh.H[*]
                                                                                                                       
                ENDIF ELSE BEGIN
                        dH[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
;###############################################################################
; Generate the time variables to plot H time series  
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
;###############################################################################
; Generate the time variables to plot kmex time series          
        k_mex    = fltarr(file_number*8)
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
;###############################################################################
; Generate the time variables to plot new kmex time series  
        new_kmex1    = fltarr(file_number*8)
        new_kmex2    = fltarr(file_number*8)
        FOR i = 0, N_ELEMENTS(exist_data_file_km)-1 DO BEGIN
                IF exist_data_file_kmn[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'                 
                        d_km = new_kmex([tmp_year, tmp_month, tmp_day])
                        
                        new_kmex1[i*8:(i+1)*8-1] = d_km.n_kmex1[*]/10.0
                        new_kmex2[i*8:(i+1)*8-1] = d_km.n_kmex2[*]/10.0
                ENDIF             
        ENDFOR
        
    k_days = findgen(file_number*8)/8.0     
;###############################################################################
;Dst Data                       
        exist_data_file_dst   = FILE_TEST(data_file_name_dst)
        capable_to_plot_dst   = N_ELEMENTS(where(exist_data_file_dst EQ 1))

        IF capable_to_plot_dst NE N_ELEMENTS(data_file_name_dst) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dst_index.',A,' impossible to plot all data.')"              
        ENDIF

        dst    = FLTARR(file_number*24)                               
        FOR i = 0, N_ELEMENTS(exist_data_file_dst)-1 DO BEGIN
                IF exist_data_file_dst[i] EQ 1 THEN BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date_2[i] = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')                
                        dat = dst_data([tmp_year, tmp_month, tmp_day])
                        
                        dst[i*24:(i+1)*24-1] = dat.Dst[*]                                                
                                                                                              
                ENDIF ELSE BEGIN
                         dst[i*24:(i+1)*24-1] = 999999.0                      
                ENDELSE                
        ENDFOR       
;###############################################################################
;Kp Data                       
        exist_data_file_kp   = FILE_TEST(data_file_name_kp)
        capable_to_plot_kp   = N_ELEMENTS(where(exist_data_file_kp EQ 1))

        IF capable_to_plot_kp NE N_ELEMENTS(data_file_name_kp) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.kp_index.',A,' impossible to plot all data.')"              
        ENDIF

        kp    = FLTARR(file_number*8)
        ap    = FLTARR(file_number*8)
        str   = STRARR(file_number*8)                              
        FOR i = 0, N_ELEMENTS(exist_data_file_kp)-1 DO BEGIN
                IF exist_data_file_kp[i] EQ 1 THEN BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date_2[i] = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')                
                        dat = kp_data([tmp_year, tmp_month, tmp_day])
                        kp[i*8:(i+1)*8-1]   = dat.Kp[*]                                               
                        ap[i*8:(i+1)*8-1]   = dat.Ap[*]
                        str[i*8:(i+1)*8-1]= dat.Kp_str[*]                                                                      
                ENDIF ELSE BEGIN
                         kp[i*8:(i+1)*8-1] = 999999.0
                         ap[i*8:(i+1)*8-1] = 999999.0                     
                ENDELSE                
        ENDFOR
        
		indexes_07 = WHERE(str EQ '-')
		indexes_03 = WHERE(str EQ '+')
		indexes_00 = WHERE(str EQ 'o')
		print, indexes_07, indexes_03
    IF indexes_03[0] NE -1 AND indexes_07[0] NE -1 THEN BEGIN
        kp[indexes_07] = kp[indexes_07]-0.3
        kp[indexes_03] = kp[indexes_03]+0.3
    ENDIF                  
;###############################################################################        
;identifying NAN values in the Time Series
        i_nan1 = where(dH eq 999999.0, ncount)
        i_nan2 = where(dH gt 100.0, n2count)
        
        prcent_nan = FLOAT(ncount+n2count)*100.0
        print,'porcentaje de valores NaN:', prcent_nan/n_elements(H),'%'
        
        for i=0, n_elements(dH)-1 do begin
            if dH[i] eq 999999.0 then begin
                dH[where(dH[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(dH)-1 do begin
            if dH[i] ge 100.0 then begin
                dH[where(dH[*] ge 100.0)] = !Values.F_NAN          
            endif
        endfor                
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo        
    dH_tmp   = dH
    dH_exist = where(finite(dH_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
  
    ; interpolate at the locations of the bad data using the good data    
    if nbaddata gt 0 then dH_tmp[baddata] = interpol(dH_tmp[dH_exist], dH_exist, baddata)
    dH = dH_tmp 
;############################################################################### 
        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(H)-1 do begin
            if H[i] EQ 99999.0 then begin
                H[where(H[*] EQ 99999.0)] = !Values.F_NAN          
            endif
        endfor  
        
        for i=0, n_elements(k_mex)-1 do begin
            if k_mex[i] GT 9.0 then begin
                k_mex[where(k_mex[*] GT 9.0)] = !Values.F_NAN          
            endif
        endfor                
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo
       
    H_tmp   = H
    H_exist = where(finite(H_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
    ; interpolate at the locations of the bad data using the good data
    
    if nbaddata gt 0 then H_tmp[baddata] = interpol(H_tmp[H_exist], H_exist, baddata)
    H = H_tmp    
    H_tmp   = H
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
; define Diono  
   mlat         = 28.06*!pi
   ld           = cos(mlat/180)
   p_a          = dst*ld
   baseline     = Bsq + p_a             
   diono   = H-baseline
   time = 3600.0
  ; print, H
   fn      = float(1.0/(2.0*time)) ; frecuencia de Nyquist      
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
; define filtering 
coeff_ddyn    = digital_filter(passband_l/fn, passband_u/fn, 50, 18)
coeff_dp2   = digital_filter(highpass_l/fn, 1.0, 50, 4)
;###############################################################################
; define disturbing DIONO effects
ddyn        = convol(diono, coeff_ddyn, /edge_wrap)
dp2         = convol(diono, coeff_dp2, /edge_wrap)        
;###############################################################################
; define device and color parameters 
        Device_bak2 = !D.Name         
        SET_PLOT, 'Z'      
        
        Xsize=fix(1600)
        Ysize=1000
        DEVICE, SET_RESOLUTION = [Xsize,Ysize],Set_Pixel_Depth=24, DECOMPOSED=1  
        DEVICE, z_buffer=4
        DEVICE, set_character_size = [10, 12] 
        
        chr_size1 = 0.9
        chr_thick1= 1.5
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
     X_label   = STRARR(file_number+1)+' '
        months    = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic']
        old_month = mh_i
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, $
                        FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
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
       days = intarr(file_number+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       days = days*24/24. 
       day_time = findgen(24)   
;############################################################################### 
    time_title = ' Tiempo Universal [dias]'
    window_title = 'TGM'+ STRING(TGM_n, FORMAT='(I01)')+', '+ $
                    STRING(old_month, yr_i, FORMAT='(A, X, I4)')                          
;###############################################################################                                    
     up = max(dH)
     down=min(dH)
    ; print, up, down
     PLOT, tot_days, dH, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = 0.9, CHARTHICK=chr_thick1, $
     POSITION=[0.1,0.55,0.95,0.9], XSTYLE = 5, XRANGE=[0, file_number], YSTYLE = 6,$
     XTICKNAME=REPLICATE(' ', file_number+1), YRANGE=[down,up], THICK=4
;###############################################################################        
     d_H = TeXtoIDL('\DeltaH') 
     dst_l = TexToIDL('Dst_\lambda')
     
     Bsq_med   = MEAN(Bsq)
     Bsq0       = Bsq-Bsq_med
     new_index  = p_a+ddyn+dp2
     OPLOT, tot_days, dst, COLOR=verde, THICK=4  
     OPLOT, tot_days, new_index, COLOR=rojo, THICK=4         
     
     tf          = N_ELEMENTS(tot_days)
     dist_func1  = ABS(dH-new_index)
     tot_diff1   = TOTAL(dist_func1)
     avr_diff1   = tot_diff1 / (tf)

     dist_func0  = ABS(dH-dst)
     tot_diff0   = TOTAL(dist_func0)
     avr_diff0   = tot_diff0 / (tf)
     
     dist_func2  = ABS(dH[34:72]-new_index[34:72])
     tot_diff2   = TOTAL(dist_func2)
     avr_diff2   = tot_diff2 / (N_ELEMENTS(dH[34:72]))   

     dist_func3  = ABS(dH[34:72]-dst[34:72])
     tot_diff3   = TOTAL(dist_func3)
     avr_diff3   = tot_diff3 / (N_ELEMENTS(dH[34:72]))          
     print, '###################################################################'  
     print, 'funcion distancia d(DH,Dst): ', tot_diff0
     print, '<Distancia total>', avr_diff0  
     print, '###################################################################'          
     print, 'funcion distancia total d(DH,Dst+PI): ', tot_diff1
     print, '<Distancia total>: ', avr_diff1      
     print, '###################################################################'
     print, 'funcion distancia d(DH,Dst) en la ventana de t: ', tot_diff3
     print, '<Distancia total>', avr_diff3     
     print, '###################################################################'          
     print, 'funcion distancia d(DH,Dst+PI) en la ventana de t: ', tot_diff2
     print, '<Distancia total>', avr_diff2     
     print, '###################################################################'               
;###############################################################################                         
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 1.2 , $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $                                            
                         XMINOR=8, $ 
                         CHARSIZE = 1.2 , $                       
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5

        AXIS, YAXIS = 0, YRANGE=[down,up], $
                         YTITLE = '', $                          
                         COLOR=negro, $
                         YSTYLE=2, $
                         CHARSIZE = 1.2,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $
                         COLOR=negro, $                                                                      
                         YSTYLE=2, $
                         CHARSIZE = 1.2,$
                         CHARTHICK=1.5                                                                                        
;###############################################################################    
     PLOT, k_days, k_mex, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = 0.9, CHARTHICK=chr_thick1, $
     POSITION=[0.1,0.1,0.95,0.45], XSTYLE = 5, XRANGE=[0, file_number], ySTYLE = 5,$
     XTICKNAME=REPLICATE(' ', file_number+1), YRANGE=[0,9], THICK=4, /NOERASE
     
     OPLOT, k_days, kp, COLOR=verde, LINESTYLE=0, THICK=4
     A = FINDGEN(17) * (!PI*2/16.)
     USERSYM, COS(A), SIN(A), /FILL     
     OPLOT, k_days, kp, COLOR=negro, PSYM=8, THICK=4, symsize=1  
FOR i=0, N_ELEMENTS(Kp)-1 DO BEGIN         
     IF Kp[i] GT k_mex[i] THEN BEGIN
        ERRPLOT,k_days[i], new_kmex1[i], Kp[i], COLOR=negro, THICK=3
        ERRPLOT, k_days[i], Kp[i], new_kmex2[i], COLOR=negro, LINESTYLE=0
    ENDIF  
        IF Kp[i] LT k_mex[i] THEN BEGIN
        ERRPLOT,k_days[i], new_kmex1[i], Kp[i], COLOR=negro, LINESTYLE=0
        ERRPLOT, k_days[i], Kp[i], new_kmex2[i], COLOR=negro, THICK=3
    ENDIF
ENDFOR                  
;###############################################################################                         
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 1.2 , $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                         
                         IF TGM_n EQ 3 THEN Ltime = 0.25
                         IF TGM_n EQ 5 THEN Ltime = 0.25
                         IF TGM_n EQ 6 THEN Ltime = 0.25
                         IF TGM_n EQ 1 THEN Ltime = 5./24.
                         IF TGM_n EQ 2 THEN Ltime = 5./24.
                         IF TGM_n EQ 4 THEN Ltime = 5./24.                                                  
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-Ltime), $
                         XTICKS=file_number, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $                                            
                         XMINOR=8, $ 
                         CHARSIZE = 1.2 , $                       
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5

        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTITLE = '', $                          
                         COLOR=negro, $
                         YSTYLE=1, $
                         CHARSIZE = 1.2,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         COLOR=negro, $                     
                         YSTYLE=1, $
                         CHARSIZE = 1.2,$
                         CHARTHICK=1.5                                                 
;###############################################################################
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.95   
   XYOUTS, X, y, window_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.85, CHARTHICK=1.5    

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.928, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.3, CHARTHICK=1.5    

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.475, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.3, CHARTHICK=1.5   
;############################################################################### 
   XYOuts, 0.14, 0.57, '(a)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.8, CHARTHICK= 4    

   
   XYOuts, 0.14, 0.41, '(b)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.8, CHARTHICK= 4          
;###############################################################################                     
;first panel legend                   
        POLYFILL, [0.79,0.82,0.82,0.79], [0.654,0.654,0.657,0.657], color = negro, /NORMAL
        POLYFILL, [0.79,0.82,0.82,0.79], [0.624,0.624,0.627,0.627], color = rojo, /NORMAL        
        POLYFILL, [0.79,0.82,0.82,0.79], [0.594,0.594,0.597,0.597], color = verde, /NORMAL   
        
        XYOUTS, 0.825, 0.65 , /NORMAL, $
                d_H, COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1 
                
        XYOUTS, 0.825, 0.62 , /NORMAL, $
               dst_l, COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1   
                
        XYOUTS, 0.825, 0.59 , /NORMAL, $
                'Dst', COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1     
                
;second panel legend  
        POLYFILL, [0.79,0.82,0.82,0.79], [0.424,0.424,0.427,0.427], color = negro, /NORMAL
        POLYFILL, [0.79,0.82,0.82,0.79], [0.394,0.394,0.397,0.397], color = verde, /NORMAL  

                XYOUTS, 0.825, 0.42 , /NORMAL, $
                'Kmex', COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1   
                
        XYOUTS, 0.825, 0.39 , /NORMAL, $
                'Kp', COLOR=negro, $
                CHARSIZE = 1.4, $
                CHARTHICK=chr_thick1   
;###############################################################################                
   y = (0.45 - 0.1) / 2. + 0.1 
   XYOUTS, 0.05, y, '[nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5   

   y = (0.9 - 0.55) / 2. + 0.55 
   XYOUTS, 0.05, y, '[nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.4, ORIENTATION=90, CHARTHICK=1.5                       
;###############################################################################
; saving png
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
                write_jpeg, path+'diono_final_V5_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'diono_final_V5_'+Date+'.png', Image, R_bak, G_bak, B_bak
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'diono_final_V5_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	
end	
