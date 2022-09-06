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


PRO km_plot, date_i, date_f
	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]      
;###############################################################################    
    file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
    tot_days= findgen(file_number*8)/8.0    
    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
        data_path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        
        string_date        = strarr(file_number)
               
        data_file_name_km  = strarr(file_number)
        data_file_name_kp  = strarr(file_number)  
                       
        string_date_2    = strarr(file_number)
             
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date_2[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')
                
		        data_file_name_kp[i] = data_path+'/kp/daily/kp_'+string_date_2[i]+'.txt' 	
                data_file_name_km[i] = data_path+'/Kmex/teo_'+string_date[i]+'.index.final'
                		       
		        file = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name_km[i] = data_path+'/Kmex/teo_'+string_date[i]+'.index.early'    
	            ENDIF 	
        	                            
        ENDFOR

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        exist_data_file_kp   = FILE_TEST(data_file_name_kp)
        capable_to_plot_kp   = N_ELEMENTS(where(exist_data_file_kp EQ 1))

        IF capable_to_plot_kp NE N_ELEMENTS(data_file_name_kp) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.kp_index.',A,' impossible to plot all data.')"              
        ENDIF

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))
        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.km_index.',A,' impossible to plot all data.')"              
        ENDIF
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
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]/10.
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF             
        ENDFOR
    k_days = findgen(file_number*8)/8.0 
    
        for i=0, n_elements(k_mex)-1 do begin
            if k_mex[i] eq 99.9 then begin
                k_mex[where(k_mex[*] eq 99.9)] = !Values.F_NAN          
            endif
        endfor    
;###############################################################################
;Kp Data                       
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

    IF indexes_03[0] NE -1 AND indexes_07[0] NE -1 THEN BEGIN
        kp[indexes_07] = kp[indexes_07]-0.3
        kp[indexes_03] = kp[indexes_03]+0.3
    ENDIF                  


    print, FORMAT='(6X,"kp",10X,"k_mex")'
    print, MAX(kp), MAX(k_mex)

 ;   WINDOW, 0,  XSIZE=800, YSIZE=800, TITLE='kmex vx kp'

 ;   PLOT, tot_days, kp, YRANGE=[0.0,9.0]
  ;  OPLOT, tot_days, k_mex, LINESTYLE=3





END















