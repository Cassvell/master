;Name:
;	indexdata_V2
;purpose:
;	leer datos de índices geomagnéticos y de contenido total de electrones para graficarlos 
;   
;author:
;	Carlos Isaac Castellanos Velazco
;	Estudiante de Maestría en Ciencias de la Tierra
;	Instituto de Geofísica, Unidad Michoacan
;	UNAM
;	ccastellanos@igeofisica.unam.mx
;
;category:
;   plotting
;
;calling sequence:
;   .r indexdata_V2.pro
;   indexdata_V2, initialdate[yyyy,mm,dd], finaldate[yyyy,mm,dd]
;parameters:
;
;
;dependencies:
;
;
;input files
;   kp data, dst data, kmex data, dh data, TEC data
;
;output files:
;   figura.

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


		r_kp = REPLICATE(DataStruct, number_of_lines-header)	
                
		READS, data[header:number_of_lines-1], r_kp, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,X,I2,I04,I02,A1,I04)'
			
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
		path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
		file_name = path+'/Kmex/'+name+'final'		
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name =  path+'/Kmex/'+name+'early'
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

  
pro indexdata_v2, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

;##############################################################################
; defining time window
;############################################################################## 
	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]

    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
    file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1    
    days_dst       = FINDGEN(file_number*24)/24.
;###############################################################################   
    tec_days= FINDGEN(file_number*12)/12.0  
;###############################################################################          
    k_days= FINDGEN(file_number*8)/8.0       
;###############################################################################
        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        string_date        = strarr(file_number)
                
        data_file_name_dh  = strarr(file_number)
        data_file_name_km  = strarr(file_number) 
        data_file_name_dst = strarr(file_number)
        data_file_name_kp  = strarr(file_number)  
        data_file_name_tec  = strarr(file_number)        
        string_date2       = strarr(file_number)            
;###############################################################################
; Generate the time variables to plot time series of DH Index    
    data_path = '/home/c-isaac/Escritorio/proyecto/master_thesis/datos'   
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date2[i]= string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')
                                
                data_file_name_dh[i] = data_path+'/dH_teo/'+'teo_'+string_date[i]+'.dst.early'

                data_file_name_dst[i]= data_path+'/dst/daily/dst_'+string_date2[i]+'.txt'
                
                data_file_name_km[i] = data_path+'/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        data_file_name_kp[i]= data_path+'/kp/daily/kp_'+string_date2[i]+'.txt'  
                data_file_name_tec[i]= data_path+'/tec/'+'tec_'+string_date2[i]+'.txt'	
                
		        file_km = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)		                
	            IF opened_files NE N_ELEMENTS(file_km) THEN begin
	                data_file_name_km[i] = data_path+'/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF    
	            
		        file_dh = FILE_SEARCH(data_file_name_dh[i], COUNT=opened_files)		                
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_dh[i] = data_path+'/dH_teo/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF  	                                                           
        ENDFOR

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
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
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.km_index.',A,' impossible to plot all data.')"              
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
                string_date2[i] = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')                
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
                string_date2[i] = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')                
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

    kp[indexes_07] = kp[indexes_07]-0.3
    kp[indexes_03] = kp[indexes_03]+0.3    
;###############################################################################
;TEC DATA        
        exist_data_file_tec   = FILE_TEST(data_file_name_tec)
        capable_to_plot_tec   = N_ELEMENTS(where(exist_data_file_tec EQ 1))        
        
        IF capable_to_plot_tec NE N_ELEMENTS(data_file_name_tec) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF 
           
        tec  = fltarr(file_number*12)
        med  = fltarr(file_number*12)
        FOR i = 0, N_ELEMENTS(exist_data_file_tec)-1 DO BEGIN
                IF exist_data_file_tec[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date2[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,X,I02,X,I02)'                 
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
;###############################################################################
; initiate the figure Device, defining colors and figure dim
;###############################################################################
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        IF file_number GT 7 THEN tmp_spam = 1.5
        IF file_number GT 15 THEN tmp_spam = 2.
        
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
;###############################################################################
; Create a time label based on the DOY initial and DOY final inputs
;###############################################################################
     X_label   = STRARR(file_number+1)+' '
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
; plot the Dst time series for each event
;###############################################################################

window_title = 'from '+string(yr_i, mh_i, dy_i, $
                FORMAT='(I4, "/", I02, "/", I02)')+' to '$
                +string(yr_f, mh_f, dy_f, $
                FORMAT='(I4, "/", I02, "/", I02)')

time_title = 'Time [UT]' 
                
MAG_source = 'Source: International Service of Geomagnetic Indices'  

    if max(H) gt max(dst) then up0 = max(H) else up0 = max(dst)
    if min(H) lt min(dst) then down0 = min(H) else down0 = min(dst)
            
    plot, days_dst, H, XTICKS=file_number, xminor = 8, POSITION=[0.07,0.66,0.95,0.96],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down0,up0],$
    title = window_title, BACKGROUND = blanco, COLOR=negro, XRANGE=[0, file_number],$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', file_number+1)

    oplot, days_dst, dst, COLOR=azul
  
    if up0-down0 gt 300 then begin
        for i = -600., 100., 100. do oplot, [0,file_number], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -600., 100., 50. do oplot, [0,file_number], [i,i], linestyle=1, COLOR=negro
    endelse
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.6, $
;                         CHARTHICK=chr_thick1, $
                         XTICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
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

;###############################################################################     
    plot, k_days, Kp, psym=6, /NoDATA, MAX_VALUE=9., XTICKS=file_number, xminor=8, $
                    TITLE = plot_title, SUBTITLE = '', XTITLE = 's', YTITLE = 's', $
                    BACKGROUND = blanco, COLOR=negro, YRANGE=[0,9], YTICKS=9, $
                    YMINOR=0, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                    POSITION=[0.07,0.36,0.95,0.64], XSTYLE = 5, ystyle = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], /NOERASE

        j = N_ELEMENTS(Kp)
        FOR i = 0, j-1 DO BEGIN
                IF Kp[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (Kp[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                Kp[i] EQ 4 : color = amarillo
                                                Kp[i] GE 4 : color = rojo
                                                ELSE       : color = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], [0,0,Kp[i]$
                                        +step,Kp[i]+step], color=color
                                      ENDIF $
                                      ELSE BEGIN                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., linestyle = 0                                         
                                         
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+k_days[i], $
                                        [0,0,9.,9.],color=morado, /LINE_FILL, $
                                         ORIENTATION=-45., linestyle = 0
                                      ENDELSE
        ENDFOR

 FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,90.], linestyle=1, COLOR=negro
        ENDFOR

        FOR i=5, 8, 1 DO OPLOT, [0,file_number], [i,i], linestyle=1, COLOR=negro

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.7, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
                         COLOR=negro, $
                         TICKLEN=0.04


        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTITLE = 'Kp index', $
                         COLOR=negro, $
                         CHARSIZE = 0.7, $
                         CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $

        XYOUTS, 0.01, 0.070 , /NORMAL, $
                MAG_source, COLOR=negro, $
                CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1

;###############################################################################
    up      = max(tec)
    down    = min(tec)
    
    plot, tec_days, tec, XTICKS=file_number, xminor = 8, POSITION=[0.07,0.1,0.95,0.34],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=6, YRANGE=[down, up],$
    BACKGROUND = blanco, COLOR=negro, XRANGE=[0, file_number], $
	XTICKNAME=REPLICATE(' ', file_number+1), /noerase

    oplot, tec_days, med, COLOR=rojo
    
    if up-down gt 120 then begin
        for i = -100., 260., 50. do oplot, [0,file_number], [i,i], linestyle=1, COLOR=negro
    endif else begin
        for i = -100., 260., 10. do oplot, [0,file_number], [i,i], linestyle=1, COLOR=negro
    endelse 
         
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         ;XTICKFORMAT='(A1)',$
                         xtitle='time [UT]', $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
;                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKFORMAT='(A1)',$                         
                         ;XTICKNAME=REPLICATE(' ', file_number+1), $
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
;###############################################################################
    if file_number gt 15 then begin
    
        XYOUTS, 0.01, 0.165 , /NORMAL, $
                'Color Code:        quiet,            disturbed,           storm,                   data not available.', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1
       
        POLYFILL, [0.07,0.10,0.10,0.07], [0.165,0.165,0.195,0.195], color = verde, /NORMAL
        POLYFILL, [0.15,0.18,0.18,0.15], [0.165,0.165,0.195,0.195], color = amarillo, /NORMAL
        POLYFILL, [0.25,0.28,0.28,0.25], [0.165,0.165,0.195,0.195], color = rojo, /NORMAL
        POLYFILL, [0.37,0.40,0.40,0.37], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., linestyle = 0
        POLYFILL, [0.37,0.40,0.40,0.37], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., linestyle = 0
    
    endif else begin
        XYOUTS, 0.01, 0.165 , /NORMAL, $
                'Color Code:       quiet,          disturbed,           storm,                 data not available.', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1
      
        POLYFILL, [0.09,0.12,0.12,0.09], [0.165,0.165,0.195,0.195], color = verde, /NORMAL
        POLYFILL, [0.19,0.22,0.22,0.19], [0.165,0.165,0.195,0.195], color = amarillo, /NORMAL
        POLYFILL, [0.32,0.35,0.35,0.32], [0.165,0.165,0.195,0.195], color = rojo, /NORMAL
        POLYFILL, [0.47,0.50,0.50,0.47], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., linestyle = 0
        POLYFILL, [0.47,0.50,0.50,0.47], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., linestyle = 0
    endelse
;first panel legend
        POLYFILL, [0.77,0.80,0.80,0.77], [0.674,0.674,0.676,0.676], color = negro, /NORMAL
        POLYFILL, [0.86,0.89,0.89,0.86], [0.674,0.674,0.676,0.676], color = azul, /NORMAL        
    if file_number gt 7 then begin
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
    if file_number gt 7 then begin
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
    if file_number gt 7 then begin
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
