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


FUNCTION ip_data, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]	

        header = 60      ; Defining number of lines of the header 
;###############################################################################
;reading data files
        date = string(year, month, format = '(I4, "-", I02)')
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'				
		file_name = path+'/ip/'+date+'.dat'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
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


PRO diono, r_ip, date_i, date_f
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
    tot_days= findgen(file_number*24)/24.0    
    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
        data_path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        
        string_date        = strarr(file_number)
        
        data_file_name_dh  = strarr(file_number)                
        data_file_name_km  = strarr(file_number)
        data_file_name_kp  = strarr(file_number)  
        data_file_name_dst = strarr(file_number) 
        data_file_name_tec  = strarr(file_number) 
                       
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
                data_file_name_dst[i]= data_path+'/dst/daily/dst_'+string_date_2[i]+'.txt'		        	        
                data_file_name_dh[i] = data_path+'/dH_teo/teo_'+string_date[i]+'.dst.early'
                data_file_name_tec[i]= data_path+'/tec/tec_'+string_date_2[i]+'.txt'
                data_file_name_km[i] = data_path+'/Kmex/teo_'+string_date[i]+'.index.final'
                		       
		        file = FILE_SEARCH(data_file_name_km[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name_km[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF 	
	            
		        file_dh = FILE_SEARCH(data_file_name_dh[i], COUNT=opened_files)
		        file_tec = FILE_SEARCH(data_file_name_tec[i], COUNT=opened_files)		        
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'    
	            ENDIF
        	                            
        ENDFOR


        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))

        exist_data_file_dst   = FILE_TEST(data_file_name_dst)
        capable_to_plot_dst   = N_ELEMENTS(where(exist_data_file_dst EQ 1))

        IF capable_to_plot_dst NE N_ELEMENTS(data_file_name_dst) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dst_index.',A,' impossible to plot all data.')"              
        ENDIF

        exist_data_file_kp   = FILE_TEST(data_file_name_kp)
        capable_to_plot_kp   = N_ELEMENTS(where(exist_data_file_kp EQ 1))

        IF capable_to_plot_kp NE N_ELEMENTS(data_file_name_kp) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.kp_index.',A,' impossible to plot all data.')"              
        ENDIF
                
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF
        
        exist_data_file_tec   = FILE_TEST(data_file_name_tec)
        capable_to_plot_tec   = N_ELEMENTS(where(exist_data_file_tec EQ 1))
        IF capable_to_plot_tec NE N_ELEMENTS(data_file_name_tec) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF

        exist_data_file_km   = FILE_TEST(data_file_name_km)
        capable_to_plot_km   = N_ELEMENTS(where(exist_data_file_km EQ 1))
        IF capable_to_plot_km NE N_ELEMENTS(data_file_name_km) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.km_index.',A,' impossible to plot all data.')"              
        ENDIF
;###############################################################################
; Generate the time variables to plot dH time series                                   
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
;###############################################################################
; Generate the time variables to plot TEC time series          
        tec  = fltarr(file_number*12)
        med  = fltarr(file_number*12)
        FOR i = 0, N_ELEMENTS(exist_data_file_tec)-1 DO BEGIN
                IF exist_data_file_tec[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date_2[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,X,I02,X,I02)'
                 
                        d_tec = tec_data([tmp_year, tmp_month, tmp_day])
                        
                        tec[i*12:(i+1)*12-1] = d_tec.tec[*]
                        med[i*12:(i+1)*12-1] = d_tec.med[*] 
                ENDIF ELSE BEGIN
                        tec[i*12:(i+1)*12-1] = 999.0
                        med[i*12:(i+1)*12-1] = 999.0
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
                        
                        k_mex[i*8:(i+1)*8-1] = d_km.k_mex[*]/10.
                        a_mex[i*8:(i+1)*8-1] = d_km.a_mex[*]
                ENDIF             
        ENDFOR
    k_days = findgen(file_number*8)/8.0 
;###############################################################################
;Dst Data                       
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
;###############################################################################        
;identifying NAN values in the Time Series
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
    tec_days= findgen(file_number*12)/12.0  
    tec_diff = tec-med    
;###############################################################################
; define ip parameters  
    ip = ip_data([yr_i, mh_i])
    
    ip_year = ip.year
    ip_doy  = ip.DOY
    ip_hour = ip.hour
    ip_Ey   = ip.E
    ip_flow = ip.flow_p

    tmp_doy = dst_doy[(idoy*24)-24:fdoy*24-1]
   
    Ey      = fltarr(n_elements(ip_doy))
    p_dyn   = fltarr(n_elements(ip_doy))
        
    tmp_doy = tmp_doy[uniq(tmp_doy, sort(tmp_doy))]

    for i=0, n_elements(tmp_doy)-1 do begin
            for j=0, n_elements(ip_doy)-1 do begin
                if tmp_doy[i] eq ip_doy[j] then begin
               p_dyn[j]    = ip_flow[j]
               Ey[j]       = ip_Ey[j]                            
                endif 
            endfor
    endfor
    
e_z = where(Ey eq 0, zcount, complement=val, ncomplement=valcount)
p_z = where(p_dyn eq 0, zcount2, complement=val_p, ncomplement=valcount_p)

Ey      = Ey[val]
p_dyn   = p_dyn[val]

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
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        Xsize=fix(800)
        Ysize=1000
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=0     
        DEVICE, set_character_size = [10, 12]
        DEVICE, DECOMPOSED=1     
        chr_size1 = 1.5
        chr_thick1= 1.5
        space     = 0.015
        rojo      = 248
        amarillo  = 220
        verde     = 150
        negro     = 0
        azul      = 70
        blanco    = 255
        gris      = 130
        morado    = 16
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT
;###############################################################################
; Time label    
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
;###############################################################################      
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
       days = intarr(file_number+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       days = days*24/24. 
       day_time = findgen(24)
;###############################################################################            
    if max(dst) gt max(H) then up = max(dst) else up = max(H)
    if min(dst) lt min(H) then down = min(dst) else down = min(H)
    
    window_title = 'TGM'+ string(TGM_n, format='(I01)')+', '+ $
                string(old_month, yr_i, format='(A, X, I4)')
;###############################################################################
; Plot data
    up_E     = max(Ey)
    down_E   = min(Ey)
    
    up_p    = max(p_dyn)
    down_p  = min(p_dyn) 
        
    plot, tot_days, Ey, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = 0.9, CHARTHICK=chr_thick1, POSITION=[0.1,0.76,0.9,0.93], $
     XSTYLE = 5, XRANGE=[0, file_number],  XTICKNAME=REPLICATE(' ', file_number+1), ySTYLE = 6,$
     YRANGE=[down_E,up_E], THICK=4  

z   	
    plot, tot_days, p_dyn, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=azul,$
     CHARSIZE = 0.9, CHARTHICK=chr_thick1, POSITION=[0.1,0.76,0.9,0.93], $
     XSTYLE = 5, XRANGE=[0, file_number], YRANGE=[down_p,up_p], $
     XTICKNAME=REPLICATE(' ', file_number+1), ySTYLE = 6, /noerase, THICK=4 ;, SUBTITLE = time_title 

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE='Tiempo Universal [dias]', $                           
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $      
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $  
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_E, up_E], $
                         ystyle=2, $  
                         YTITLE = 'Ey [mV/m]', $                          
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down_p,up_p], $
                         ystyle=2, $  
                         YTITLE = 'P [nPa]', $                          
                         COLOR=azul, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=3 
;###############################################################################
    up      = MAX(H+10)
    down    = MIN(H-10)    
    plot, tot_days, H, XTICKS=file_number, xminor = 8, POSITION=[0.1,0.52,0.9,0.69],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=5, YRANGE=[down, up],$
    title = '', BACKGROUND = blanco, COLOR=negro, XRANGE=[0, file_number],$
	ytitle = '',  XTICKNAME=REPLICATE(' ', file_number+1), /noerase, THICK=4

    oplot, tot_days, dst, COLOR=verde, linestyle=0, THICK=4   
     dH = TeXtoIDL('\DeltaH') 
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTITLE='Tiempo Universal [dias]', $  
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $      
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $  
                         CHARSIZE = 0.9, $                         
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down,up], $
                         YTITLE = dH+' y Dst [nT]', $
                         ystyle=1, $                          
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $       
                         COLOR=negro, $
                         ystyle=1, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5
;###############################################################################   
    plot, k_days, k_mex, XTICKS=file_number, xminor = 8, POSITION=[0.1,0.28,0.9,0.45],$
    XTICKFORMAT='LABEL_DATE', xstyle = 5, ystyle=5, YRANGE=[0,9],$
    BACKGROUND = blanco, COLOR=negro, XRANGE=[0, file_number],$
	ytitle = '',  XTICKNAME=REPLICATE(' ', file_number+1), /noerase, THICK=4

     OPLOT, k_days, kp, COLOR=verde, LINESTYLE=0, THICK=4

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9 , $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04
                         
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
                         CHARSIZE = 0.9 , $                       
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTITLE = 'Kp y Kmex', $                          
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         YSTYLE=1, $
                         CHARSIZE = 0.9 ;, $
                        
        AXIS, YAXIS = 1, YRANGE=[0,9], $        
                         COLOR=negro, $                         
                         YSTYLE=1, $
                         CHARTHICK=1.5,$
                         CHARSIZE = 0.9 ;, $  
;###############################################################################
    IF MAX(tec) GT MAX(med) THEN up_tecdiff = MAX(tec) ELSE up_tecdiff = MAX(med)
    IF MIN(tec) LT MIN(med) THEN down_tecdiff = MIN(tec) ELSE down_tecdiff = MIN(med)
        
    plot, tec_days, tec, XTICKS=file_number, xminor=8, BACKGROUND = blanco, $
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.1,0.04,0.9,0.21], $
     XSTYLE = 5, XRANGE=[0, file_number], XTICKNAME=REPLICATE(' ', file_number+1), ySTYLE = 6,$
     /noerase, YRANGE=[down_tecdiff, up_tecdiff], THICK=4, COLOR=rojo

        oplot, tec_days, med, color=negro, linestyle=0, THICK=4
    
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTITLE='Tiempo Universal [dias]', $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $  
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5,$                         
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down_tecdiff, up_tecdiff], $
                         ystyle=2, $  
                         YTITLE = 'TEC y <TEC> [TECu]', $                          
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         CHARSIZE = 0.9;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down_tecdiff, up_tecdiff], $
                         ystyle=2, $                    
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=1.5   
;###############################################################################
;second panel legend                   
        POLYFILL, [0.77,0.8,0.8,0.77], [0.564,0.564,0.567,0.567], color = negro, /NORMAL
        POLYFILL, [0.77,0.8,0.8,0.77], [0.534,0.534,0.537,0.537], color = verde, /NORMAL                
                
        XYOUTS, 0.81, 0.56 , /NORMAL, $
                dH, COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
                
        XYOUTS, 0.81, 0.53 , /NORMAL, $
                'Dst', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1     

;third panel legend  
        POLYFILL, [0.77,0.8,0.8,0.77], [0.424,0.424,0.427,0.427], color = negro, /NORMAL
        POLYFILL, [0.77,0.8,0.8,0.77], [0.394,0.394,0.397,0.397], color = verde, /NORMAL  

                XYOUTS, 0.81, 0.42 , /NORMAL, $
                'Kmex', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
                
        XYOUTS, 0.81, 0.39 , /NORMAL, $
                'Kp', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1  
                
;fourth panel legend  
        POLYFILL, [0.77,0.8,0.8,0.77], [0.184,0.184,0.187,0.187], color = negro, /NORMAL
        POLYFILL, [0.77,0.8,0.8,0.77], [0.154,0.154,0.157,0.157], color = rojo, /NORMAL  

                XYOUTS, 0.81, 0.18 , /NORMAL, $
                '<TEC>', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
                
        XYOUTS, 0.81, 0.15 , /NORMAL, $
                'TEC', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
;###############################################################################
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.97
   XYOuts, X, y, window_title, /Normal, color=negro, Alignment=0.5,$
   Charsize=2, CHARTHICK = 1.5  
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.949, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.8, CHARTHICK = 1.5   

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.71, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.8, CHARTHICK = 1.5     
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.468, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.8, CHARTHICK = 1.5     
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.231, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.8, CHARTHICK = 1.5              
;###############################################################################
   XYOuts, 0.14, 0.9, '(a)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.2, CHARTHICK= 3    
      
   XYOuts, 0.14, 0.55, '(b)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.2, CHARTHICK= 3 
   
   XYOuts, 0.14, 0.42, '(c)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.2, CHARTHICK= 3    
      
   XYOuts, 0.14, 0.17, '(d)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.2, CHARTHICK= 3    
;############################################################################### 
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak   
;###############################################################################
; open the post stript device
    path = '../rutidl/output/eventos_tgm/'
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'mag_tec_V2_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'mag_tec_V2_'+Date+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'mag_tec_V2_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	                    
END
