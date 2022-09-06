;+
; NAME:
;       scatter_plot.pro
;
;
; PURPOSE:
;
;       statistic data analysis
;
; AUTHOR:
;
;       Pedro Corona Romero
;       C. Isaac Castellanos Velazco
;       Instituto de Geofisica, UNAM
;       UNAM Campus Morelia, Antigua Carretera a Patzcuaron,
;       Exhacienda de San Jose del Cerrito, Morelia, Michoacan, Mexico
;       piter.cr@gmail.com
;       ccastellanos@igeofisica.unam.mx
;       06.iv.mmxxii
;
; CATEGORY:
;
;       Numerical Data Analize
;
; CALLING SEQUENCE:
;       .r scatter_plot
;       scatter_plot
;
;       Description:
;       ???????????????????????
;
;
; PARAMETERS:
;       Dst index    : dH index
;
; KEYWORD PARAMETERS:
;
;       /????????? : ?????????????
;
; DEPENDENCIAS:
;       ?????????? : ????????????
;
; ARCHIVOS ANALIZADOS:
;       ??????????
;
; ARCHIVOS DE SALIDA: grafica.png de dispersión
;
; HISTORIA:
;   versión 1.1, agosto 2022






PRO scatter_plot;, idate, fdate

        input_dir  = '/home/c-isaac/Escritorio/proyecto/master_thesis/rutidl/output/'
        path       = input_dir

	;iyr	= idate[0]
	;imh	= idate[1]
	;idy = idate[2]
	
	;fyr	= fdate[0]
	;fmh	= fdate[1]
	;fdy = fdate[2]					

       ; header = 1      ; Defining number of lines of the header 

       ; idate = string(iyr, imh, idy, format = '(I4, I02, I02)')
      ;  fdate = string(fyr, fmh, fdy, format = '(I4, I02, I02)')
                
        
        data_files         = FILE_SEARCH(input_dir+'tgmdata'+'?????????????????'+'.txt')
         
        files_lines_number = FILE_LINES(data_files)
        MX_latitude        = 28.06*!Pi/180.
       ; corellation_index  = FLTARR(N_ELEMENTS(data_files));<<<<<<<<<<<<<<<<preguntar
        
        limits = [-100, -400]
        
;print, data_files

        dat_str0 = { doy: INTARR(MAX(files_lines_number)), $
                    hour: FLTARR(MAX(files_lines_number)), $
                     dst: INTARR(MAX(files_lines_number)), $
                     dh: FLTARR(MAX(files_lines_number)), $
                     number_of_lines: 0 , $
                     correlation: FLTARR(N_ELEMENTS(limits)), $
                     correlation_points : INTARR(N_ELEMENTS(limits))}
        
        data = REPLICATE(dat_str0, N_ELEMENTS(data_files))
        data[*].number_of_lines = files_lines_number
;print, data[*].number_of_lines
        data[*].number_of_lines = files_lines_number-1
        data[*].correlation[*]  = 999.
;print, data[*].number_of_lines
;Return
         dh = FLTARR(1)
         
         dst = FLTARR(1)
        FOR event=0, N_ELEMENTS(data_files)-1 DO BEGIN
                data_strings  = STRARR(files_lines_number[event])
                
                OPENR, file_lun, data_files[event], /GET_LUN
                READF, file_lun, data_strings;, FORMAT = '(A)'
                FREE_LUN, file_lun
;print, data_strings
;return
                dat_str1 = { doy:0, hour: 0., dst: 0, dh: 0 }
                tmp_data = REPLICATE(dat_str1, data[event].number_of_lines)
                
                READS, data_strings[1:*], tmp_data, $
                FORMAT='(I03,X,F4.1,X,I04,F6.1)'
                
                data[event].doy[*] = 999
;print, data[event].number_of_lines, N_ELEMENTS(tmp_data[*].doy)
;RETURN
                data[event].doy[0:data[event].number_of_lines-1]=tmp_data[*].doy
                
                data[event].hour[*] = 999.
                data[event].hour[0:data[event].number_of_lines-1]=tmp_data[*].hour
                
                data[event].dst[*] = 999.
                data[event].dst[0:data[event].number_of_lines-1]=tmp_data[*].dst*cos(MX_latitude)
                
                data[event].dh[*] = 999.
                data[event].dh[0:data[event].number_of_lines-1]=tmp_data[*].dh
                              
                dh  = [dh, tmp_data[*].dh]
                dst = [dst, tmp_data[*].dst*cos(MX_latitude)]
               ; print, correlate(dh, dst)^2
                tmp_index = WHERE(data[event].dh[*] GE 999.)
                data[event].dst[tmp_index] = 999.

        ENDFOR

      ;  print, CORRELATE(dh, dst)^2
        
        idx = WHERE(dst GE -100)      
        j =   WHERE(dst GE -400 AND dst LT -100)
        
        corr1 = CORRELATE(dh[idx], dst[idx])^2
        corr2 =  CORRELATE(dh[j], dst[j])^2
        
      ;  print, corr1
        ;print, N_ELEMENTS(dh[j])
        
    correlation = FLTARR(N_ELEMENTS(data_files))
        FOR event=0, N_ELEMENTS(data_files)-1 DO BEGIN
            correlation[event] = CORRELATE(data[event].dst[0:data[event].number_of_lines-1], $
            data[event].dh[0:data[event].number_of_lines-1])
            
         ;   print, N_ELEMENTS(data[*].dst[0:data[event].number_of_lines-1])
                FOR i=0, N_ELEMENTS(limits)-1 DO BEGIN
                        IF i EQ 0 THEN index_tmp = WHERE(data[event].dst[*] GE limits[i] AND data[event].dst[*] LT 200 ) $
                                  ELSE index_tmp = WHERE(data[event].dst[*] GE limits[i] AND data[event].dst[*] LT limits[i-1] )
      
                        IF N_ELEMENTS(index_tmp) GT 2 THEN data[event].correlation[i] = CORRELATE(data[event].dst[index_tmp], data[event].dh[index_tmp])
                        data[event].correlation_points[i] = N_ELEMENTS(index_tmp)
                ENDFOR
        ENDFOR
        
   ; corr_med    = MEAN(correlation)^2
   ; corr_std    = STDDEV(correlation)^2
   ; print,  corr_med,  corr_std      
        means=FLTARR(N_ELEMENTS(limits))
        devs =FLTARR(N_ELEMENTS(limits))
        maxs =FLTARR(N_ELEMENTS(limits))
        mins =FLTARR(N_ELEMENTS(limits))
        ;print, means
        for i=0, N_ELEMENTS(limits)-1 DO BEGIN
                ;print, n_elements(limits), limits
            ;    print, data[*].correlation[i]
                index_tmp = WHERE(data[*].correlation[i] LT 100)
                print, index_tmp
                ;print, data[index_tmp].correlation[i]
                means[i] = Mean( data[index_tmp].correlation[i]^2 )
                devs[i]  = STDDEV( data[index_tmp].correlation[i]^2, /NAN)
                maxs[i]  = MAX( data[index_tmp].correlation[i] )
                mins[i]  = MIN( data[index_tmp].correlation[i] )
                ;print, means[i], devs[i];, maxs[i], mins[i]
              ;  print, data[index_tmp].correlation[i]^2;, data[index_tmp].doy
        ENDFOR

;###############################################################################          
; define device and color parameters 
;###############################################################################      
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        Xsize=fix(600)
        Ysize=600
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
             
        chr_size1 = 0.9
        chr_thick1= 1.1
        space     = 0.015
        rojo      = 248
        amarillo  = 200
        naranja   = 220
        verde     = 100
        negro     = 0
        azul      = 90
        blanco    = 255
        ;gris      = 130
        morado    = 16
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 0, /SILENT

    ;path = '../rutidl/output/eventos_tgm/'
    H = TeXtoIDL('\DeltaH')
        event = 0
   up0  =  100
   down0=-500        
   slope = findgen(601)-500
   x1= slope
   y1= slope
   plot, x1,y1, xstyle=5, ystyle=5, color=negro, background=blanco, YRANGE=[down0,up0], $
   XRANGE=[down0,up0], CHARSIZE = 1.2, CHARTHICK=1.2, Xtitle='Dst', Ytitle='H',$
   POSITION=[0.17,0.1,0.9,0.9], TITLE='Dispersion entre Dst y '+H, /NODATA                 
    

        POLYFILL, [10.,-210.,-210.,10.], [0.,0.,0.5,0.5], color=amarillo
    
  ;  POLYFILL, [!X.CRANGE[0], 0, 0, !X.CRANGE[0]], [!Y.CRANGE[0], !Y.CRANGE[0], 0, 0], COLOR=verde
    
    POLYFILL, [!X.CRANGE[0], -100, -100, !X.CRANGE[0]], [!Y.CRANGE[0], !Y.CRANGE[0], -100, -100], COLOR=amarillo            

        AXIS, Xaxis=0, Xtitle='Dst [nT]', CHARSIZE=1.1, color=negro, XRANGE=[down0,up0], Xstyle=1, $
        CHARTHICK=2
        AXIS, Xaxis=1, color=negro, XTICKNAME=REPLICATE( ' ', 8 ), XRANGE=[down0,up0], $
        CHARTHICK=2
        
        AXIS, Yaxis=0, Ytitle=H+' [nT]', CHARSIZE=1.1, color=negro, Ystyle=1, $
        CHARTHICK=2
        AXIS, yaxis=1, color=negro, YTICKNAME=REPLICATE( ' ', 7 ), $
        CHARTHICK=2


    OPLOT, dst, dh, PSYM=4, COLOR=azul, SYMSIZE=0.5        
        
    OPLOT, x1, y1, THICK=1.0, COLOR=negro 
;###############################################################################                             
;###############################################################################
    med = texToidl('R^2 = ')
        XYOUTS, 0.2, 0.75 , /NORMAL, $
             ;   med, COLOR=negro, $
                string(med, 0.76, FORMAT='(A, F4.2)') , COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
                
    desv = texToidl('(100 \leq Dst \leq -100)')                
        XYOUTS, 0.2, 0.80 , /NORMAL, $
                 desv, COLOR=negro, $
                ;desv+' '+corr_std, COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1 
;###############################################################################
    med = texToidl('R^2 = ')
        XYOUTS, 0.2, 0.55 , /NORMAL, $
             ;   med, COLOR=negro, $
                string(med, 0.63, FORMAT='(A, F4.2)') , COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1   
                
    desv = texToidl('(-100 \leq Dst \leq -400)')                
        XYOUTS, 0.2, 0.60 , /NORMAL, $
                desv , COLOR=negro, $
                ;desv+' '+corr_std, COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=chr_thick1 

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
                write_jpeg, path+'dispersion_plot.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'dispersion_plot.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'dispersion_plot.png'
                print, ''
        ENDIF
        RETURN 	


RETURN
END




