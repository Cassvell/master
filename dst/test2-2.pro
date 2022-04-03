;file = DIALOG_PICKFILE(/READ, GET_PATH=dir)
file = DIALOG_PICKFILE(FILTER='*.dat')
file_name = FILE_BASENAME(file)

readcol, file, DATE, TIME, DOY, Dst, FORMAT = 'A,A,A,F', SKIPLINE = 25

DATETIME = DATE + ' ' + TIME

year = strmid(DATE, 0,4)
month = strmid(DATE, 5,2)
day = strmid(DATE, 8,2)


hour = strmid(TIME, 10, 2)
minute = strmid(TIME, 14, 2)
second = strmid(TIME, 17, 6)


i = where(Dst LE -150, count)
T = n_elements(DATETIME)

DUM = LABEL_DATE( DATE_FORMAT= ['%Y/%N/%D'])

tiempo = TIMEGEN(T, START=julday(month[0], day[0], year[0]), UNITS='Hours')

plot, tiempo[2863:3271], Dst[2863:3271],  XTICKFORMAT='LABEL_DATE', $
	title = 'tormenta geomagnetica intensa', $
	ytitle = 'Indice DST [nT]', $
	xstyle = 1


;PRINT, day[7000]

;pro graph

;for i=0, count-1 do begin
;	t = tiempo[i-72:i+240]
;	dst_i = Dst[i-72:i+240]
;	plot, t, dst_i, XTICKFORMAT='LABEL_DATE', $
;	title = 'TGM intensa', $
;	ytitle = 'Dst index [nT]', $
;	xstyle = 1
;endfor

;end


Date = strmid(DATETIME(0),0,4)

SET_PLOT, 'ps'
path = '../dst/output1'
psname = path+Date+'GMS_S.ps'
page_height = 27.94
page_width = 21.59
plot_left = 5.
plot_bottom = 5
xsize = 14
ysize = 10









