file = DIALOG_PICKFILE(/READ, GET_PATH=dir)
file_name = FILE_BASENAME(file)

readcol, file, DATE, TIME, DOY, Dst, FORMAT = 'A,A,A,F,F', SKIPLINE = 25

DATETIME = DATE+' '+TIME

year = (strmid(DATETIME, 0,4))
month = (strmid(DATETIME, 5,2))
day = (strmid(DATETIME, 8,2))
hour = (strmid(DATETIME, 10, 2))
minute = (strmid(DATETIME, 14, 2))
second = (strmid(DATETIME, 17, 6))

Date = strmid(DATETIME(0),0,4)
OutFile = '../dst/output1/'+Date+'GSM2.txt'
i = where(Dst LE -150, count)

;PRINT, Dst[i], N_ELEMENTS(i)
;OPENW, 1, Outfile
;PRINTF,1, DATETIME[i], Dst[i], i,  format = "(A23,'   ' ,A24, '   ', A24)"
;close,1

PRINT, Date, Outfile





	
