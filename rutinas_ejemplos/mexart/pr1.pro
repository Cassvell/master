dataStruct = { dt:0.0, x:0.0}

file = 'M20190509c48'
nrows = File_Lines(file)
data = Replicate(dataStruct, nrows)

OpenR, lun, file, /GET_LUN
ReadF, lun, data
Free_Lun, lun

;print, data.dt
;print, data.x
;plot, data.x

;dt_t = ((data.dt, 2), 0, 10) )
;fx = Fix( StrMid( StrTrim(data.x, 2), 4, 2) )

y = fft(data.x)
print, y
plot, y
