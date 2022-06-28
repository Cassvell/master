file = DIALOG_PICKFILE(/READ, GET_PATH=dir)
file_name = FILE_BASENAME(file)

readcol, file, dt, fx

n_int=12
T = 0.02
secs=10.0
n=LONG(secs/T)
interval=n_int*n

plot, dt, fx, title=file_name
print, "Click cursor at ON-SOURCE start time"
	cursor, x, y, /down
i0 = 0L

while (dt[i0] Lt x) do i0=i0+1
t1=dt[i0]
t2=dt[i0+interval-1]

dt1 = dt[i0:i0+interval-1]
fx1 = fx[i0:i0+interval-1]
;print, dt[i0:i0+interval-1]
;help, dt[i0:i0+interval-1]


y = fft(fx1)

plot, y
;print, y


f_k = findgen(n_elements(fx1))/(n_elements(fx1)*dt1)
;print, f_k

psw_s = abs(y)^2
psw_s = smooth(psw_s, 25)



plot, f_k, psw_s, /xlog, /ylog, title = 'power spectrum vs freq', xrange=[1e-4,1],yrange=[1e-10,1e-4], psym=3
;print, f_k, psw_s

END

