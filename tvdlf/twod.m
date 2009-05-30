rd	3	#
		da $1
		lines 1 1
		read {_nx 1 _ny 2 _Lx 3 _Ly 4 _t 5}
 		define $2 (_Lx/2)
 		define $3 (_Ly/2)
 		define nx (_nx)
 		define ny (_ny)
 		define dx (1./_nx)
 		define dy (1./_ny)
		echo $nx $ny
		lines 2 10000000
		read {x 1 y 2 rho 3 ux 4 uy 5 uz 6 Bx 7 By 8 Bz 9 u 10}
		#
		set gam = 5./3.
                set p = u*(gam - 1.)
                set lp = lg(p) 
pl	2	#
		rd $1 Lx Ly
		define ma 1.1
		if($Lx > $Ly) {define L ($ma*$Lx)} else {define L ($ma*$Ly)}
		set i=1,$nx*$ny
                set y=(i - $ny*int((i-0.5)/$ny) - 0.5)*$dy
		set x=(int((i-0.5)/$ny) + 0.5)*$dx
                #set x=(i - $nx*int((i-0.5)/$nx) - 0.5)*$dx
		#set y=(int((i-0.5)/$nx) + 0.5)*$dy
		image($nx,$ny) -$Lx $Lx -$Ly $Ly
		set ix = int(x/$dx)
		set iy = int(y/$dy)
		set image[ix,iy] = $2[i-1]
		#
		limits -$L $L -$L $L
		erase
		box
		minmax min max echo $min $max
		if($min*$max < 0.) {\
			define delta (($max-$min)/10.)
			set lev=$min,-$delta,$delta
			levels lev
			ltype 2
			contour
			#
			set lev=$delta,$max,$delta
			levels lev
			ltype 0
			contour
		} \
		else {\
			set lev=$min,$max,($max-$min)/10.
			levels lev
			ltype 0
			contour
		}
		#
		#
		#
vpl     3       #
                rd $1 Lx Ly
                define ma 1.1
                if($Lx > $Ly) {define L ($ma*$Lx)} else {define L ($ma*$Ly)}
                set i=1,$nx*$ny
                define dx (2.*$Lx/$nx)
                define dy (2.*$Ly/$ny)
                set y=(i - $ny*int((i-0.5)/$ny) - 0.5)*$dy - $Ly
                set x=(int((i-0.5)/$ny) + 0.5)*$dx - $Lx
                #
                limits -$L $L -$L $L
                erase
                box
                #
                set VVx = $2x
                set VVy = $2y
                set ang=180.*atan2(VVy,VVx)/pi
                set len=$3*sqrt(VVx**2 + VVy**2)
                vfield x y len ang
