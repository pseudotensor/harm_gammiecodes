pl	2	#
		rd $1 Lx Ly
		echo $Lx $Ly
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
cuta	0	#
		set x1 = -.5,.5,.05
		set y1 = image(x1,0.)
		limits x1 y1 
		erase
		box
		connect x1 y1
		#
cutb	2	#
		rd $1
		set i=1,$nx*$ny
                set x=(i - $nx*int((i-0.5)/$nx) - 0.5)*$dx
                set y=(int((i-0.5)/$nx) + 0.5)*$dy
                image($nx,$ny) -.5 .5 -.5 .5
                set ix = int(x/$dx)
                set iy = int(y/$dy)
                set image[ix,iy] = $2[i-1]
		set x1 = -.5,.5,.05
		set y1 = image(x1,0.)
		connect x1 y1
vpl	3	#
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
vplp	4	#
		rd $1 Lx Ly
		set i=1,$nx*$ny
		set ix = int((i-0.5)/$ny)
		set iy = i - $ny*ix
		#set iy = int(i/$nx) + 1
		#set ix = i - $nx*(iy - 1)
		#
		set use = (int(iy/$4) - iy/$4 == 0 && \
			int(ix/$4) - ix/$4 == 0) ? 1 : 0
		#
		define dx (2.*$Lx/$nx)
		define dy (2.*$Ly/$ny)
		#
                set y=((i - $ny*int((i-0.5)/$ny) - 0.5)*$dy - $Ly) if (use)
		set x=((int((i-0.5)/$ny) + 0.5)*$dx - $Lx) if (use)
		#
		set VVx = $2x if(use)
		set VVy = $2y if(use)
		set ang=180.*atan2(VVy,VVx)/pi
		set len=$3*sqrt(VVx**2 + VVy**2)
		vfield x y len ang
rd	3	#
		da $1
		lines 1 1
		read {_nx 1 _ny 2 _Lx 3 _Ly 4 _W 5 _q 6 _t 7}
		define $2 (_Lx/2)
		define $3 (_Ly/2)
		define nx (_nx)
		define ny (_ny)
		define dx (_Lx/_nx)
		define dy (_Ly/_ny)
		lines 2 1000000
		#
		read {X 1 Y 2 r 3 p 4 vx 5 vy 6 pot 7 pv 8 dv 9}
		set g = 3. - 2./1.4
		set ek = 0.5*r*(vx*vx + vy*vy)
		set dvy = vy + 1.5*1*X
		set dvx = vx
		set dpy = r*dvy
		set dpx = r*dvx
		set T = p/r
		set px = r*vx
		set py = r*vy
		set eth = p/(g - 1)
		set ent = p/r**(1.7)
		set enth = g*eth/r
		set Be = 0.5*(vx*vx + vy*vy) + pot + enth
		set pv = 2./r + pv
		set cv = -dv
		set dr = r - 1. 
		#
		##read {X 1 Y 2 r 3}
plx	2	#
		rd $1
		limits X $2
		erase
		box
		ptype 4 0
		points X $2
		connect X $2
		xla X
		yla $2
ply	2	#
		rd $1
		limits Y $2
		erase
		box
		ptype 4 0
		points Y $2
		connect Y $2
		xla Y
		yla $2
pla	2	#
		rd $1 Lx Ly
		set i=1,$nx*$ny
                set y=(i - $ny*int((i-0.5)/$ny) - 0.5)*$dy
		set x=(int((i-0.5)/$ny) + 0.5)*$dx
		image($nx,$ny) -$Lx $Lx -$Ly $Ly
		set ix = int(x/$dx)
		set iy = int(y/$dy)
		set image[ix,iy] = $2[i-1]
		#
		minmax min max echo $min $max
		#
		set lev=$min,$max,($max-$min)/10.
		levels lev contour
		contour
		#
anim	1	#
		rd dump000 Lx Ly
		limits X $1
		define nl ($fy1 - 0.2*($fy2 - $fy1))
		define nh ($fy2 + 0.8*($fy2 - $fy1))
		limits X $nl $nh
		#limits X 0 8
		box
		#
		pls dump001 $1
		pls dump002 $1
		pls dump003 $1
		pls dump004 $1
		pls dump005 $1
		pls dump006 $1
		pls dump007 $1
		pls dump008 $1
		pls dump009 $1
		pls dump010 $1
		pls dump011 $1
		pls dump012 $1
		pls dump013 $1
		pls dump014 $1
		pls dump015 $1
		pls dump016 $1
		pls dump017 $1
		pls dump018 $1
		pls dump019 $1
		pls dump020 $1
pls	2	#
		rd $1 Lx Ly
		connect X $2
		#
plb	2	#
		rd $1 Lx Ly
		define ma 1.1
		if($Lx > $Ly) {define L ($ma*$Lx)} else {define L ($ma*$Ly)}
		set i=1,$nx*$ny
                set y=(i - $ny*int((i-0.5)/$ny) - 0.5)*$dy
		set x=(int((i-0.5)/$ny) + 0.5)*$dx
		image($nx,$ny) -$Lx $Lx -$Ly $Ly
		set ix = int(x/$dx)
		set iy = int(y/$dy)
		set image[ix,iy] = $2[i-1]
		#
		erase
		box
		contour
		#
labe	1	#
		limits 0 1 0 1
		relocate 0.5 1.02
		putlabel 5 $1
